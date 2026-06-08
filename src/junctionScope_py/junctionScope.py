#!/usr/bin/env python3
"""
junctionScope.py — splice-junction calling pipeline (single-cell aware).

Supports two config modes:
  1. Single junction:  JXNCOORD + GENE (+ optional NTSEQ) in config
  2. Junction table:   JXN_TABLE pointing to a TSV with columns:
					   junction_name, gene, junction, nt_seq

Per-sample output files are named <sample>.<junction_name>.*
Output columns (one line per sample × junction):
  sample, gene, target_gene_coord,
  targetJxn_coord, samtools_jxn_coord, regtools_jxn_coord, nt_sequence,
  targetJxn_read_count, targetJxn_seq_read_count, pct_reads_seq,
  targetJxn_cell_count, targetJxn_seq_cell_count, pct_cells_seq,
  geneJxn_count_mean, geneJxn_count_median, geneJxn_count_mode,
  geneJxn_count_sum, geneJxn_counts, geneJxn_transcripts

Dependencies:
	pysam, pandas  (pip install pysam pandas)
	regtools       (system / module load — used for junctions annotate only)

Usage:
	python junctionScope.py -c project.conf [--setup-only]
	python junctionScope.py -c project.conf --sample SampleA --bam /path/to/A.cram
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO
from collections import defaultdict
import shutil




# =============================================================================
# Config / junction table parsing
# =============================================================================

def parse_config(config_path: str) -> dict:
	"""Parse a shell-style KEY="value" config file into a dict."""
	config = {}
	with open(config_path) as fh:
		for line in fh:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			line = re.sub(r'\s*#.*$', '', line)
			m = re.match(r'^(\w+)\s*=\s*"?([^"]*)"?\s*$', line)
			if m:
				config[m.group(1)] = m.group(2)
	config.setdefault("BUFFER", "200")
	config.setdefault("THREADS", "1")
	config.setdefault("NTSEQ", "NULL")
	config.setdefault("GENE", "NULL")
	config.setdefault("JXNCOORD", "NULL")
	config.setdefault("JXN_TABLE", "NULL")
	config.setdefault("REGTOOLS", False)
	return config


def validate_config(config: dict):
	"""Require INPUT/FASTA/GTF/OUTPUT, and either JXNCOORD+GENE or JXN_TABLE."""
	required_always = ["INPUT", "FASTA", "GTF", "OUTPUT"]
	missing = [k for k in required_always if not config.get(k)]
	if missing:
		sys.exit(f"[junctionScope] Missing required config keys: {', '.join(missing)}")
	has_single = config["JXNCOORD"] != "NULL" and config["GENE"] != "NULL"
	has_table  = config["JXN_TABLE"] != "NULL"
	if not has_single and not has_table:
		sys.exit(
			"[junctionScope] Config must supply either:\n"
			"  JXNCOORD + GENE   (single-junction mode)\n"
			"  JXN_TABLE         (multi-junction table mode)"
		)


def validate_regtools(verbose: bool = False) -> bool:
	"""
	Validate that regtools is available and that the
	junctions extract and annotate commands exist.
	Parameters
	----------
	verbose : bool
		Print diagnostic messages.
	Returns
	-------
	bool
		True if regtools appears usable, otherwise False.
	"""
	# Check executable exists
	regtools_path = shutil.which("regtools")
	if regtools_path is None:
		if verbose:
			print("regtools not found in PATH")
		return False
	try:
		# Check main executable works
		result = subprocess.run(
			["regtools"],
			capture_output=True,
			text=True,
			timeout=10
		)
		output = (
			result.stdout +
			result.stderr
		).lower()
		if "junctions" not in output:
			if verbose:
				print("regtools executable found but output unexpected")
			return False
		# Check extract command
		extract_result = subprocess.run(
			["regtools", "junctions", "extract"],
			capture_output=True,
			text=True,
			timeout=10
		)
		extract_output = (
			extract_result.stdout +
			extract_result.stderr
		).lower()
		# Check annotate command
		annotate_result = subprocess.run(
			["regtools", "junctions", "annotate"],
			capture_output=True,
			text=True,
			timeout=10
		)
		annotate_output = (
			annotate_result.stdout +
			annotate_result.stderr
		).lower()
		extract_ok = (
			"usage" in extract_output
			or "bam" in extract_output
		)
		annotate_ok = (
			"usage" in annotate_output
			or "gtf" in annotate_output
		)
		if not extract_ok:
			if verbose:
				print("regtools junctions extract unavailable")
		if not annotate_ok:
			if verbose:
				print("regtools junctions annotate unavailable")
		return extract_ok and annotate_ok
	except Exception as e:
		if verbose:
			print(f"regtools validation failed: {e}")
		return False

def parse_jxn_table(path: str) -> list[dict]:
	"""
	Load a junction table TSV.  Accepted column names (case-insensitive):
		junction_name / name
		gene
		junction / jxncoord / coord
		nt_seq / ntseq / sequence
	Returns a list of dicts with normalised keys:
		junction_name, gene, junction, nt_seq
	"""
	df = pd.read_csv(path, sep="\t", comment="#")
	df.columns = [c.strip().lower() for c in df.columns]
	col_map = {
		"junction_name": ["junction_name", "name", "jxn_name"],
		"gene":          ["gene", "gene_symbol"],
		"junction":      ["junction", "jxncoord", "coord", "junction_coord"],
		"nt_seq":        ["nt_seq", "ntseq", "sequence", "seq"],
	}
	rename = {}
	for canonical, aliases in col_map.items():
		for alias in aliases:
			if alias in df.columns and canonical not in df.columns:
				rename[alias] = canonical
				break
	df = df.rename(columns=rename)
	for col in ["junction_name", "gene", "junction"]:
		if col not in df.columns:
			sys.exit(f"[parse_jxn_table] Required column '{col}' not found in {path}")
	if "nt_seq" not in df.columns:
		df["nt_seq"] = "NULL"
	df["nt_seq"] = df["nt_seq"].fillna("NULL").astype(str)
	return df[["junction_name", "gene", "junction", "nt_seq"]].to_dict("records")


def config_to_jxn_list(config: dict) -> list[dict]:
	"""
	Return a normalised list of junction dicts regardless of config mode.
	Each dict has: junction_name, gene, junction, nt_seq
	"""
	if config["JXN_TABLE"] != "NULL":
		return parse_jxn_table(config["JXN_TABLE"])
	else:
		return [{
			"junction_name": config.get("JXN_NAME", config["GENE"]),
			"gene":          config["GENE"],
			"junction":      config["JXNCOORD"],
			"nt_seq":        config["NTSEQ"],
		}]


# =============================================================================
# Region helpers
# =============================================================================

def parse_region(region: str) -> tuple:
	"""'chr1:1000-2000' → ('chr1', 1000, 2000)"""
	chrom, coords = region.split(":")
	start, end = coords.split("-")
	return chrom, int(start), int(end)


def add_buffer(region: str, buffer: int) -> str:
	"""Return a buffered region string chr:start-end."""
	if buffer <= 0:
		return region
	chrom, start, end = parse_region(region)
	return f"{chrom}:{max(0, start - buffer)}-{end + buffer}"


# =============================================================================
# NT flanking sequence  (pysam.FastaFile replaces samtools faidx)
# =============================================================================

def get_nt_seq(region: str, fasta: str, window: int = 5, rev_comp: bool = False) -> str:
	"""
	Extract <window> bases flanking each side of a junction and concatenate.
	pysam.FastaFile uses 0-based half-open coordinates.
	"""
	chrom, start, end = parse_region(region)
	# Convert from 1-based closed (shell script convention) to 0-based half-open
	left_start  = start - window - 1
	left_end    = start - 1
	right_start = end
	right_end   = end + window

	with pysam.FastaFile(fasta) as fa:
		left_seq  = fa.fetch(chrom, left_start, left_end)
		right_seq = fa.fetch(chrom, right_start, right_end)

	seq = (left_seq + right_seq).upper()
	if rev_comp:
		seq = seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
	return seq

def load_fasta_kmers(
	fasta_path: str,
	revC: bool = True,
	k: int = 31
):
	reference_sequences = {}
	kmer_index = defaultdict(set)
	for record in SeqIO.parse(fasta_path, "fasta"):
		seq = str(record.seq).upper()
		reference_sequences[record.id] = seq
		# forward kmers
		for i in range(len(seq) - k + 1):
			kmer = seq[i:i+k]
			kmer_index[kmer].add(record.id)
		# reverse complement
		if revC:
			rc_name = record.id + "_revC"
			rc_seq = str(
				record.seq.reverse_complement()
			).upper()
			reference_sequences[rc_name] = rc_seq
			for i in range(len(rc_seq) - k + 1):
				kmer = rc_seq[i:i+k]
				kmer_index[kmer].add(rc_name)
	return reference_sequences, kmer_index


# =============================================================================
# GTF helpers
# =============================================================================

def subset_gtf(gtf_file: str, gene: str, out_gtf: str):
	"""Write a GTF subset containing only rows that mention the gene symbol."""
	with open(gtf_file) as fh, open(out_gtf, "w") as out:
		for line in fh:
			if f'"{gene}"' in line:
				out.write(line)


def gene_region_from_gtf(gtf_file: str, chrom: str) -> str:
	"""Return chr:min_start-max_end from all features in a (subsetted) GTF."""
	starts, ends = [], []
	with open(gtf_file) as fh:
		for line in fh:
			if line.startswith("#"):
				continue
			fields = line.split("\t")
			if len(fields) >= 5:
				try:
					starts.append(int(fields[3]))
					ends.append(int(fields[4]))
				except ValueError:
					pass
	if not starts:
		return None
	return f"{chrom}:{min(starts)}-{max(ends)}"


# =============================================================================
# Region extraction  (pysam replaces samtools view)
# =============================================================================

def extract_jxn_region(
	input_file: str,
	region: str,
	output_sam: str,
	buffer: int = 200,
	threads: int = 1,
):
	"""
	Subset a BAM/CRAM to a buffered region and write a SAM file.
	Auto-detects BAM vs CRAM from file suffix.
	"""
	suff = Path(input_file).suffix.lower()
	mode = "rb" if suff == ".bam" else "rc"
	buffered = add_buffer(region, buffer)
	with pysam.AlignmentFile(input_file, mode, threads=threads) as bam, \
		 pysam.AlignmentFile(output_sam, "w", header=bam.header) as out:
		for read in bam.fetch(region=buffered):
			out.write(read)


# =============================================================================
# SAM parser  (produces one row per junction per read, with barcode)
# =============================================================================

def parse_sam_per_junction(sam_file: str) -> pd.DataFrame:
	"""
	Parse a SAM file into a per-junction DataFrame.

	One row per (read × junction).  Reads with no N in CIGAR get one row
	with junction_start/end = None.  Columns include cell_barcode (CB:Z: tag).
	"""
	records = []
	with open(sam_file) as f:
		for line in f:
			if line.startswith("@"):
				continue
			fields = line.rstrip().split("\t")
			if len(fields) < 11:
				continue
			read_id    = fields[0]
			chrom      = fields[2]
			read_start = int(fields[3])   # 1-based SAM POS
			cigar      = fields[5]
			seq        = fields[9]
			read_len   = len(seq)
			# Extract cell barcode from auxiliary tags
			barcode = None
			umi = None
			for tag in fields[11:]:
				if tag.startswith("CB:Z:"):
					barcode = tag[5:]
				if tag.startswith("UB:Z:"):
					umi = tag[5:]
			# Parse CIGAR for N operations
			cigar_ops = [
				(int(n), op)
				for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
			]
			ref_pos   = read_start
			ref_len   = 0
			junctions = []
			for length, op in cigar_ops:
				if op in {"M", "D", "N", "=", "X"}:
					ref_len += length
				if op == "N":
					junction_start = ref_pos
					junction_end = ref_pos + length - 1
					junctions.append(
						(junction_start, junction_end)
					)
				# move genomic position
				if op in {"M", "D", "N", "=", "X"}:
					ref_pos += length
			query_length = ref_len
			query_end    = read_start + ref_len - 1
			n_junctions  = len(junctions)
			base_row = {
				"read_id":       read_id,
				"chromosome":    chrom,
				"read_start":    read_start,
				"reference_end": query_end,
				"query_length":  ref_len,
				"read_length":   read_len,
				"seq":           seq,
				"cell_barcode":  barcode,
				"umi":           umi,
				"barcode_umi":   (str(barcode)+'_'+str(umi)),
				"n_junctions":   n_junctions,
				"multi_junction":n_junctions > 1,
			}
			if n_junctions == 0:
				records.append({**base_row,
					"junction_start": None, "junction_end": None, "junction_index": None})
			else:
				for idx, (j_start, j_end) in enumerate(junctions, start=1):
					records.append({**base_row,
						"junction_start": j_start,
						"junction_end":   j_end,
						"junction_index": idx})
	df = pd.DataFrame(records)
	# Downcast float64 nullable int columns (junction coords can be None)
	for col in df.select_dtypes("float64").columns:
		df[col] = df[col].astype("Int64")
	return df


# =============================================================================
# Filtering helpers
# =============================================================================

def filter_to_junction(
	parsed_sam: pd.DataFrame,
	coord: str,
	buffer: int = 2,
	filter_barcodes: bool = False,
	chr_col: str = 'chromosome',
	start_col: str = 'junction_start',
	end_col: str = 'junction_end',
	barcode_col: str = 'cell_barcode'
	):
	"""Filter parsed sam file to junction coordinates of interest."""
	chrom, start, end = parse_region(coord)
	parsed_sam_coord = parsed_sam[
								(parsed_sam[chr_col] == chrom) &
								(parsed_sam[start_col].between(start-buffer,start+buffer)) &
								(parsed_sam[end_col].between(end-buffer,end+buffer))]
	if filter_barcodes:
		return pd.DataFrame(parsed_sam_coord)
	all_barcodes = (parsed_sam[[barcode_col]].drop_duplicates())
	matched_barcodes = (parsed_sam_coord[[barcode_col]].drop_duplicates())
	missing_barcodes = (all_barcodes[~all_barcodes[barcode_col].isin(matched_barcodes[barcode_col])].copy())
	if not missing_barcodes.empty:
		missing_barcodes[chr_col] = chrom
		missing_barcodes[start_col] = start
		missing_barcodes[end_col] = end
		# optional placeholder cols
		extra_cols = [
			c for c in parsed_sam.columns
			if c not in missing_barcodes.columns
			]
		for c in extra_cols:
			missing_barcodes[c] = pd.NA
		missing_barcodes = missing_barcodes[parsed_sam_coord.columns]
		parsed_sam_coord = pd.concat([parsed_sam_coord,missing_barcodes],ignore_index=True)
	return pd.DataFrame(parsed_sam_coord)


def seq_check(
	parsed_sam: pd.DataFrame,
	seq: str,
	seq_col: str = "seq",
	flag_col: str = "junction_seq_found",
	filter_rows: bool = False,
	case_sensitive: bool = False
):
	parsed_sam = parsed_sam.copy()
	if case_sensitive:
		parsed_sam[flag_col] = parsed_sam[seq_col].str.contains(
			seq,
			regex=False,
			na=False
		)
	else:
		parsed_sam[flag_col] = (
			parsed_sam[seq_col]
			.str.upper()
			.str.contains(seq.upper(), regex=False, na=False)
		)
	if filter_rows:
		parsed_sam = parsed_sam[parsed_sam[flag_col]]
	return parsed_sam


def query_bam_to_fasta(
	input_file: str,
	fasta_path: str,
	revC: bool = True,
	threads: int = 8,
	k: int = 31
):
	total_reads = 0
	matched_reads = 0
	ref_fa, kmer_index = load_fasta_kmers(
		fasta_path,
		revC,
		k
	)
	suff = Path(input_file).suffix
	if suff == ".bam":
		mode = "rb"
	elif suff == ".cram":
		mode = "rc"
	else:
		raise ValueError("Input must be BAM or CRAM")
	bam = pysam.AlignmentFile(
		input_file,
		mode,
		threads=threads
	)
	matches = []
	for read in bam.fetch(until_eof=True):
		total_reads += 1
		if read.query_sequence is None:
			continue
		query = read.query_sequence.upper()
		# skip short reads
		if len(query) < k:
			continue
		# choose seed kmer
		seed = query[:k]
		# fast lookup
		candidate_refs = kmer_index.get(seed, set())
		# no possible match
		if not candidate_refs:
			continue
		bc = read.get_tag("CB") if read.has_tag("CB") else None
		umi = read.get_tag("UB") if read.has_tag("UB") else None
		# full confirmation
		for ref_name in candidate_refs:
			ref_seq = ref_fa[ref_name]
			if query in ref_seq:
				matched_reads += 1
				strand = (
					"-"
					if ref_name.endswith("_revC")
					else "+"
				)
				clean_ref = re.sub(
					r"_revC$",
					"",
					ref_name
				)
				matches.append({
					"read_name": read.query_name,
					"barcode": bc,
					"umi": umi,
					"reference": clean_ref,
					"strand": strand,
					"sequence": query
				})
	bam.close()
	matches_df = pd.DataFrame(matches)
	print(f"{matched_reads}/{total_reads} matched")
	return matches_df


# =============================================================================
# Junction summarisation
# =============================================================================

def summ_junc(
	parsed_sam: pd.DataFrame,
	):
	parsed_sam_summ = (
		parsed_sam
		.dropna(subset=["junction_start", "junction_end"])
		.groupby(
			["chromosome", "junction_start", "junction_end"]
		)
		.agg(
			# total junction support
			sam_read_count=("read_id", "nunique"),
			sam_n_cells=("cell_barcode", "nunique"),
			sam_n_umi=("barcode_umi", "nunique"),
			multi_junction_reads=("multi_junction", "sum"),
			# sequence-supported reads
			seq_read_count=(
				"read_id",
				lambda x: (parsed_sam.loc[x.index]
					.query("junction_seq_found == True")["read_id"]
					.nunique()
					)
			),
			# sequence-supported cells
			seq_n_cells=(
				"cell_barcode",
				lambda x: (parsed_sam.loc[x.index]
					.query("junction_seq_found == True")["cell_barcode"]
					.nunique()
					)
			),
			seq_n_umi=(
				"barcode_umi",
				lambda x: (parsed_sam.loc[x.index]
					.query("junction_seq_found == True")["barcode_umi"]
					.nunique()
				)
			),
			multi_junction_reads_uniq=(
				"read_id",
				lambda x: (parsed_sam.loc[x.index]
					.query("multi_junction == True")["read_id"]
					.nunique()
				)
			)
		)
		.reset_index()
	)
	parsed_sam_summ["pct_reads_seq"] = (
		parsed_sam_summ["seq_read_count"]
		/ parsed_sam_summ["sam_read_count"]
	)
	parsed_sam_summ["pct_cells_seq"] = (
		parsed_sam_summ["seq_n_cells"]
		/ parsed_sam_summ["sam_n_cells"]
	)
	parsed_sam_summ["pct_umi_seq"] = (
		parsed_sam_summ["seq_n_umi"]
		/ parsed_sam_summ["sam_n_umi"]
	)
	return parsed_sam_summ


def summ_junc_bybcd(parsed_sam: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize junction support by chromosome/start/end/barcode.

    Treat <NA> values in boolean columns as not counted, without modifying
    the original parsed_sam table.
    """
    def n_unique_where_true(group_idx, value_col, flag_col):
        mask = parsed_sam.loc[group_idx, flag_col].fillna(False).astype(bool)
        return parsed_sam.loc[group_idx, value_col][mask].nunique()
    parsed_sam_summ = (
        parsed_sam
        .groupby(
            ["chromosome", "junction_start", "junction_end", "cell_barcode"],
            dropna=False
        )
        .agg(
            reads_per_barcode=("read_id", "nunique"),
            seq_reads_per_barcode=(
                "read_id",
                lambda x: n_unique_where_true(
                    x.index,
                    value_col="read_id",
                    flag_col="junction_seq_found"
                )
            ),
            umis_per_barcode=("barcode_umi", "nunique"),
            seq_umis_per_barcode=(
                "barcode_umi",
                lambda x: n_unique_where_true(
                    x.index,
                    value_col="barcode_umi",
                    flag_col="junction_seq_found"
                )
            ),
            multi_junction_reads=(
                "multi_junction",
                lambda x: x.fillna(False).astype(bool).sum()
            ),
            multi_junction_reads_uniq=(
                "read_id",
                lambda x: n_unique_where_true(
                    x.index,
                    value_col="read_id",
                    flag_col="multi_junction"
                )
            )
        )
        .reset_index()
    )
    parsed_sam_summ["pct_reads_seq_bcd"] = (
        parsed_sam_summ["seq_reads_per_barcode"]
        / parsed_sam_summ["reads_per_barcode"]
    )
    parsed_sam_summ["pct_umis_seq_bcd"] = (
        parsed_sam_summ["seq_umis_per_barcode"]
        / parsed_sam_summ["umis_per_barcode"]
    )
    return parsed_sam_summ

def select_best_junction(
	junct_summ: pd.DataFrame,
	min_reads: int = 5,
	min_seq_reads: int = 1
	) -> pd.DataFrame:
	"""Select highest-confidence junction from summary table."""
	df = junct_summ.copy()
	# optional pre-filter
	df = df[
		(df["sam_read_count"] >= min_reads) &
		(df["seq_read_count"] >= min_seq_reads)
	]
	if df.empty:
		raise ValueError(
			"No junctions pass filtering criteria."
		)
	df = df.sort_values(
		by=[
			"seq_read_count",
			"pct_reads_seq",
			"sam_read_count",
			"sam_n_cells",
			"sam_n_umi"
		],
		ascending=False
	)
	return df.head(1)


# =============================================================================
# Total junction count across whole BAM  (pysam replaces samtools view | wc)
# =============================================================================

def count_total_junctions(input_file: str, threads: int = 1) -> dict:
	"""
	Count total splice junctions and spliced reads across a whole BAM/CRAM.
	Returns {'read_count': n_spliced_reads, 'junction_count': n_junctions}.
	"""
	mode = "rb" if Path(input_file).suffix.lower() == ".bam" else "rc"
	j_count = r_count = 0
	with pysam.AlignmentFile(input_file, mode, threads=threads) as bam:
		for read in bam.fetch(until_eof=True):
			if read.is_secondary or read.is_supplementary or not read.cigartuples:
				continue
			n = sum(1 for op, _ in read.cigartuples if op == 3)
			if n:
				r_count += 1
				j_count += n
	return {"read_count": r_count, "junction_count": j_count}


# =============================================================================
# regtools extract + annotate  (subprocess — no pure-Python equivalent for annotate)
# =============================================================================

def extract_jxn_regtools(
	input_bam: str,
	region: str,
	output_bed: str,
	buffer: int = 1000,
):
	"""
	Run regtools junctions extract for all three strandedness modes and merge.
	Handles regtools ≥1.x strand argument change.
	"""
	buffered = add_buffer(region, buffer)
	try:
		ver = subprocess.run(
			["regtools", "--version"], capture_output=True, text=True
		)
		lines = ver.stderr.splitlines() + ver.stdout.splitlines()
		ver_str = next((l.split()[1] for l in lines if l.strip()), "0.0.0")
		major = int(ver_str.split(".")[0])
	except Exception:
		major = 0
	strand_args = ["XS", "RF", "FR"] if major >= 1 else ["0", "1", "2"]
	first = True
	with open(output_bed, "w") as out:
		for strand in strand_args:
			result = subprocess.run(
				["regtools", "junctions", "extract",
				 "-r", buffered, "-s", strand, input_bam],
				capture_output=True, text=True, check=True,
			)
			lines = result.stdout.splitlines(keepends=True)
			if first:
				out.writelines(lines)
				first = False
			else:
				for line in lines:
					if not line.startswith("chrom"):
						out.write(line)


def annotate_jxn_regtools(
	input_bed: str,
	fasta: str,
	gtf: str,
	output_anno: str,
):
	"""Run regtools junctions annotate and deduplicate output rows."""
	result = subprocess.run(
		["regtools", "junctions", "annotate", input_bed, fasta, gtf],
		capture_output=True, text=True, check=True,
	)
	lines = result.stdout.splitlines()
	if not lines:
		Path(output_anno).write_text("")
		return
	header     = lines[0]
	data_lines = sorted(set(lines[1:]))
	with open(output_anno, "w") as out:
		out.write(header + "\n")
		out.write("\n".join(data_lines))
		if data_lines:
			out.write("\n")


def read_anno_as_df(anno_file: str) -> pd.DataFrame:
	"""Read a regtools annotate output file into a DataFrame with tidy column names."""
	df = pd.read_csv(anno_file, sep="\t")
	# Normalise column names produced by different regtools versions
	rename = {
		"chrom": "chromosome", "start": "junction_start", "end": "junction_end",
		"score": "score_regtools",
		"known_donor": "known_donor", "known_acceptor": "known_acceptor",
		"known_junction": "known_junction",
		"genes": "genes", "transcripts": "transcripts",
	}
	df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
	return df


# =============================================================================
# Gene-level junction summary from regtools anno
# =============================================================================

def get_jxn_counts_gene(reg_anno: pd.DataFrame, gene: str) -> dict:
	"""
	Summarise per-junction regtools counts for known junctions within a gene.
	Filters to rows where known_donor == known_acceptor == known_junction == 1.
	"""
	mask = (
		(reg_anno.get("known_donor",    1) == 1) &
		(reg_anno.get("known_acceptor", 1) == 1) &
		(reg_anno.get("known_junction", 1) == 1) &
		(reg_anno["genes"] == gene)
	)
	sub = reg_anno[mask]
	scores = sub["score_regtools"] if "score_regtools" in sub.columns else sub["score"]
	mode_val = ",".join(scores.mode().astype(str)) if not scores.empty else "NA"
	transcripts = (
		",".join(sub["transcripts"].dropna().unique())
		if "transcripts" in sub.columns else "NA"
	)
	return {
		"mean_val":       round(scores.mean(),   4) if not scores.empty else 0,
		"median_val":     round(scores.median(), 4) if not scores.empty else 0,
		"mode_val":       mode_val,
		"sum_val":        int(scores.sum()),
		"counts":         ",".join(scores.astype(str)),
		"transcript_ids": transcripts,
	}


# =============================================================================
# SAM ↔ regtools coordinate matching
# =============================================================================

def match_sam_to_regtools(
	sam_summ: pd.DataFrame,
	reg_anno: pd.DataFrame,
	coord_tol: int = 1,
) -> pd.DataFrame:
	"""
	Merge the SAM junction summary with the regtools annotation on chromosome,
	compute coordinate distances, and flag matches within coord_tol bases.
	Returns the merged DataFrame sorted by coord_diff_distance.
	"""
	sam = sam_summ.copy()
	reg = reg_anno.copy()

	merged = sam.merge(reg, on="chromosome", suffixes=("_sam", "_reg"))
	merged["junction_start_diff"] = (
		(merged["junction_start_sam"] - merged["junction_start_reg"]).abs().astype(int)
	)
	merged["junction_end_diff"] = (
		(merged["junction_end_sam"] - merged["junction_end_reg"]).abs().astype(int)
	)
	merged["coord_diff_distance"] = (
		merged["junction_start_diff"] + merged["junction_end_diff"]
	)
	merged["junction_start_match"] = merged["junction_start_diff"] <= coord_tol
	merged["junction_end_match"]   = merged["junction_end_diff"]   <= coord_tol
	merged["junction_match"] = (
		merged["junction_start_match"] & merged["junction_end_match"]
	)
	# Put gene/transcript annotation columns at the end
	for col in ["genes", "transcripts"]:
		if col in merged.columns:
			merged = merged[[c for c in merged.columns if c != col] + [col]]

	return merged.sort_values("coord_diff_distance").reset_index(drop=True)


# =============================================================================
# Single-junction single-sample runner
# =============================================================================

OUTPUT_COLUMNS = [
	"sample", "gene", "jxn_name", "target_gene_coord",
	"targetJxn_coord", "samtools_jxn_coord", "regtools_jxn_coord", "nt_sequence",
	"targetJxn_read_count", "targetJxn_seq_read_count", "pct_reads_seq",
	"targetJxn_cell_count", "targetJxn_seq_cell_count", "pct_cells_seq",
	"targetJxn_umi_count", "targetJxn_seq_umi_count", "pct_umi_seq",
	"geneJxn_count_mean", "geneJxn_count_median", "geneJxn_count_mode",
	"geneJxn_count_sum", "geneJxn_counts", "geneJxn_transcripts",
]

OUTPUT_COLUMNS_noReg = [
	"sample", "gene", "jxn_name", "target_gene_coord",
	"targetJxn_coord", "samtools_jxn_coord", "regtools_jxn_coord", "nt_sequence",
	"targetJxn_read_count", "targetJxn_seq_read_count", "pct_reads_seq",
	"targetJxn_cell_count", "targetJxn_seq_cell_count", "pct_cells_seq",
	"targetJxn_umi_count", "targetJxn_seq_umi_count", "pct_umi_seq"
]



step_n = 1

def run_sample_junction(
	sample: str,
	bam_file: str,
	jxn_entry: dict,
	gene_region: str,
	fasta: str,
	gtf: str,
	output_dir: str,
	buffer: int = 200,
	threads: int = 1,
	regtools: bool = False
) -> dict:
	"""
	Run the full pipeline for one sample × one junction entry.

	Parameters
	----------
	sample      : sample name (used for file naming)
	bam_file    : path to BAM/CRAM
	jxn_entry   : dict with keys junction_name, gene, junction, nt_seq
	gene_region : chr:start-end for the full gene (for regtools)
	fasta       : reference FASTA
	gtf         : (subsetted) GTF for this gene
	output_dir  : Intermediate/<sample> directory
	buffer      : nt buffer for region extraction
	threads     : samtools threads
	regtools    : bool == True to run regtools steps

	Returns a dict matching OUTPUT_COLUMNS.
	"""
	jxn_name = jxn_entry["junction_name"]
	gene     = jxn_entry["gene"]
	jxncoord = jxn_entry["junction"]
	nt_seq   = jxn_entry["nt_seq"]
	step_count = 6 if regtools and validate_regtools() else 5
	global step_n
	step_n = 1

	# Derive NT sequence if not provided
	if nt_seq == "NULL" or not nt_seq:
		nt_seq = get_nt_seq(jxncoord, fasta, window=5)
		print(f"  [NT seq] derived: {nt_seq}")

	# Resolve output_dir to absolute so file I/O works regardless of cwd
	output_dir = str(Path(output_dir).resolve())
	os.makedirs(output_dir, exist_ok=True)
	# Base path prefix for all intermediate files for this sample × junction
	base = os.path.join(output_dir, f"{sample}.{jxn_name}")

	# ------------------------------------------------------------------
	# 1. Extract region → SAM
	# ------------------------------------------------------------------
	region_sam = f"{base}.region.sam"
	print(f"  [{step_n}/{step_count}] Extracting region {jxncoord} ...")
	step_n += 1
	extract_jxn_region(bam_file, jxncoord, region_sam, buffer=buffer, threads=threads)

	# ------------------------------------------------------------------
	# 2. Parse SAM to DataFrame
	# ------------------------------------------------------------------
	print(f"  [{step_n}/{step_count}] Parsing SAM ...")
	step_n += 1
	parsed = parse_sam_per_junction(region_sam)
	parsed.to_csv(f"{base}.region.parsed.tsv", sep="\t", index=False)

	# ------------------------------------------------------------------
	# 3. Filter to junction coordinates + sequence check
	# ------------------------------------------------------------------
	print(f"  [{step_n}/{step_count}] Filtering to junction coordinates ...")
	step_n += 1
	filtered = filter_to_junction(parsed, jxncoord)
	#filtered = filter_to_junction(parsed, jxncoord, filter_barcodes = True)
	filtered = seq_check(filtered, nt_seq)
	filtered.to_csv(f"{base}.region.filtered.tsv", sep="\t", index=False)

	# ------------------------------------------------------------------
	# 4. Summarise — sample level and per-barcode level
	# ------------------------------------------------------------------
	print(f"  [{step_n}/{step_count}] Summarising junctions ...")
	step_n += 1
	sam_summ    = summ_junc(filtered)
	bcd_summ    = summ_junc_bybcd(filtered)
	sam_summ.to_csv(f"{base}.sam.jxn.summ.tsv",     sep="\t", index=False)
	bcd_summ.to_csv(f"{base}.sam.jxn.bcd.summ.tsv", sep="\t", index=False)

	# Best junction from SAM method
	best = select_best_junction(sam_summ)
	if best.empty:
		sam_jxn_coord   = "NA"
		read_count      = 0
		seq_read_count  = 0
		pct_reads_seq   = 0.0
		cell_count      = 0
		seq_cell_count  = 0
		pct_cells_seq   = 0.0
		umi_count       = 0
		seq_umi_count   = 0
		pct_umi_seq     = 0.0
		bcd_summ_best   = pd.DataFrame()
	else:
		top = best.iloc[0]
		sam_jxn_coord  = f"{top['chromosome']}:{int(top['junction_start'])}-{int(top['junction_end'])}"
		read_count     = int(top["sam_read_count"])
		seq_read_count = int(top["seq_read_count"])
		pct_reads_seq  = round(float(top["pct_reads_seq"]), 4)
		cell_count     = int(top["sam_n_cells"])
		seq_cell_count = int(top["seq_n_cells"])
		pct_cells_seq  = round(float(top["pct_cells_seq"]), 4)
		umi_count      = int(top["sam_n_umi"])
		seq_umi_count  = int(top["seq_n_umi"])
		pct_umi_seq    = round(float(top["pct_umi_seq"]), 4)

		# Per-barcode subset for the best junction
		bcd_summ_best = bcd_summ[
			(bcd_summ["chromosome"]     == top["chromosome"]) &
			(bcd_summ["junction_start"] == top["junction_start"]) &
			(bcd_summ["junction_end"]   == top["junction_end"])
		].copy()
		bcd_summ_best.to_csv(f"{base}.sam.jxn.bcd.best.tsv", sep="\t", index=False)

	# ------------------------------------------------------------------
	# 5. regtools extract + annotate (gene-level)
	# ------------------------------------------------------------------
	if regtools and validate_regtools():
		print(f"  [{step_n}/{step_count}] Running regtools ...")
		step_n += 1
		reg_bed  = f"{base}.gene.bed"
		reg_anno = f"{base}.gene.bed.anno"
		extract_jxn_regtools(bam_file, gene_region, reg_bed, buffer=1000)
		annotate_jxn_regtools(reg_bed, fasta, gtf, reg_anno)

		reg_df = read_anno_as_df(reg_anno)

		# Match SAM summary to regtools annotation
		reg_jxn_coord = "NA"
		if not sam_summ.empty and not reg_df.empty:
			merged = match_sam_to_regtools(sam_summ, reg_df)
			merged.to_csv(f"{base}.sam_reg.merged.tsv", sep="\t", index=False)

			# Best regtools match: closest coordinate to the SAM best junction
			matched = merged[merged["junction_match"]]

		# Gene-level junction summary from regtools
		gene_summ = get_jxn_counts_gene(reg_df, gene)

	# ------------------------------------------------------------------
	# 6. Assemble result row
	# ------------------------------------------------------------------
	print(f"  [{step_n}/{step_count}] Writing result ...")
	step_n += 1
	base_result = {
		"sample":                   sample,
		"gene":                     gene,
		"jxn_name":                 jxn_name,
		"target_gene_coord":        gene_region,
		"targetJxn_coord":          jxncoord,
		"samtools_jxn_coord":       sam_jxn_coord,
		"nt_sequence":              nt_seq,
		"targetJxn_read_count":     read_count,
		"targetJxn_seq_read_count": seq_read_count,
		"pct_reads_seq":            pct_reads_seq,
		"targetJxn_cell_count":     cell_count,
		"targetJxn_seq_cell_count": seq_cell_count,
		"pct_cells_seq":            pct_cells_seq,
		"targetJxn_umi_count":      umi_count,
		"targetJxn_seq_umi_count":  seq_umi_count,
		"pct_umi_seq":              pct_umi_seq,
	}


	if regtools and validate_regtools():
		base_result.update({
		"geneJxn_count_mean":       gene_summ["mean_val"],
		"geneJxn_count_median":     gene_summ["median_val"],
		"geneJxn_count_mode":       gene_summ["mode_val"],
		"geneJxn_count_sum":        gene_summ["sum_val"],
		"geneJxn_counts":           gene_summ["counts"],
		"geneJxn_transcripts":      gene_summ["transcript_ids"]})
		# matched may be unbound if sam_summ or reg_df was empty
		if 'matched' not in dir():
			matched = pd.DataFrame()

		if matched.empty:
			# No coordinate match — one row with NA regtools coord
			results = [{**base_result, "regtools_jxn_coord": "NA"}]
		else:
			results = []
			for _, r in matched.iterrows():
				reg_coord = (
					f"{r['chromosome']}:"
					f"{int(r['junction_start_reg'])}-{int(r['junction_end_reg'])}"
				)
				# Pull in extra regtools columns not already in base_result
				extra = {
					col: r[col]
					for col in matched.columns
					if col not in base_result and col != "regtools_jxn_coord"
				}
				results.append({**base_result, "regtools_jxn_coord": reg_coord, **extra})
			# for _, r in matched.iterrows():
			# 	result = base_result.copy()
			# 	result["regtools_jxn_coord"] = (
			# 		f"{r['chromosome']}:"
			# 		f"{int(r['junction_start_reg'])}-{int(r['junction_end_reg'])}"
			# 	)
			# 	# Optional: keep extra matched/regtools columns too, if they exist
			# 	for col in matched.columns:
			# 		if col not in result:
			# 			result[col] = r[col]
			# 	results.append(result)
	else:
		# regtools unavailable — single row, no gene-level stats
		results = [{**base_result, "regtools_jxn_coord": "NA"}]

	result_df = pd.DataFrame(results)

	ordered_cols = [c for c in OUTPUT_COLUMNS if c in result_df.columns]
	extra_cols = [c for c in result_df.columns if c not in ordered_cols]

	#result_df[ordered_cols + extra_cols].to_csv(
	result_df[ordered_cols].to_csv(
		f"{base}.junctScope.txt",
		sep="\t",
		index=False
	)
	#result_line = "\t".join(str(result[c]) for c in OUTPUT_COLUMNS)
	#with open(f"{base}.junctScope.txt", "w") as out:
	#	out.write(result_line + "\n")
	return result_df


# =============================================================================
# Setup  — derive gene regions, write per-sample execution scripts
# =============================================================================

def setup(config: dict) -> tuple[list[dict], list[tuple[str, str]]]:
	"""
	Prepare project directory, subset GTFs, derive gene regions, derive NT
	sequences where missing, and write per-sample-per-junction execution scripts.

	Returns (jxn_list, samples) so callers can iterate directly.
	"""
	# Resolve every path to absolute so generated scripts and file I/O work
	# regardless of which directory the caller is in.
	output     = str(Path(config["OUTPUT"]).resolve())
	fasta      = str(Path(config["FASTA"]).resolve())
	gtf        = str(Path(config["GTF"]).resolve())
	input_list = str(Path(config["INPUT"]).resolve())
	buffer     = int(config["BUFFER"])
	threads    = int(config["THREADS"])
	regtools   = bool(config["REGTOOLS"])

	# Store resolved paths back so callers that access config directly get abs paths
	config["OUTPUT"]     = output
	config["FASTA"]      = fasta
	config["GTF"]        = gtf
	config["INPUT"]      = input_list

	os.makedirs(output, exist_ok=True)

	jxn_list = config_to_jxn_list(config)

	# ── Per-gene GTF subsets and gene regions ────────────────────────────────
	genes = list(dict.fromkeys(j["gene"] for j in jxn_list))   # unique, ordered
	gene_gtfs    = {}
	gene_regions = {}

	for gene in genes:
		gene_gtf_path = str(Path(output) / f"temp_{gene}.gtf")
		subset_gtf(gtf, gene, gene_gtf_path)
		gene_gtfs[gene] = gene_gtf_path

		chrom = jxn_list[0]["junction"].split(":")[0]   # assume all same chrom
		region = gene_region_from_gtf(gene_gtf_path, chrom)
		if not region:
			sys.exit(f"[setup] No GTF entries found for gene '{gene}'")
		gene_regions[gene] = region
		print(f"[setup] {gene} region: {region}")

	# ── Derive missing NT sequences ──────────────────────────────────────────
	for jxn in jxn_list:
		if jxn["nt_seq"] == "NULL" or not jxn["nt_seq"]:
			jxn["nt_seq"] = get_nt_seq(jxn["junction"], fasta, window=5)
			print(f"[setup] NT seq for {jxn['junction_name']}: {jxn['nt_seq']}")

	# Store derived values back on config for direct run_sample_junction calls
	config["_gene_gtfs"]    = gene_gtfs
	config["_gene_regions"] = gene_regions
	config["_jxn_list"]     = jxn_list

	# ── Parse sample list ────────────────────────────────────────────────────
	samples = []
	with open(input_list) as fh:
		for line in fh:
			line = line.strip().strip('"').replace("\r", "")
			if not line:
				continue
			parts = line.split("\t")
			if len(parts) == 1:
				bam_file = parts[0]
				sample   = Path(bam_file).stem
			else:
				sample   = parts[0].strip('"')
				bam_file = parts[1].strip('"')
			samples.append((sample.replace(".", ""), bam_file))

	# ── Per-sample execution scripts ─────────────────────────────────────────
	output_abs  = str(Path(output).resolve())   # already absolute; defensive copy
	master_exc  = str(Path(output_abs) / f"{Path(output_abs).name}_exc.sh")
	script_path = str(Path(__file__).resolve())

	with open(master_exc, "w") as master:
		for sample, bam_file in samples:
			# Resolve bam_file to absolute so generated scripts are location-independent
			bam_abs = str(Path(bam_file).resolve())
			if not os.path.isfile(bam_abs):
				print(f"[setup] WARNING: {bam_abs} not found — skipping {sample}")
				continue
			sample_dir = str(Path(output_abs) / "Intermediate" / sample)
			os.makedirs(sample_dir, exist_ok=True)

			exc_script = str(Path(sample_dir) / f"{sample}_exc.py")
			with open(exc_script, "w") as exc:
				exc.write(_render_exc_script(
					script_path, sample, bam_abs,
					jxn_list, gene_regions, gene_gtfs,
					fasta, buffer, threads,
					output_abs, regtools
				))
			# Write absolute path so the master script runs from any directory
			master.write(f"python3 {exc_script}\n")

	print(f"[setup] Master execution script : {master_exc}")

	# ── Summary script ───────────────────────────────────────────────────────
	summary_script = str(Path(output_abs) / "junctScopeSummarize.py")
	with open(summary_script, "w") as summ:
		summ.write(_render_summary_script(output_abs, jxn_list, regtools))
	print(f"[setup] After all samples finish: python3 {summary_script}")

	return jxn_list, samples


def _render_exc_script(
	script_path, sample, bam_file,
	jxn_list, gene_regions, gene_gtfs,
	fasta, buffer, threads, output, regtools
) -> str:
	"""Render the Python source for a per-sample execution script."""
	return f"""#!/usr/bin/env python3
\"\"\"Auto-generated runner for sample {sample}.\"\"\"
import sys, os
sys.path.insert(0, {repr(str(Path(script_path).parent))})
from junctionScope import run_sample_junction
from pathlib import Path

os.chdir({repr(str(Path(output).resolve()))})

sample   = {repr(sample)}
bam_file = {repr(bam_file)}
fasta    = {repr(fasta)}
buffer   = {buffer}
threads  = {threads}
regtools = {regtools}

jxn_list     = {repr(jxn_list)}
gene_regions = {repr(gene_regions)}
gene_gtfs    = {repr(gene_gtfs)}

for jxn in jxn_list:
	gene = jxn["gene"]
	print(f"[{{sample}}] junction: {{jxn['junction_name']}}")
	run_sample_junction(
		sample      = sample,
		bam_file    = bam_file,
		jxn_entry   = jxn,
		gene_region = gene_regions[gene],
		fasta       = fasta,
		gtf         = gene_gtfs[gene],
		output_dir  = os.path.join({repr(output)}, "Intermediate", sample, jxn["junction_name"]),
		buffer      = buffer,
		threads     = threads,
		regtools    = regtools
	)

print(f"[{{sample}}] all junctions complete.")
"""


def _render_summary_script(output: str, jxn_list: list, regtools: bool = False) -> str:
	"""Render the Python source for the cross-sample summary script."""
	if regtools and validate_regtools():
		header = "\t".join(OUTPUT_COLUMNS)
	else:
		header = "\t".join(OUTPUT_COLUMNS_noReg)
	# output is already absolute (resolved in setup())
	output_abs = str(Path(output).resolve())
	return f"""#!/usr/bin/env python3
\"\"\"Collect per-sample per-junction results into a single summary table
and merge per-barcode best-junction files across junctions and samples.\"\"\"
import glob, os, sys
import numpy as np
from pathlib import Path

# Absolute paths baked in at setup time — script runs correctly from any directory
OUTPUT_DIR  = {repr(output_abs)}
output_file = os.path.join(OUTPUT_DIR, Path(OUTPUT_DIR).name + "_output.txt")
header      = {repr(header)}
jxn_list    = {repr(jxn_list)}

sys.path.insert(0, {repr(str(Path(__file__).resolve().parent))})
from junctionScope import merge_bcd_best

# ── Sample-level summary ────────────────────────────────────────────────────
with open(output_file, "w") as out:
	out.write(header + "\\n")
	pattern = os.path.join(OUTPUT_DIR, "Intermediate", "**", "*.junctScope.txt")
	for txt in sorted(glob.glob(pattern, recursive=True)):
		with open(txt) as fh:
			# Skip the header row written by result_df.to_csv
			_hdr = fh.readline()
			out.write(fh.read())

print(f"Sample summary written to {{output_file}}")

# ── Per-barcode merge ────────────────────────────────────────────────────────
merge_bcd_best(OUTPUT_DIR, jxn_list)
"""



# =============================================================================
# CLI
# =============================================================================

def main():
	parser = argparse.ArgumentParser(
		description="junctionScope — splice junction caller (single-cell aware)",
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog="""
Examples:
  # Setup project (write per-sample scripts, do not run)
  python junctionScope.py -c project.conf --setup-only

  # Setup + run all samples sequentially
  python junctionScope.py -c project.conf

  # Run a single sample (e.g. from a Slurm array)
  python junctionScope.py -c project.conf --sample SampleA --bam /path/to/A.cram
		""",
	)
	parser.add_argument("-c", "--config",  required=True, help="Path to .conf file")
	parser.add_argument("--setup-only",    action="store_true",
						help="Write execution scripts but do not run samples")
	parser.add_argument("--sample",        help="Run a single named sample")
	parser.add_argument("--bam",           help="BAM/CRAM path (required with --sample)")
	parser.add_argument("--merge-barcodes", action="store_true",
						help="Merge per-barcode bcd.best files across junctions/samples and exit")
	parser.add_argument("--run-regtools", action="store_true", default = False,
						help="Run regtools junction extraction and annotation steps")
	args = parser.parse_args()

	config = parse_config(args.config)
	validate_config(config)
	jxn_list, samples = setup(config)

	if args.setup_only:
		return

	if args.merge_barcodes:
		merge_bcd_best(config["OUTPUT"], jxn_list)
		return

	fasta   = config["FASTA"]
	buffer  = int(config["BUFFER"])
	threads = int(config["THREADS"])
	output  = config["OUTPUT"]
	gene_regions = config["_gene_regions"]
	gene_gtfs    = config["_gene_gtfs"]

	if "REGTOOLS" not in config:
		config["REGTOOLS"] = args.run_regtools

	regtools = config["REGTOOLS"]

	def _run_one(sample, bam_file):
		for jxn in jxn_list:
			sample_dir = os.path.join(output, "Intermediate", sample, jxn['junction_name'])
			gene = jxn["gene"]
			print(f"[{sample}] junction: {jxn['junction_name']}")
			run_sample_junction(
				sample      = sample,
				bam_file    = bam_file,
				jxn_entry   = jxn,
				gene_region = gene_regions[gene],
				fasta       = fasta,
				gtf         = gene_gtfs[gene],
				output_dir  = sample_dir,
				buffer      = buffer,
				threads     = threads,
				regtools    = regtools
			)

	if args.sample:
		if not args.bam:
			sys.exit("--sample requires --bam")
		_run_one(args.sample.replace(".", ""), args.bam)
		return

	for sample, bam_file in samples:
		if os.path.isfile(bam_file):
			_run_one(sample, bam_file)
		else:
			print(f"[main] WARNING: {bam_file} not found — skipping {sample}")


# =============================================================================
# Per-barcode merge across junctions and samples
# =============================================================================

def merge_bcd_best(
	output_dir: str,
	jxn_list: list,
	output_file: str = None,
) -> pd.DataFrame:
	"""
	Merge all *.sam.jxn.bcd.best.tsv files across junctions and samples into
	a single wide per-barcode table.

	For each junction, the per-barcode file is loaded and the following
	transformations are applied before merging:
	  - cell_barcode is suffixed with the sample name:
			ACGTACGT-1  →  ACGTACGT-1__<sample>
	  - numeric count/pct columns are prefixed with the junction name:
			reads_per_barcode  →  <junction_name>__reads_per_barcode

	Files are then row-bound across samples into a long table per junction,
	and finally the junction tables are outer-merged on cell_barcode so that
	every barcode appears once with one column group per junction.
	Parameters
	----------
	output_dir  : top-level project output directory (absolute or relative)
	jxn_list    : list of junction dicts (junction_name, gene, junction, nt_seq)
	output_file : path for the merged TSV
				  (default: <output_dir>/<output_dir_name>_bcd_merged.tsv)

	Returns the merged DataFrame (also written to output_file).
	"""
	# Resolve immediately so every downstream path call is stable
	output_dir = str(Path(output_dir).resolve())
	if output_file is None:
		base_name = Path(output_dir).name
		output_file = os.path.join(output_dir, f"{base_name}_bcd_merged.tsv")
	else:
		output_file = str(Path(output_file).resolve())
	# Columns that identify the barcode — kept as-is, not prefixed
	ID_COLS = {"cell_barcode"}
	# Collect one wide-ish DataFrame per junction (rows = barcodes from all samples)
	jxn_frames = {}
	for jxn in jxn_list:
		jxn_name = jxn["junction_name"]
		# Glob pattern: Intermediate/<sample>/<jxn_name>/<sample>.<jxn_name>.sam.jxn.bcd.best.tsv
		pattern = str(
			Path(output_dir) / "Intermediate" / "*" / jxn_name
			/ f"*.{jxn_name}.sam.jxn.bcd.best.tsv"
		)
		files = sorted(glob_files(pattern))
		if not files:
			print(f"[merge_bcd_best] WARNING: no bcd.best files found for "
				  f"junction '{jxn_name}'\n  pattern: {pattern}")
			continue
		per_sample_dfs = []
		for fpath in files:
			fpath = str(Path(fpath).resolve())
			df = pd.read_csv(fpath, sep="\t")
			if df.empty:
				continue
			# Derive sample name from directory structure:
			# Intermediate/<sample>/<jxn_name>/<file>
			sample_name = Path(fpath).parts[-3]
			# Suffix cell_barcode with sample name
			df["cell_barcode"] = (
				df["cell_barcode"].astype(str) + "__" + sample_name
			)
			df = df.drop(columns=['chromosome','junction_start','junction_end'])
			df.columns = df.columns.str.replace(r'_per_barcode|_bcd', '', regex=True)
			# Prefix all value columns with junction name
			rename = {
				c: f"{jxn_name}__{c}"
				#for c in df.columns
				for c in df.columns
				if c not in ID_COLS
			}
			df = df.rename(columns=rename)
			per_sample_dfs.append(df)
		if not per_sample_dfs:
			continue
		# Row-bind all samples for this junction
		jxn_long = pd.concat(per_sample_dfs, ignore_index=True)
		jxn_frames[jxn_name] = jxn_long
		print(f"[merge_bcd_best] {jxn_name}: {len(jxn_long)} barcode rows "
			  f"from {len(per_sample_dfs)} sample(s)")
	if not jxn_frames:
		print("[merge_bcd_best] No data found — check that samples have run "
			  f"and that output_dir is correct:\n  {output_dir}")
		return pd.DataFrame()
	# Outer-merge all junctions on cell_barcode so every barcode appears once
	merged = None
	for jxn_name, df in jxn_frames.items():
		if merged is None:
			merged = df
		else:
			merged = merged.merge(df, on="cell_barcode", how="outer")
	# Sort by first junction's read count descending (if present)
	first_jxn = jxn_list[0]["junction_name"]
	sort_col = f"{first_jxn}__reads_per_barcode"
	if sort_col in merged.columns:
		merged = merged.sort_values(sort_col, ascending=False)
	merged = merged.reset_index(drop=True)
	numeric_cols = merged.select_dtypes(include=[np.number]).columns
	merged[numeric_cols] = merged[numeric_cols].fillna(0)
	merged.to_csv(output_file, sep="\t", index=False)
	print(f"[merge_bcd_best] Merged table written to: {output_file} "
		  f"({len(merged)} barcodes × {len(merged.columns)} columns)")
	return merged


def glob_files(pattern: str) -> list:
	"""Wrapper around glob.glob for testability."""
	import glob
	return glob.glob(pattern, recursive=True)



if __name__ == "__main__":
	main()


