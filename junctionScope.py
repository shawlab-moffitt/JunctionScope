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
import shlex
import json
import copy
import tempfile
from functools import reduce
from datetime import datetime
import scanpy as sc
import scvelo as scv
import cellrank as cr
import anndata as ad
import pyranges as pr
import importlib
import warnings




# =============================================================================
# Config / junction table parsing
# =============================================================================
config_default = {
  "main_config": {
	"input": "",
	"proj_name": "JunctionScope_Output",
	"jxn_table": "",
	"gene": "",
	"jxncoord": "",
	"jxnsource":"",
	"ntseq": "",
	"sc_method":"",
	"bcd_tag":"",
	"umi_tag":"",
	"fasta": "",
	"gtf": "",
	"mode": "sc",
	"regtools": False,
	"buffer": 200,
	"threads": 16,
	"qc_step": False,
	"velocyto": False
  },
  "function_config": {
	"get_nt_seq": {"window": 5,"rev_comp": False},
	"load_fasta_kmers": {"revC": True,"k": 31},
	"extract_jxn_region": {"buffer": 200},
	"filter_to_junction": {"buffer": 5,"filter_barcodes": False},
	"seq_check": {"filter_rows": False,"case_sensitive": False},
	"query_bam_to_fasta": {"revC": True,"k": 31},
	"select_best_junction": {"min_reads": 5,"min_seq_reads": 0},
	"extract_jxn_regtools": {"buffer": 1000},
	"match_sam_to_regtools": {"coord_tol": 1}
  }
}


def deep_merge(d1, d2):
	"""
	Recursively merges two dictionaries. 
	Values in d2 will overwrite or add to values in d1.
	"""
	merged = copy.deepcopy(d1)
	for key, value in d2.items():
		if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
			# Recursively merge sub-dictionaries
			merged[key] = deep_merge(merged[key], value)
		else:
			# Overwrite or add the new key/value
			merged[key] = copy.deepcopy(value)
	return merged

def parse_config(config_path: str) -> dict:
	"""Parse a shell-style KEY="value" config file into a dict."""
	ext = Path(config_path).suffix
	if ext == '.json':
		config = json.load(open(config_path,'r'))
		config_out = deep_merge(config_default,config)
		return config_out
	else:
		config = {}
		with open(config_path) as fh:
			for line in fh:
				line = line.strip()
				if not line or line.startswith("#"):
					continue
				line = re.sub(r'\s*#.*$', '', line)
				m = re.match(r'^(\w+)\s*=\s*"?([^"]*)"?\s*$', line)
				if m:
					config[m.group(1).lower()] = m.group(2)
		if 'proj_name' not in config and 'output' in config:
			config["proj_name"] = config.pop("output")
		config_out = deep_merge(config_default,config)
		return config_out


def validate_config(config: dict, require_input: bool = True):
	"""
	Validate main_config dict.
	- proj_name always required
	- input required unless running single-sample mode
	- fasta required when regtools=True OR any junction is missing an NT sequence
	- gtf required when regtools=True
	"""
	required_always = ["proj_name"]
	if require_input:
		required_always.append("input")
	# fasta needed for regtools AND for deriving missing NT sequences
	need_fasta = (
		config.get("regtools", False) or
		config.get("ntseq", "") == "" and config.get("jxn_table", "") == ""
	)
	if need_fasta:
		required_always.append("fasta")
	if config.get("regtools", False) or config.get("velocyto", False):
		required_always.append("gtf")
	missing = [k for k in required_always if not config.get(k)]
	if missing:
		sys.exit(f"[junctionScope] Missing required config keys: {', '.join(missing)}")
	has_single = config.get("jxncoord", "") != "" and config.get("gene", "") != ""
	has_table  = config.get("jxn_table", "") != ""
	if not has_single and not has_table:
		sys.exit(
			"[junctionScope] Config must supply either:\n"
			"  jxncoord + gene   (single-junction mode)\n"
			"  jxn_table         (multi-junction table mode)"
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
	_hdr = open(path).readline()
	if ':' in _hdr and '-' in _hdr:
		df = pd.read_csv(path, sep="\t", comment="#", names = ['junction_name','gene','junction','nt_seq'])
	else:
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
		df["nt_seq"] = ""
	df["nt_seq"] = df["nt_seq"].fillna("").astype(str)
	return df[["junction_name", "gene", "junction", "nt_seq"]].to_dict("records")


def config_to_jxn_list(config: dict) -> list[dict]:
	"""
	Return a normalised list of junction dicts regardless of config mode.
	Each dict has: junction_name, gene, junction, nt_seq
	"""
	if config["jxn_table"] != "":
		return parse_jxn_table(config["jxn_table"])
	else:
		return [{
			"junction_name": config.get("jxn_name", config["gene"]),
			"gene":          config["gene"],
			"junction":      config["jxncoord"],
			"nt_seq":        config["ntseq"],
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
# Check for format of chromosome name - with or without 'chr''
# =============================================================================

def contig_alias_key(contig: str) -> str:
	"""
	Return a normalized key for common chromosome aliases.

	Examples
	--------
	chr1 and 1       -> 1
	chrX and X       -> X
	chrM, M, and MT  -> MT
	"""
	if contig is None:
		return None

	try:
		if pd.isna(contig):
			return None
	except (TypeError, ValueError):
		pass

	name = str(contig).strip()

	if name.lower().startswith("chr"):
		name = name[3:]

	name = name.upper()

	if name in {"M", "MT"}:
		return "MT"

	return name


def resolve_contig_name(
	contig: str,
	references,
	source_label: str = "reference",
) -> str:
	"""
	Resolve a requested chromosome name against contigs found in a
	BAM, CRAM, or FASTA header.

	Exact matches are preferred. Common aliases such as chrX/X and
	chr1/1 are used only when necessary.
	"""
	requested = str(contig).strip()
	refs = tuple(str(ref) for ref in references)

	# Best case: exact match
	if requested in refs:
		return requested

	# Case-insensitive exact match
	case_matches = [
		ref for ref in refs
		if ref.lower() == requested.lower()
	]

	if len(case_matches) == 1:
		return case_matches[0]

	# Direct chr-prefix counterpart
	if requested.lower().startswith("chr"):
		direct_alias = requested[3:]
	else:
		direct_alias = f"chr{requested}"

	if direct_alias in refs:
		return direct_alias

	# General alias comparison
	target_key = contig_alias_key(requested)

	alias_matches = [
		ref for ref in refs
		if contig_alias_key(ref) == target_key
	]

	if len(alias_matches) == 1:
		return alias_matches[0]

	if len(alias_matches) > 1:
		raise ValueError(
			f"Ambiguous contig alias {requested!r} for {source_label}; "
			f"possible matches: {alias_matches}"
		)

	preview = ", ".join(refs[:12])

	if len(refs) > 12:
		preview += ", ..."

	raise ValueError(
		f"Contig {requested!r} was not found in {source_label}. "
		f"Available contigs begin with: {preview}"
	)



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
	#left_start  = start - window - 1
	left_start  = start - window
	#left_end    = start - 1
	left_end    = start
	right_start = end
	right_end   = end + window
	with pysam.FastaFile(fasta) as fa:
		resolved_chrom = resolve_contig_name(
			chrom,
			fa.references,
			source_label=f"FASTA {fasta}",
		)

		left_seq = fa.fetch(
			resolved_chrom,
			left_start,
			left_end,
		)

		right_seq = fa.fetch(
			resolved_chrom,
			right_start,
			right_end,
		)
	# with pysam.FastaFile(fasta) as fa:
	# 	left_seq  = fa.fetch(chrom, left_start, left_end)
	# 	right_seq = fa.fetch(chrom, right_start, right_end)
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
		if len(seq) < k:
			k = len(seq)
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

# def extract_jxn_region(
# 	input_file: str,
# 	region: str,
# 	output_sam: str,
# 	buffer: int = 200,
# 	threads: int = 1,
# ):
# 	"""
# 	Subset a BAM/CRAM to a buffered region and write a SAM file.
# 	Auto-detects BAM vs CRAM from file suffix.
# 	"""
# 	suff = Path(input_file).suffix.lower()
# 	mode = "rb" if suff == ".bam" else "rc"
# 	buffered = add_buffer(region, buffer)
# 	with pysam.AlignmentFile(input_file, mode, threads=threads) as bam, \
# 		 pysam.AlignmentFile(output_sam, "w", header=bam.header) as out:
# 		for read in bam.fetch(region=buffered):
# 			out.write(read)


def extract_jxn_region(
	input_file: str,
	region: str,
	output_sam: str,
	buffer: int = 200,
	threads: int = 1,
) -> str:
	"""
	Subset a BAM/CRAM to a buffered region and write a SAM file.

	The requested chromosome name is resolved against the alignment header,
	allowing common aliases such as chrX/X and chr1/1.

	Returns
	-------
	str
		The unbuffered region using the chromosome name present in the BAM.
	"""
	suff = Path(input_file).suffix.lower()

	if suff == ".bam":
		mode = "rb"
	elif suff == ".cram":
		mode = "rc"
	else:
		raise ValueError("input_file must end in .bam or .cram")

	requested_chrom, start, end = parse_region(region)

	with pysam.AlignmentFile(
		input_file,
		mode,
		threads=threads,
	) as bam, pysam.AlignmentFile(
		output_sam,
		"w",
		header=bam.header,
	) as out:

		resolved_chrom = resolve_contig_name(
			requested_chrom,
			bam.references,
			source_label=f"alignment {input_file}",
		)

		if resolved_chrom != requested_chrom:
			print(
				f"  [contig] BAM uses {resolved_chrom!r}; "
				f"requested coordinate used {requested_chrom!r}"
			)

		# JunctionScope coordinates are treated as 1-based.
		# pysam fetch coordinates are 0-based, half-open.
		fetch_start = max(0, start - buffer - 1)
		fetch_end = end + buffer

		for read in bam.fetch(
			resolved_chrom,
			fetch_start,
			fetch_end,
		):
			out.write(read)

	return f"{resolved_chrom}:{start}-{end}"



# =============================================================================
# SAM parser  (produces one row per junction per read, with barcode)
# =============================================================================

# bug3
BARCODE_TAGS = ["CB", "XC"]
UMI_TAGS = ["UB", "XM"]
def get_first_tag(read, tags):
	if isinstance(read, pysam.AlignedSegment):
		for tag in tags:
			if read.has_tag(tag):
				return tag, read.get_tag(tag)
		return None
	else:
		if isinstance(read, (list, tuple)):
			optional_fields = read
		elif not read or read.startswith("@"):
			return None
		else:
			fields = read.rstrip("\n").split("\t")
			# SAM optional fields start at column 12, index 11
			optional_fields = fields[11:]
		for wanted_tag in tags:
			prefix = f"{wanted_tag}:"
			for field in optional_fields:
				if field.startswith(prefix):
					# Format is TAG:TYPE:VALUE, e.g. CB:Z:AAAC...
					parts = field.split(":", 2)
					if len(parts) == 3:
						return wanted_tag, parts[2]


def parse_sam_per_junction(sam_file: str,
	mode: str,
	bcd_tag: str = None,
	umi_tag: str = None
) -> pd.DataFrame:
	"""
	Parse a SAM file into a per-junction DataFrame.

	One row per (read × junction).  Reads with no N in CIGAR get one row
	with junction_start/end = None.  Columns include cell_barcode (CB:Z: tag).
	"""
	if mode != 'sc' and mode != 'bulk':
		print("[parse_sam_per_junction] 'mode' argument not found ...")
		print("[parse_sam_per_junction] mode = 'sc' for single-cell data with barcodes")
		sys.exit("[parse_sam_per_junction] mode = 'bulk' for bulk RNASeq data")
		return
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
			if mode == 'sc':
				# Extract cell barcode from auxiliary tags
				_bcd_tags  = ["CB", "XC"] if bcd_tag is None else [bcd_tag]
				_umi_tags  = ["UB", "XM"] if umi_tag is None else [umi_tag]
				_bcd_result = get_first_tag(fields[11:], _bcd_tags)
				_umi_result = get_first_tag(fields[11:], _umi_tags)
				_, barcode = _bcd_result if _bcd_result else (None, None)
				_, umi     = _umi_result if _umi_result else (None, None)
				# Does not count read if now processed barcode or umi
				if barcode == None or umi == None:
					continue
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
				"n_junctions":   n_junctions,
				"multi_junction":n_junctions > 1,
			}
			if mode == 'sc':
				base_row['cell_barcode'] = barcode
				base_row['umi']          = umi
				base_row['barcode_umi']  = str(barcode) + '_' + str(umi)
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
	mode: str,
	buffer: int = 2,
	filter_barcodes: bool = False,
	chr_col: str = 'chromosome',
	start_col: str = 'junction_start',
	end_col: str = 'junction_end',
	barcode_col: str = 'cell_barcode'
	):
	"""Filter parsed sam file to junction coordinates of interest."""
	chrom, start, end = parse_region(coord)
	col_check = parsed_sam.columns
	if chr_col in col_check and start_col in col_check and end_col in col_check:
		target_contig_key = contig_alias_key(chrom)
		chrom_mask = (
			parsed_sam[chr_col]
			.map(contig_alias_key)
			.eq(target_contig_key)
		)
		parsed_sam_coord = parsed_sam[
			chrom_mask
			& parsed_sam[start_col].between(
				start - buffer,
				start + buffer,
			)
			& parsed_sam[end_col].between(
				end - buffer,
				end + buffer,
			)
		]
		# parsed_sam_coord = parsed_sam[
		# 							(parsed_sam[chr_col] == chrom) &
		# 							(parsed_sam[start_col].between(start-buffer,start+buffer)) &
		# 							(parsed_sam[end_col].between(end-buffer,end+buffer))]
		if filter_barcodes or mode == 'bulk':
			return pd.DataFrame(parsed_sam_coord)
		else:
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

def query_seq_to_bam(
	input_file: str,
	mode: str,
	fasta_path: str = None,
	seq: str = None,
	revC: bool = True,
	threads: int = 8,
	k: int = 31,
):
	total_reads = 0
	matched_reads = 0
	if len(seq) > 0:
		with tempfile.NamedTemporaryFile(mode='w+t', delete=False) as fasta_temp:
			fasta_temp.write(">InputSeq\n")
			fasta_temp.write(seq+'\n')
			fasta_temp.seek(0)
			fasta_temp.flush()
			print(f"Writing temporary fasta: {fasta_temp.name}")
		fasta_path = fasta_temp.name
	ref_fa, kmer_index = load_fasta_kmers(
		fasta_path,
		revC,
		k
	)
	os.remove(fasta_temp.name)
	suff = Path(input_file).suffix
	if suff == ".bam":
		input_mode = "rb"
	elif suff == ".cram":
		input_mode = "rc"
	else:
		raise ValueError("Input must be BAM or CRAM")
	bam = pysam.AlignmentFile(
		input_file,
		input_mode,
		threads=threads
	)
	matches = []
	for read in bam.fetch(until_eof=True):
		total_reads += 1
		if read.query_sequence is None:
			continue
		if mode == 'sc':
			bc = read.get_tag("CB") if read.has_tag("CB") else None
			umi = read.get_tag("UB") if read.has_tag("UB") else None
		query = read.query_sequence.upper()
		matching_key = next((k for k, v in ref_fa.items() if v in query), None)
		if matching_key == None:
			continue
		else:
			matched_reads += 1
			matching_val = ref_fa[matching_key]
			strand = (
				"-"
				if matching_key.endswith("_revC")
				else "+"
			)
			clean_ref = re.sub(
				r"_revC$",
				"",
				matching_key
			)
			if mode == 'sc':
				matches.append({
					"read_name": read.query_name,
					"mapped_read": read.is_mapped,
					"barcode": bc,
					"umi": umi,
					"ref_query_name": clean_ref,
					"ref_query": matching_val,
					"strand": strand,
					"read_sequence": query
				})
			else:
				matches.append({
					"read_name": read.query_name,
					"mapped_read": read.is_mapped,
					"ref_query_name": clean_ref,
					"ref_query": matching_val,
					"strand": strand,
					"read_sequence": query
				})
	bam.close()
	matches_df = pd.DataFrame(matches)
	print(f"{matched_reads}/{total_reads} matched")
	return matches_df


def query_bam_to_fasta(
	input_file: str,
	fasta_path: str,
	mode: str,
	revC: bool = True,
	threads: int = 8,
	k: int = 31,
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
		input_mode = "rb"
	elif suff == ".cram":
		input_mode = "rc"
	else:
		raise ValueError("Input must be BAM or CRAM")
	bam = pysam.AlignmentFile(
		input_file,
		input_mode,
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
		if mode == 'sc':
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
				if mode == 'sc':
					matches.append({
						"read_name": read.query_name,
						"barcode": bc,
						"umi": umi,
						"reference": clean_ref,
						"strand": strand,
						"sequence": query
					})
				else:
					matches.append({
						"read_name": read.query_name,
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
	mode: str = "sc"
):
	"""
	Summarize junction support from parsed SAM data.

	Parameters
	----------
	parsed_sam : pd.DataFrame
		Parsed junction-level SAM dataframe.
	mode : str
		"sc"   = single-cell RNA-seq
		"bulk" = bulk RNA-seq
	"""
	if mode not in ["sc", "bulk"]:
		raise ValueError("mode must be either 'sc' or 'bulk'")

	# Treat NA multi_junction values as False for counting
	parsed_sam = parsed_sam.copy()
	parsed_sam["multi_junction"] = (
		parsed_sam["multi_junction"]
		.fillna(False)
		.astype(bool)
	)

	agg_dict = {
		# total junction support
		"sam_read_count": ("read_id", "nunique"),

		# sequence-supported reads
		"seq_read_count": (
			"read_id",
			lambda x: (
				parsed_sam.loc[x.index]
				.query("junction_seq_found == True")["read_id"]
				.nunique()
			)
		),

		# total multi-junction observations
		"multi_junction_reads": ("multi_junction", "sum"),

		# unique reads marked multi-junction
		"multi_junction_reads_uniq": (
			"read_id",
			lambda x: (
				parsed_sam.loc[x.index]
				.query("multi_junction == True")["read_id"]
				.nunique()
			)
		)
	}

	if mode == "sc":
		agg_dict.update({
			"sam_n_cells": ("cell_barcode", "nunique"), # bug6
			"sam_n_umi": ("barcode_umi", "nunique"),

			"seq_n_cells": (
				"cell_barcode",
				lambda x: (
					parsed_sam.loc[x.index]
					.query("junction_seq_found == True")["cell_barcode"]
					.nunique()
				)
			),

			"seq_n_umi": (
				"barcode_umi",
				lambda x: (
					parsed_sam.loc[x.index]
					.query("junction_seq_found == True")["barcode_umi"]
					.nunique()
				)
			)
		})

	parsed_sam_summ = (
		parsed_sam
		.dropna(subset=["junction_start", "junction_end"])
		.groupby(
			["chromosome", "junction_start", "junction_end"]
		)
		.agg(**agg_dict)
		.reset_index()
	)

	# Read-level percentages
	parsed_sam_summ["pct_reads_seq"] = (
		parsed_sam_summ["seq_read_count"]
		/ parsed_sam_summ["sam_read_count"]
	)

	# Single-cell only metrics
	if mode == "sc":
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
	mode: str,
	min_reads: int = 5,
	min_seq_reads: int = 1,
	) -> pd.DataFrame:
	"""Select highest-confidence junction from summary table."""
	df = junct_summ.copy()
	# optional pre-filter
	df = df[
		(df["sam_read_count"] >= min_reads) &
		(df["seq_read_count"] >= min_seq_reads)
	]
	if df.empty:
		df = pd.DataFrame(columns=df.columns)
		#raise ValueError(
		#	"No junctions pass filtering criteria."
		#)
	else:
		if mode == 'sc':
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
		else:
			df = df.sort_values(
				by=[
					"seq_read_count",
					"pct_reads_seq",
					"sam_read_count"
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


def _valid_tag_value(value) -> bool:
	"""
	Return True if a barcode/UMI tag value is usable.
	Avoids counting missing placeholders as real barcodes/UMIs,
	because computers love making 'None' look like data.
	"""
	if value is None:
		return False
	value = str(value).strip()
	if value == "":
		return False
	return value.lower() not in {
		"none", "nan", "na", "n/a", "<na>", "null", "*"
	}

def calc_pct(n,d,r: int = 4):
	if r > 0:
		return(np.round(np.array(n)/np.array(d)*100, r))
	else:
		return((np.array(n)/np.array(d))*100)

def count_sample_stats(
	input_file: str,
	mode: str = 'sc',
	verbose: bool = True,
	threads: int = 1,
	bcd_tag: str = None,
	umi_tag: str = None,
	reference_fasta: str = None,
) -> dict:
	if (verbose):
		print(datetime.now())
	suff = Path(input_file).suffix.lower()
	if suff == ".bam":
		mode_input = "rb"
	elif suff == ".cram":
		mode_input = "rc"
	else:
		raise ValueError("input_file must be a BAM or CRAM file")
	read_count = 0
	junction_read_count = 0
	junction_count = 0
	if mode == 'sc':
		_bcd_tags = ["CB", "XC"] if bcd_tag is None else [bcd_tag]
		_umi_tags = ["UB", "XM"] if umi_tag is None else [umi_tag]
	barcodes = set()
	umis = set()
	barcode_umis = set()
	barcodes_wJunc = set()
	umis_wJunc = set()
	barcode_umis_wJunc = set()
	open_kwargs = {
		"threads": threads
	}
	if reference_fasta is not None:
		open_kwargs["reference_filename"] = reference_fasta
	with pysam.AlignmentFile(input_file, mode_input, **open_kwargs) as bam:
		for read in bam.fetch(until_eof=True):
			if read.is_secondary or read.is_supplementary or not read.cigartuples:
				continue
			read_count += 1
			n_junctions = sum(1 for op, _ in read.cigartuples if op == 3)
			has_junction = n_junctions > 0
			if has_junction:
				junction_read_count += 1
				junction_count += n_junctions
			if mode == 'sc':
				_bcd_result = get_first_tag(read, _bcd_tags)
				_umi_result = get_first_tag(read, _umi_tags)
				_, barcode = _bcd_result if _bcd_result else (None, None)
				_, umi = _umi_result if _umi_result else (None, None)
				barcode_ok = _valid_tag_value(barcode)
				umi_ok = _valid_tag_value(umi)
				if barcode_ok:
					barcodes.add(str(barcode))
				if umi_ok:
					umis.add(str(umi))
				if barcode_ok and umi_ok:
					barcode_umis.add(f"{barcode}_{umi}")
				if has_junction:
					if barcode_ok:
						barcodes_wJunc.add(str(barcode))
					if umi_ok:
						umis_wJunc.add(str(umi))
					if barcode_ok and umi_ok:
						barcode_umis_wJunc.add(f"{barcode}_{umi}")
	if (verbose):
		print(datetime.now())
	if mode == 'sc':
		count_col = ["Read Count","Junction Count","Cell Barcode Count","UMI Count"]
		total_col = [read_count,junction_count,len(barcodes),len(umis)]
		wJunc_col = [junction_read_count,junction_count,len(barcodes_wJunc),len(umis_wJunc)]
	else:
		count_col = ["Read Count","Junction Count"]
		total_col = [read_count,junction_count]
		wJunc_col = [junction_read_count,junction_count]
	pct_junc_col = calc_pct(wJunc_col,total_col)
	qc_dict = {
		'Count':count_col,
		'Total':total_col,
		'With Junction':wJunc_col,
		'Percent with Junction':pct_junc_col
	}
	return(qc_dict)



def sample_junc_qc(
	jxn_summ: str,
	qc_dict: dict,
	output: str = None,
	mode: str = 'sc',
	verbose: bool = True
	):
	total_col = qc_dict.get("Total")
	if isinstance(jxn_summ,pd.DataFrame):
		summ = jxn_summ
	else:
		summ=pd.read_csv(jxn_summ,sep = '\t')
		if is_empty(output):
			base='.'.join(jxn_summ.split('.')[:2])
			output=f"{base}.jxn.qc.tsv"
	junc_reads = summ['targetJxn_read_count'].iloc[0]
	junc_seq_reads = summ['targetJxn_seq_read_count'].iloc[0]
	if mode == 'sc':
		cell_reads = summ['targetJxn_cell_count'].iloc[0]
		cell_seq_reads = summ['targetJxn_seq_cell_count'].iloc[0]
		umi_reads = summ['targetJxn_umi_count'].iloc[0]
		umi_seq_reads = summ['targetJxn_seq_umi_count'].iloc[0]
		wTJunc_col = [junc_reads,junc_reads,cell_reads,umi_reads]
		wTSJunc_col = [junc_seq_reads,junc_seq_reads,cell_seq_reads,umi_seq_reads]
	else:
		wTJunc_col = [junc_reads,junc_reads]
		wTSJunc_col = [junc_seq_reads,junc_seq_reads]
	pct_Tjunc_col = calc_pct(wTJunc_col,total_col,0)
	pct_TSjunc_col = calc_pct(wTSJunc_col,total_col,0)
	qc_dict.update({
		'With Target Junction':wTJunc_col,
		'Percent with Target Junction':pct_Tjunc_col,
		'With Target Junction Sequence':wTSJunc_col,
		'Percent with Target Junction Sequence':pct_TSjunc_col
	})
	qc_dict_df = pd.DataFrame(qc_dict)
	if verbose:
		print(qc_dict_df)
	qc_dict_df.to_csv(output, sep = '\t', index = False)
	#return(qc_dict)



QC_ROWS = ["Read Count", "Junction Count", "Cell Barcode Count", "UMI Count"]

# filename stub used for each of the 4 output files
_ROW_TO_STUB = {
	"Read Count": "read_qc_summary",
	"Junction Count": "junction_qc_summary",
	"Cell Barcode Count": "cell_qc_summary",
	"UMI Count": "umi_qc_summary",
}


def _read_qc_table(qc_path: Path) -> pd.DataFrame:
	"""Read a single *.jxn.qc.tsv file, indexed by the Count column."""
	return pd.read_csv(qc_path, sep="\t", index_col=0)


def _write_qc_summaries(row_frames: dict, index_name: str, output_dir: Path, prefix: str):
	"""
	row_frames: dict of {qc_row_name: pd.DataFrame} where each DataFrame has
				rows = the merge dimension (samples or junctions) and
				columns = the 8 QC metrics.
	Writes 4 tsv files named <prefix>.<stub>.tsv
	"""
	output_dir.mkdir(parents=True, exist_ok=True)
	written = []
	for qc_row, df in row_frames.items():
		df.index.name = index_name
		out_path = output_dir / f"{prefix}.{_ROW_TO_STUB[qc_row]}.tsv"
		df.to_csv(out_path, sep="\t")
		written.append(out_path)
	return written


def merge_qc_across_samples(base_dir: str, junction: str, output_dir: str):
	"""
	Merge QC tables for a single junction across all samples that contain it.

	base_dir/<sample>/<junction>/*.jxn.qc.tsv is searched for every sample.
	Output rows = sample names.
	"""
	base_dir = Path(base_dir)
	output_dir = Path(output_dir)
	qc_files = sorted(base_dir.glob(f"*/{junction}/*.jxn.qc.tsv"))
	if not qc_files:
		raise FileNotFoundError(
			f"No qc.tsv files found under {base_dir}/*/{junction}/"
		)
	row_frames = {qc_row: {} for qc_row in QC_ROWS}
	for qc_path in qc_files:
		sample = qc_path.parent.parent.name  # base_dir/<sample>/<junction>/file
		table = _read_qc_table(qc_path)
		for qc_row in QC_ROWS:
			if qc_row in table.index:
				row_frames[qc_row][sample] = table.loc[qc_row]
			else:
				print(f"Warning: '{qc_row}' missing in {qc_path}, skipping for {sample}")
	merged = {
		qc_row: pd.DataFrame.from_dict(samples, orient="index")
		for qc_row, samples in row_frames.items()
	}
	return _write_qc_summaries(merged, index_name="sample", output_dir=output_dir, prefix=junction)


def merge_qc_within_sample(base_dir: str, sample: str, output_dir: str):
	"""
	Merge QC tables across all junctions found within a single sample.

	base_dir/<sample>/<junction>/*.jxn.qc.tsv is searched for every junction.
	Output rows = junction names.
	"""
	base_dir = Path(base_dir)
	output_dir = Path(output_dir)
	qc_files = sorted((base_dir / sample).glob("*/*.jxn.qc.tsv"))
	if not qc_files:
		raise FileNotFoundError(f"No qc.tsv files found under {base_dir}/{sample}/*/")
	row_frames = {qc_row: {} for qc_row in QC_ROWS}
	for qc_path in qc_files:
		junction = qc_path.parent.name  # base_dir/sample/<junction>/file
		table = _read_qc_table(qc_path)
		for qc_row in QC_ROWS:
			if qc_row in table.index:
				row_frames[qc_row][junction] = table.loc[qc_row]
			else:
				print(f"Warning: '{qc_row}' missing in {qc_path}, skipping for {junction}")
	merged = {
		qc_row: pd.DataFrame.from_dict(junctions, orient="index")
		for qc_row, junctions in row_frames.items()
	}
	return _write_qc_summaries(merged, index_name="junction", output_dir=output_dir, prefix=sample)




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
# Velocyto
# =============================================================================


def validate_velocyto(verbose: bool = False) -> bool:
	# Check executable exists
	velocyto_path = shutil.which("velocyto")
	if velocyto_path is None:
		if verbose:
			print("velocyto not found in PATH")
		return False
	try:
		# Check main executable works
		result = subprocess.run(
			["velocyto"],
			capture_output=True,
			text=True,
			timeout=10
		)
		output = (
			result.stdout +
			result.stderr
		).lower()
		if "velocity" not in output:
			if verbose:
				print("velocyto executable found but output unexpected")
			return False
		# Check extract command
		extract_result = subprocess.run(
			["velocyto", "run"],
			capture_output=True,
			text=True,
			timeout=10
		)
		extract_output = (
			extract_result.stdout +
			extract_result.stderr
		).lower()
		extract_ok = (
			"usage" in extract_output
			or "bam" in extract_output
		)
		if not extract_ok:
			if verbose:
				print("velocyto unavailable")
		return extract_ok
	except Exception as e:
		if verbose:
			print(f"velocyto validation failed: {e}")
		return False





def bam_to_loom(
	bam_file: str,
	gtf_ref: str,
	out_dir: str,
	threads: int = 16,
) -> Path:
	"""
	Generate a velocyto loom file from a BAM file.
	Equivalent to:
		velocyto run -@ 16 -o OUT_DIR BAM_FILE GTF_REF
	"""
	bam_file = Path(bam_file).expanduser().resolve()
	gtf_ref = Path(gtf_ref).expanduser().resolve()
	out_dir = Path(out_dir).expanduser().resolve()
	if not bam_file.is_file():
		raise FileNotFoundError(f"BAM file not found: {bam_file}")
	if not gtf_ref.is_file():
		raise FileNotFoundError(f"GTF file not found: {gtf_ref}")
	out_dir.mkdir(parents=True, exist_ok=True)
	command = [
		"velocyto",
		"run",
		"-@",
		str(threads),
		"-o",
		str(out_dir),
		str(bam_file),
		str(gtf_ref),
	]
	subprocess.run(
		command,
		check=True,
	)
	expected_loom = out_dir / f"{bam_file.stem}.loom"
	if expected_loom.is_file():
		return expected_loom
	loom_files = sorted(
		out_dir.glob("*.loom"),
		key=lambda path: path.stat().st_mtime,
		reverse=True,
	)
	if len(loom_files) == 1:
		print(
			"[bam_to_loom] Velocyto output filename differed "
			f"from the expected name: {loom_files[0]}"
		)
		return loom_files[0]
	if not loom_files:
		raise FileNotFoundError(
			"Velocyto finished without producing a loom file in "
			f"{out_dir}"
		)
	raise RuntimeError(
		"Multiple loom files were found and the expected output "
		f"{expected_loom.name} was not present:\n"
		+ "\n".join(str(path) for path in loom_files)
	)





	# subprocess.run(command, check=True)
	# # Without -e/--sampleid, velocyto generally derives the sample ID
	# # from the BAM filename.
	# expected_loom = out_dir / f"{bam_file.stem}.loom"
	# return expected_loom




def add_loom_intron_exon_ratio(
	loom_file: str,
	barcode_summary=None,
	sample_name: str = None,
	output_file: str = None,
	barcode_col: str = "cell_barcode",
	sample_col: str = "sample_name",
	pseudocount: float = 0.0,
) -> pd.DataFrame:
	"""
	Summarize velocyto spliced/unspliced counts per barcode and optionally
	append them to an existing barcode summary table.
	The reported intron_exon_ratio is:
		(unspliced_molecules + pseudocount)
		/ (spliced_molecules + pseudocount)
	Ambiguous molecules, when present, are reported but are not included
	in the ratio.
	Parameters
	----------
	loom_file : str
		Path to a velocyto-compatible .loom file.
	barcode_summary : pandas.DataFrame, str, Path, or None
		Existing barcode summary DataFrame or TSV path. If omitted, only
		the loom-derived barcode table is returned.
	sample_name : str, optional
		Sample name to attach to the loom barcodes. If barcode_summary
		contains exactly one sample, its sample name is used automatically.
	output_file : str, optional
		Optional TSV output path.
	barcode_col : str
		Barcode column in barcode_summary.
	sample_col : str
		Sample-name column in barcode_summary.
	pseudocount : float
		Optional pseudocount used only for intron_exon_ratio.
		With the default of 0, cells with zero spliced molecules receive NA
		rather than infinity.
	Returns
	-------
	pandas.DataFrame
		Loom-derived barcode summary, or the supplied barcode summary with
		loom metrics appended.
	"""
	loom_path = Path(loom_file).resolve()
	if not loom_path.is_file():
		raise FileNotFoundError(f"Loom file not found: {loom_path}")
	if pseudocount < 0:
		raise ValueError("pseudocount must be >= 0")
	# Read as cells × genes.
	#
	# X_name="spliced" places the spliced matrix in adata.X in some Scanpy
	# versions, while other readers may retain it in adata.layers. The
	# helper below accommodates either arrangement.
	with warnings.catch_warnings():
		warnings.filterwarnings(
			"ignore",
			message="Variable names are not unique.*",
		)
		adata = sc.read_loom(
			str(loom_path),
			sparse=True,
			X_name="spliced",
			obs_names="CellID",
			var_names="Gene",
		)
	# Preserve the original gene names before AnnData edits them.
	adata.var["gene_symbol_original"] = adata.var_names.astype(str)
	if not adata.var_names.is_unique:
		n_duplicate_gene_rows = int(adata.var_names.duplicated().sum())
		print(
			"[add_loom_intron_exon_ratio] "
			f"Found {n_duplicate_gene_rows} duplicated gene-name rows in loom; "
			"making adata.var_names unique for internal AnnData compatibility."
		)
		adata.var_names_make_unique()
	# adata = sc.read_loom(
	# 	str(loom_path),
	# 	sparse=True,
	# 	X_name="spliced",
	# 	obs_names="CellID",
	# 	var_names="Gene",
	# )
	layer_map = {
		str(layer_name).lower(): layer_name
		for layer_name in adata.layers.keys()
	}
	def _sum_layer(
		layer_name: str,
		allow_x: bool = False,
	) -> np.ndarray:
		"""
		Sum one loom layer across genes for each barcode.
		"""
		actual_layer_name = layer_map.get(layer_name.lower())
		if actual_layer_name is not None:
			matrix = adata.layers[actual_layer_name]
		elif allow_x:
			matrix = adata.X
		else:
			available = ", ".join(map(str, adata.layers.keys()))
			raise ValueError(
				f"Required loom layer '{layer_name}' was not found. "
				f"Available layers: {available or '[none]'}"
			)
		return (
			np.asarray(matrix.sum(axis=1))
			.ravel()
			.astype(float)
		)
	spliced = _sum_layer(
		"spliced",
		allow_x=True,
	)
	unspliced = _sum_layer(
		"unspliced",
		allow_x=False,
	)
	# --------------------------------------------------------------
	# Ratio: unspliced / spliced
	# --------------------------------------------------------------
	ratio_numerator = unspliced + pseudocount
	ratio_denominator = spliced + pseudocount
	intron_exon_ratio = np.divide(
		ratio_numerator,
		ratio_denominator,
		out=np.full(
			spliced.shape,
			np.nan,
			dtype=float,
		),
		where=ratio_denominator > 0,
	)
	# --------------------------------------------------------------
	# Extract and clean loom CellID values
	# --------------------------------------------------------------
	loom_cell_ids = pd.Series(
		adata.obs_names.astype(str),
		dtype="string",
	)
	# Typical velocyto IDs can resemble:
	#     SampleA:AAACCTGAGAGCTGCA-1
	#     SampleA:AAACCTGAGAGCTGCAx
	#
	# Keep the original CellID but make a barcode-only version.
	loom_barcodes = (
		loom_cell_ids
		.str.rsplit(":", n=1)
		.str[-1]
		.str.replace(r"x$", "", regex=True)
	)
	ratio_df = pd.DataFrame({
		"loom_cell_id": loom_cell_ids,
		barcode_col: loom_barcodes,
		"spliced_molecules": spliced,
		"unspliced_molecules": unspliced,
		"intron_exon_ratio": intron_exon_ratio
	})
	# --------------------------------------------------------------
	# Barcode normalization used only for matching
	# --------------------------------------------------------------
	def _barcode_key(values: pd.Series) -> pd.Series:
		"""
		Normalize barcode forms for matching within a sample.
		Handles:
			BARCODE-1
			BARCODE
			BARCODEx
			sample:BARCODE-1
		The original barcode column is not altered.
		"""
		return (
			values.astype("string")
			.str.strip()
			.str.rsplit(":", n=1)
			.str[-1]
			.str.replace(r"x$", "", regex=True)
			.str.replace(r"-\d+$", "", regex=True)
			.str.upper()
		)
	ratio_df["_barcode_key"] = _barcode_key(
		ratio_df[barcode_col]
	)
	# --------------------------------------------------------------
	# Return loom summary without merging
	# --------------------------------------------------------------
	if barcode_summary is None:
		if sample_name is not None:
			ratio_df.insert(
				1,
				sample_col,
				sample_name,
			)
		result = ratio_df.drop(
			columns="_barcode_key"
		)
	# --------------------------------------------------------------
	# Append metrics to an existing barcode summary
	# --------------------------------------------------------------
	else:
		if isinstance(barcode_summary, (str, Path)):
			summary = pd.read_csv(
				barcode_summary,
				sep="\t",
			)
		elif isinstance(barcode_summary, pd.DataFrame):
			summary = barcode_summary.copy()
		else:
			raise TypeError(
				"barcode_summary must be a pandas DataFrame, "
				"TSV path, or None"
			)
		if barcode_col not in summary.columns:
			raise ValueError(
				f"barcode_summary is missing barcode column "
				f"'{barcode_col}'"
			)
		merge_cols = ["_barcode_key"]
		# When the summary has a sample column, merge by both sample
		# and barcode to prevent accidental cross-sample matches.
		if sample_col in summary.columns:
			summary_samples = (
				summary[sample_col]
				.dropna()
				.astype(str)
				.unique()
			)
			if sample_name is None:
				if len(summary_samples) == 1:
					sample_name = summary_samples[0]
				else:
					raise ValueError(
						"barcode_summary contains multiple samples. "
						"Provide sample_name for the loom file."
					)
			ratio_df[sample_col] = sample_name
			merge_cols = [
				sample_col,
				"_barcode_key",
			]
		elif sample_name is not None:
			ratio_df[sample_col] = sample_name
		summary["_barcode_key"] = _barcode_key(
			summary[barcode_col]
		)
		# Removing -1 or x can theoretically collapse different GEM
		# groups. Refuse to guess silently if that happens.
		duplicate_mask = ratio_df.duplicated(
			merge_cols,
			keep=False,
		)
		if duplicate_mask.any():
			duplicates = (
				ratio_df.loc[
					duplicate_mask,
					merge_cols,
				]
				.drop_duplicates()
				.head(10)
				.to_dict("records")
			)
			raise ValueError(
				"Loom barcodes are not unique after normalization. "
				f"Examples: {duplicates}"
			)
		value_cols = [
			*merge_cols,
			"loom_cell_id",
			"spliced_molecules",
			"unspliced_molecules",
			"intron_exon_ratio"
		]
		result = (
			summary
			.merge(
				ratio_df[value_cols],
				on=merge_cols,
				how="left",
				validate="many_to_one",
			)
			.drop(columns="_barcode_key")
		)
		matched = int(
			result["loom_cell_id"]
			.notna()
			.sum()
		)
		print(
			"[add_loom_intron_exon_ratio] "
			f"Matched {matched}/{len(result)} barcode rows "
			f"from {loom_path.name}"
		)
	# --------------------------------------------------------------
	# Optional output
	# --------------------------------------------------------------
	if output_file is not None:
		output_path = Path(output_file).resolve()
		output_path.parent.mkdir(
			parents=True,
			exist_ok=True,
		)
		result.to_csv(
			output_path,
			sep="\t",
			index=False,
		)
		print(
			"[add_loom_intron_exon_ratio] "
			f"Wrote: {output_path}"
		)
	return result





# =============================================================================
# Single-junction single-sample runner
# =============================================================================

def is_empty(value='', zero_is_empty=False):
	"""
	Returns True if value is empty/null.
	Returns False if value contains usable data.
	"""
	if value is None:
		return True
	try:
		if pd.isna(value) and not isinstance(
			value, (list, tuple, set, dict, pd.Series, pd.DataFrame)
		):
			return True
	except Exception:
		pass
	if isinstance(value, str):
		return value.strip() == ""
	if isinstance(value, (list, tuple, set, dict)):
		return len(value) == 0
	if isinstance(value, pd.DataFrame):
		return value.shape[0] == 0 or value.shape[1] == 0
	if isinstance(value, pd.Series):
		return len(value) == 0
	if zero_is_empty and isinstance(value, (int, float, np.number)):
		return value == 0
	return False


OUTPUT_COLUMNS_sc = [
	"sample", "gene", "jxn_name", "target_gene_coord",
	"targetJxn_coord", "samtools_jxn_coord", "regtools_jxn_coord", "nt_sequence",
	"targetJxn_read_count", "targetJxn_seq_read_count", "pct_reads_seq",
	"targetJxn_cell_count", "targetJxn_seq_cell_count", "pct_cells_seq",
	"targetJxn_umi_count", "targetJxn_seq_umi_count", "pct_umi_seq",
	"geneJxn_count_mean", "geneJxn_count_median", "geneJxn_count_mode",
	"geneJxn_count_sum", "geneJxn_counts", "geneJxn_transcripts"
]

OUTPUT_COLUMNS_noReg_sc = [
	"sample", "gene", "jxn_name", "target_gene_coord",
	"targetJxn_coord", "samtools_jxn_coord", "nt_sequence",
	"targetJxn_read_count", "targetJxn_seq_read_count", "pct_reads_seq",
	"targetJxn_cell_count", "targetJxn_seq_cell_count", "pct_cells_seq",
	"targetJxn_umi_count", "targetJxn_seq_umi_count", "pct_umi_seq"
]

OUTPUT_COLUMNS_bulk = [
	"sample", "gene", "jxn_name", "target_gene_coord",
	"targetJxn_coord", "samtools_jxn_coord", "regtools_jxn_coord", "nt_sequence",
	"targetJxn_read_count", "targetJxn_seq_read_count", "pct_reads_seq",
	"geneJxn_count_mean", "geneJxn_count_median", "geneJxn_count_mode",
	"geneJxn_count_sum", "geneJxn_counts", "geneJxn_transcripts"
]

OUTPUT_COLUMNS_noReg_bulk = [
	"sample", "gene", "jxn_name", "target_gene_coord",
	"targetJxn_coord", "samtools_jxn_coord", "nt_sequence",
	"targetJxn_read_count", "targetJxn_seq_read_count", "pct_reads_seq"
]




def get_output_columns(mode: str, use_regtools: bool) -> list:
	if mode == "sc" and use_regtools:
		return OUTPUT_COLUMNS_sc
	if mode == "sc" and not use_regtools:
		return OUTPUT_COLUMNS_noReg_sc
	if mode == "bulk" and use_regtools:
		return OUTPUT_COLUMNS_bulk
	if mode == "bulk" and not use_regtools:
		return OUTPUT_COLUMNS_noReg_bulk
	raise ValueError("mode must be 'sc' or 'bulk'")


def write_empty_result(
	sample: str,
	gene: str,
	jxn_name: str,
	gene_region: str,
	jxncoord: str,
	nt_seq: str,
	base: str,
	mode: str,
	use_regtools: bool
) -> pd.DataFrame:
	cols = get_output_columns(mode, use_regtools)
	result = {col: pd.NA for col in cols}
	result.update({
		"sample": sample,
		"gene": gene,
		"jxn_name": jxn_name,
		"target_gene_coord": gene_region,
		"targetJxn_coord": jxncoord,
		"nt_sequence": nt_seq,
	})
	result_df = pd.DataFrame([result], columns=cols)
	result_df.to_csv(
		f"{base}.junctScope.txt",
		sep="\t",
		index=False
	)
	return result_df



step_n = 1

def run_sample_junction(
	sample: str,
	bam_file: str,
	jxn_entry: dict,
	output_dir: str,
	mode: str,
	config_func: dict,
	#qc_step: bool = False,
	gene_region: str = None,
	fasta: str = None,
	gtf: str = None,
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
	#qc_step     : bool == True to run sample level read QC
	threads     : samtools threads
	regtools    : bool == True to run regtools steps
	mode        : 'sc' if single cell data, 'bulk' if bulk RNASeq data

	Returns a dict matching OUTPUT_COLUMNS.
	"""
	use_regtools = regtools and validate_regtools()
	jxn_name = jxn_entry["junction_name"]
	gene     = jxn_entry["gene"]
	jxncoord = jxn_entry["junction"]
	nt_seq   = jxn_entry["nt_seq"]
	step_count = 5
	step_count = step_count+1 if use_regtools else step_count
	#step_count = step_count+1 if qc_step else step_count
	global step_n
	step_n = 1
	if mode != 'sc' and mode != 'bulk':
		print("[run_sample_junction] 'mode' argument not found ...")
		print("[run_sample_junction] mode = 'sc' for single-cell data with barcodes")
		sys.exit("[run_sample_junction] mode = 'bulk' for bulk RNASeq data")
		return
	# Derive NT sequence if not provided
	if nt_seq == "NULL" or not nt_seq:
		nt_seq = get_nt_seq(jxncoord, fasta,
			window=config_func['get_nt_seq'].get('window',5))
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
	extract_jxn_region(bam_file, jxncoord, region_sam,
		buffer=config_func['extract_jxn_region'].get('buffer',200),
		threads=threads) # bug7
	# ------------------------------------------------------------------
	# 2. Parse SAM to DataFrame
	# ------------------------------------------------------------------
	print(f"  [{step_n}/{step_count}] Parsing SAM ...")
	step_n += 1
	parsed = parse_sam_per_junction(region_sam, mode)
	if is_empty(parsed):
		print("[EMPTY] Parsed SAM is empty. Writing NA result.")
		return write_empty_result(
			sample, gene, jxn_name, gene_region, jxncoord,
			nt_seq, base, mode, use_regtools
		)
	parsed.to_csv(f"{base}.region.parsed.tsv", sep="\t", index=False)
	# ------------------------------------------------------------------
	# 3. Filter to junction coordinates + sequence check
	# ------------------------------------------------------------------
	print(f"  [{step_n}/{step_count}] Filtering to junction coordinates ...")
	step_n += 1
	filtered = filter_to_junction(parsed, jxncoord, mode,
		buffer=config_func['filter_to_junction'].get('buffer',2),
		filter_barcodes=config_func['filter_to_junction'].get('filter_barcodes',False))
	filtered = seq_check(filtered, nt_seq,
		filter_rows=config_func['seq_check'].get('filter_rows',False),
		case_sensitive=config_func['seq_check'].get('case_sensitive',False))
	if is_empty(filtered):
		print("[EMPTY] Filtered junction table is empty. Writing NA result.")
		return write_empty_result(
			sample, gene, jxn_name, gene_region, jxncoord,
			nt_seq, base, mode, use_regtools
		)
	filtered.to_csv(f"{base}.region.filtered.tsv", sep="\t", index=False)
	# ------------------------------------------------------------------
	# 4. Summarise — sample level and per-barcode level
	# ------------------------------------------------------------------
	print(f"  [{step_n}/{step_count}] Summarising junctions ...")
	step_n += 1
	sam_summ    = summ_junc(filtered, mode)
	if is_empty(sam_summ):
		print("[EMPTY] Junction summary is empty. Writing NA result.")
		return write_empty_result(
			sample, gene, jxn_name, gene_region, jxncoord,
			nt_seq, base, mode, use_regtools
		)
	sam_summ.to_csv(f"{base}.sam.jxn.summ.tsv", sep="\t", index=False)
	#if len(sam_summ) > 1:
	# now get the 'best' junction no matter what, even if only 1 option and print out result and save to best'
	best = select_best_junction(sam_summ, mode = mode,
		min_reads=config_func['select_best_junction'].get('min_reads',5),
		min_seq_reads=config_func['select_best_junction'].get('min_seq_reads',1))
	if best.empty:
		sam_jxn_coord = "NA"
		read_count     = 0
		seq_read_count = 0
		pct_reads_seq  = 0.0
	else:
		top_bulk = best.iloc[0]
		best_chr = top_bulk["chromosome"]
		best_str = top_bulk["junction_start"]
		best_end = top_bulk["junction_end"]
		sam_jxn_coord = f"{best_chr}:{int(best_str)}-{int(best_end)}"
		read_count     = int(top_bulk["sam_read_count"])
		seq_read_count = int(top_bulk["seq_read_count"])
		pct_reads_seq  = round(float(top_bulk["pct_reads_seq"]), 4)
	best.to_csv(f"{base}.sam.jxn.best.summ.tsv", sep="\t", index=False)
	if mode == 'sc':
		bcd_summ    = summ_junc_bybcd(filtered)
		bcd_summ.to_csv(f"{base}.sam.jxn.bcd.summ.tsv", sep="\t", index=False)
		# Best junction from SAM method
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
	# QC --------------------
	# if qc_step:
	# 	print(f"  [{step_n}/{step_count}] Sample QC {jxncoord} ...")
	# 	step_n += 1
	# 	qc_dict = count_sample_stats(
	# 		input_file = bam_file,
	# 		mode = mode,
	# 		verbose = True,
	# 		threads = threads,
	# 		reference_fasta = fasta
	# 	)
	# 	sample_junc_qc(best,qc_dict,f"{base}.jxn.qc.tsv",mode)
	# ------------------------------------------------------------------
	# 5. regtools extract + annotate (gene-level)
	# ------------------------------------------------------------------
	if regtools and validate_regtools():
		print(f"  [{step_n}/{step_count}] Running regtools ...")
		step_n += 1
		reg_bed  = f"{base}.gene.bed"
		reg_anno = f"{base}.gene.bed.anno"
		extract_jxn_regtools(bam_file, gene_region, reg_bed,
			buffer=config_func['extract_jxn_regtools'].get('buffer',1000))
		annotate_jxn_regtools(reg_bed, fasta, gtf, reg_anno)
		reg_df = read_anno_as_df(reg_anno)
		# Match SAM summary to regtools annotation
		reg_jxn_coord = "NA"
		if not sam_summ.empty and not reg_df.empty:
			merged = match_sam_to_regtools(sam_summ, reg_df,
				coord_tol=config_func['match_sam_to_regtools'].get('coord_tol',1))
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
		"pct_reads_seq":            pct_reads_seq
	}
	if mode == 'sc':
		base_result.update({
		"targetJxn_cell_count":     cell_count,
		"targetJxn_seq_cell_count": seq_cell_count,
		"pct_cells_seq":            pct_cells_seq,
		"targetJxn_umi_count":      umi_count,
		"targetJxn_seq_umi_count":  seq_umi_count,
		"pct_umi_seq":              pct_umi_seq})
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
		if matched.empty and regtools and validate_regtools():
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
	else:
		# regtools unavailable — single row, no gene-level stats
		results = [{**base_result}]
	result_df = pd.DataFrame(results)
	output_cols = get_output_columns(mode, use_regtools)
	ordered_cols = [c for c in output_cols if c in result_df.columns]
	#ordered_cols = [c for c in OUTPUT_COLUMNS_sc if c in result_df.columns]
	extra_cols = [c for c in result_df.columns if c not in ordered_cols]
	result_df[ordered_cols].to_csv(
		f"{base}.junctScope.txt",
		sep="\t",
		index=False
	)
	return result_df


# =============================================================================
# Per-barcode merge across junctions and samples
# =============================================================================



def merge_sample_best_jxn(
	sample_dir: str,
	output_file: str = None,
	mode: str = 'sc',
	use_regtools: bool = False
) -> pd.DataFrame:
	sample_dir = str(Path(sample_dir).resolve())
	sample_name = Path(sample_dir).name
	sampSumm_header = "\t".join(get_output_columns(mode,use_regtools))
	if output_file is None:
		base_name = Path(sample_dir).name
		output_file = os.path.join(sample_dir, f"{base_name}.junctScope_SampleSummary.txt")
	else:
		output_file = str(Path(output_file).resolve())
	pattern = str(Path(sample_dir) / "*" /f"*.*.junctScope.txt")
	files = sorted(glob_files(pattern))
	with open(output_file, "w") as out:
		out.write(sampSumm_header + "\n")
		for fpath in files:
			with open(fpath) as fh:
				_hdr = fh.readline()
				out.write(fh.read())
	print(f"Sample summary for {sample_name} written to {output_file}")



def merge_sample_best_bcd_jxn(
	sample_dir: str,
	output_file: str = None,
	loom_file: str = None
) -> pd.DataFrame:
	sample_dir = str(Path(sample_dir).resolve())
	sample_name = Path(sample_dir).name
	if output_file is None:
		base_name = Path(sample_dir).name
		output_file = os.path.join(sample_dir, f"{base_name}_fullJxn_bcd_summary.tsv")
	else:
		output_file = str(Path(output_file).resolve())
	pattern = str(Path(sample_dir) / "*" /f"*.*.sam.jxn.bcd.best.tsv")
	files = sorted(glob_files(pattern))
	per_jxn_dfs = []
	ID_COLS = {"cell_barcode","sample_name"}
	for fpath in files:
		fpath = str(Path(fpath).resolve())
		#df = pd.read_csv(fpath, sep="\t") ####
		df = pd.read_csv(fpath, sep="\t", dtype={"cell_barcode": "string"})
		if df.empty:
			continue
		df["cell_barcode"] = df["cell_barcode"].fillna("NO_BARCODE").astype(str) ####
		jxn_name = Path(fpath).parts[-2]
		df.insert(loc = 4, column = 'sample_name', value = sample_name)
		df = df.drop(columns=['chromosome','junction_start','junction_end'])
		df.columns = df.columns.str.replace(r'_per_barcode|_bcd', '', regex=True)
		rename = {
			c: f"{jxn_name}_{c}"
			for c in df.columns
			if c not in ID_COLS
		}
		df = df.rename(columns=rename)
		per_jxn_dfs.append(df)
	if not per_jxn_dfs:
		print("[merge_sample_best_bcd_jxn] No data found — check that sample has been run "
			  f"and that sample_dir is correct:\n  {sample_dir}")
		return pd.DataFrame()
	#sample_jxn_summ = reduce(lambda left, right: pd.merge(left, right, how='outer'), per_jxn_dfs)
	sample_jxn_summ = reduce(
		lambda left, right: pd.merge(
			left,
			right,
			on=["cell_barcode", "sample_name"],
			how="outer"
		),
		per_jxn_dfs
	)
	numeric_cols = sample_jxn_summ.select_dtypes(include=[np.number]).columns
	sample_jxn_summ[numeric_cols] = sample_jxn_summ[numeric_cols].fillna(0)
	if loom_file is not None:
		sample_jxn_summ = add_loom_intron_exon_ratio(
			loom_file=loom_file,
			barcode_summary=sample_jxn_summ,
			sample_name=sample_name,
		)
	sample_jxn_summ.to_csv(output_file, sep="\t", index=False)
	print(f"[merge_sample_best_bcd_jxn] {sample_name}: {len(sample_jxn_summ)} barcode rows "
		f"from {len(per_jxn_dfs)} junctions(s) → {output_file}")
	return sample_jxn_summ ####


def merge_fullJxn_bcd_best(
	output_dir: str,
	output_file: str = None,
) -> pd.DataFrame:
	"""
	Row-bind sample-level full-junction barcode summaries.
	Samples may contain different junction columns. Missing columns are added
	and filled with NA in the combined output.
	"""
	output_dir = str(Path(output_dir).resolve())
	base_name = Path(output_dir).name
	if output_file is None:
		output_file = os.path.join(
			output_dir,
			f"{base_name}_fullJxn_bcd_summary.tsv"
		)
	else:
		output_file = str(Path(output_file).resolve())
	pattern = str(Path(output_dir) / "Intermediate" / "*" / "*_fullJxn_bcd_summary.tsv")
	files = sorted(glob_files(pattern))
	if not files:
		print(
			"[merge_fullJxn_bcd_best] No sample barcode summary files found.\n"
			f"  Pattern: {pattern}"
		)
		return pd.DataFrame()
	sample_bcd_summ_dfs = []
	for fpath in files:
		fpath = str(Path(fpath).resolve())
		df = pd.read_csv(
			fpath,
			sep="\t",
			dtype={
				"cell_barcode": "string",
				"sample_name": "string",
			}
		)
		if df.empty:
			print(
				f"[merge_fullJxn_bcd_best] WARNING: empty file skipped: "
				f"{fpath}"
			)
			continue
		# Retain missing barcodes as an explicit label if desired.
		df["cell_barcode"] = df["cell_barcode"].fillna("NO_BARCODE")
		sample_bcd_summ_dfs.append(df)
	if not sample_bcd_summ_dfs:
		print(
			"[merge_fullJxn_bcd_best] Files were found, but all were empty."
		)
		return pd.DataFrame()
	# Preserve first-seen column order across all sample DataFrames.
	column_order = list(dict.fromkeys(
		column
		for df in sample_bcd_summ_dfs
		for column in df.columns
	))
	# Reindex adds absent columns and fills them with NA.
	aligned_dfs = [
		df.reindex(columns=column_order)
		for df in sample_bcd_summ_dfs
	]
	# Row-bind samples. Missing junction columns remain NA.
	merged = pd.concat(
		aligned_dfs,
		axis=0,
		ignore_index=True,
		sort=False,
	)
	merged.to_csv(
		output_file,
		sep="\t",
		index=False,
	)
	print(
		f"[merge_fullJxn_bcd_best] Merged "
		f"{len(sample_bcd_summ_dfs)} samples into "
		f"{len(merged)} barcode rows × {len(merged.columns)} columns "
		f"→ {output_file}"
	)
	return merged

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
			ACGTACGT-1  →  ACGTACGT-1_<sample>
	  - numeric count/pct columns are prefixed with the junction name:
			reads_per_barcode  →  <junction_name>_reads_per_barcode

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
	ID_COLS = {"cell_barcode","sample_name"}
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
			#df = pd.read_csv(fpath, sep="\t")
			df = pd.read_csv(fpath, sep="\t", dtype={"cell_barcode": "string"}) ####
			df["cell_barcode"] = df["cell_barcode"].fillna("NO_BARCODE").astype(str)
			if df.empty:
				continue
			# Derive sample name from directory structure:
			# Intermediate/<sample>/<jxn_name>/<file>
			sample_name = Path(fpath).parts[-3]
			# Suffix cell_barcode with sample name
			df.insert(loc = 4, column = 'sample_name', value = sample_name)
			# df["cell_barcode"] = (
			# 	df["cell_barcode"].astype(str) + "_" + sample_name
			# )
			df = df.drop(columns=['chromosome','junction_start','junction_end'])
			df.columns = df.columns.str.replace(r'_per_barcode|_bcd', '', regex=True)
			# Prefix all value columns with junction name
			rename = {
				c: f"{jxn_name}_{c}"
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
		# Sort by read count descending
		sort_col = f"{jxn_name}_reads"
		if sort_col in jxn_long.columns:
			jxn_long = jxn_long.sort_values(sort_col, ascending=False)
		jxn_long = jxn_long.reset_index(drop=True)
		# Fill numeric NAs with 0
		numeric_cols = jxn_long.select_dtypes(include=[np.number]).columns
		jxn_long[numeric_cols] = jxn_long[numeric_cols].fillna(0)
		# Write per-junction barcode summary file
		jxn_out = os.path.join(output_dir, f"{base_name}_{jxn_name}_bcd_summary.tsv")
		jxn_long.to_csv(jxn_out, sep="\t", index=False)
		jxn_frames[jxn_name] = jxn_long
		print(f"[merge_bcd_best] {jxn_name}: {len(jxn_long)} barcode rows "
			  f"from {len(per_sample_dfs)} sample(s) → {jxn_out}")
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
			merged = merged.merge(df, on=["cell_barcode","sample_name"], how="outer")
	# Sort by first junction's read count descending (if present)
	first_jxn = jxn_list[0]["junction_name"]
	sort_col = f"{first_jxn}_reads_per_barcode"
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


# =============================================================================
# Setup  — derive gene regions, write per-sample execution scripts
# =============================================================================

def _setup_references(config: dict) -> tuple[list[dict], dict, dict]:
	"""
	Minimal setup shared by both full-pipeline and --sample-only modes.

	Resolves FASTA / GTF / OUTPUT to absolute paths, builds the junction list,
	subsets the GTF per gene, derives gene regions, and fills any missing NT
	sequences.  Does NOT read INPUT or write any execution scripts.

	Stores derived values on config under the same keys that setup() uses
	(_gene_gtfs, _gene_regions, _jxn_list) so _run_one() works identically
	in both modes.

	Returns (jxn_list, gene_regions, gene_gtfs).
	"""
	# Resolve shared paths to absolute
	config["proj_name"] = str(Path(config["proj_name"]).resolve())
	if config.get("fasta", ""):
		config["fasta"] = str(Path(config["fasta"]).resolve())
	if config.get("gtf", ""):
		config["gtf"]   = str(Path(config["gtf"]).resolve())
	regtools         = bool(config["regtools"])
	output = config["proj_name"]
	fasta  = config["fasta"]
	gtf    = config["gtf"]
	os.makedirs(output, exist_ok=True)
	jxn_list = config_to_jxn_list(config)
	gene_gtfs    = {}
	gene_regions = {}
	if regtools: # bug8
		# Per-gene GTF subsets and gene regions
		genes = list(dict.fromkeys(j["gene"] for j in jxn_list))
		for gene in genes:
			gene_gtf_path = str(Path(output) / f"temp_{gene}.gtf")
			subset_gtf(gtf, gene, gene_gtf_path)
			gene_gtfs[gene] = gene_gtf_path
			chrom  = jxn_list[0]["junction"].split(":")[0]
			region = gene_region_from_gtf(gene_gtf_path, chrom)
			if not region:
				sys.exit(f"[_setup_references] No GTF entries found for gene '{gene}'")
			gene_regions[gene] = region
			print(f"[setup] {gene} region: {region}")
	# Store on config so callers that read config directly get the derived values
	config["_gene_gtfs"]    = gene_gtfs
	config["_gene_regions"] = gene_regions
	# Fill missing NT sequences
	for jxn in jxn_list: # bug9
		if not jxn["nt_seq"] or jxn["nt_seq"] == "NULL":
			if is_empty(fasta):
				sys.exit("[setup] FASTA file required to derive NT sequence (ntseq not provided)")
			jxn["nt_seq"] = get_nt_seq(jxn["junction"], fasta, window=5)
			print(f"[setup] NT seq for {jxn['junction_name']}: {jxn['nt_seq']}")
	config["_jxn_list"]     = jxn_list
	return jxn_list, gene_regions, gene_gtfs




def write_slurm_array_script(
	project_dir: str,
	sample_jobs: list[tuple[str, str]],
	threads: int = 8,
	mem_per_cpu: str = "64G",
	time_limit: str = "1-00:00:00",
	anaconda_module: str = "Anaconda3/2022.05",
	conda_env: str = "juncScope_velocyto_scvelo",
) -> Path:
	"""
	Write an optional Slurm array submission script for the generated
	per-sample junctionScope runners.

	Parameters
	----------
	project_dir : str
		Top-level junctionScope project directory.
	sample_jobs : list of (sample_name, runner_path) tuples
		One entry per sample that has a generated execution script.
	threads : int
		Number of CPUs requested for each array task. This is written to
		``#SBATCH --cpus-per-task`` and should match the pipeline thread count.
	mem_per_cpu : str
		Slurm memory request per CPU. Default is 32G.
	time_limit : str
		Slurm walltime in D-HH:MM:SS format. Default is three days.
	anaconda_module : str
		Environment module loaded at the beginning of each task.
	conda_env : str
		Conda environment activated before running the sample script.

	Returns
	-------
	pathlib.Path
		Path to the generated Slurm shell script.
	"""
	project_dir = Path(project_dir).expanduser().resolve()
	project_dir.mkdir(parents=True, exist_ok=True)

	if not sample_jobs:
		raise ValueError(
			"Cannot create a Slurm array script because no sample execution "
			"scripts were generated. Check the input BAM/CRAM paths."
		)
	if int(threads) < 1:
		raise ValueError("threads must be at least 1")

	# Resolve paths now so the submitted job does not depend on its launch cwd.
	resolved_jobs = [
		(str(sample), str(Path(runner).expanduser().resolve()))
		for sample, runner in sample_jobs
	]

	project_name = project_dir.name
	job_stub = re.sub(r"[^A-Za-z0-9_.-]+", "_", project_name).strip("_")
	job_name = f"jScope_{job_stub or 'project'}"
	slurm_script = project_dir / f"{project_name}_slurm_array.sh"
	log_dir = project_dir / "slurm_logs"
	log_dir.mkdir(parents=True, exist_ok=True)

	sample_names = "\n".join(
		f"\t{shlex.quote(sample)}"
		for sample, _ in resolved_jobs
	)
	runner_paths = "\n".join(
		f"\t{shlex.quote(runner)}"
		for _, runner in resolved_jobs
	)
	array_end = len(resolved_jobs) - 1

	script_text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --cpus-per-task={int(threads)}
#SBATCH --mem={mem_per_cpu}
#SBATCH --time={time_limit}
#SBATCH --array=0-{array_end}%5
#SBATCH --output={log_dir}/%x_%A_%a.out
#SBATCH --error={log_dir}/%x_%A_%a.err
#SBATCH --export=NONE

set -eo pipefail

module purge
ml {shlex.quote(anaconda_module)}

eval "$(conda shell.bash hook)"
conda activate {shlex.quote(conda_env)}

SAMPLE_NAMES=(
{sample_names}
)

SAMPLE_SCRIPTS=(
{runner_paths}
)

TASK_ID="${{SLURM_ARRAY_TASK_ID}}"
if (( TASK_ID < 0 || TASK_ID >= ${{#SAMPLE_SCRIPTS[@]}} )); then
	echo "Invalid SLURM_ARRAY_TASK_ID: $TASK_ID" >&2
	exit 1
fi

SAMPLE="${{SAMPLE_NAMES[$TASK_ID]}}"
RUNNER="${{SAMPLE_SCRIPTS[$TASK_ID]}}"

echo "[$(date)] Starting junctionScope sample: $SAMPLE"
echo "[$(date)] Runner: $RUNNER"
echo "[$(date)] CPUs: $SLURM_CPUS_PER_TASK"

python3 "$RUNNER"

echo "[$(date)] Completed junctionScope sample: $SAMPLE"
"""

	slurm_script.write_text(script_text)
	os.chmod(slurm_script, 0o750)
	return slurm_script


def setup(config: dict) -> tuple[list[dict], list[tuple[str, str, str]]]:
	"""
	Prepare project directory, subset GTFs, derive gene regions, derive NT
	sequences where missing, and write per-sample-per-junction execution scripts.

	Returns (jxn_list, samples) so callers can iterate directly.
	"""
	# Resolve every path to absolute so generated scripts and file I/O work
	# regardless of which directory the caller is in.
	config_main = config['main_config']
	config_func = config['function_config']
	jxn_list, gene_regions, gene_gtfs = _setup_references(config_main)
	output     = config_main["proj_name"]
	input_list = str(Path(config_main["input"]).resolve())
	config_main["input"] = input_list
	buffer     = int(config_main["buffer"])
	threads    = int(config_main["threads"])
	qc_step    = bool(config_main["qc_step"])
	velocyto    = bool(config_main["velocyto"])
	regtools   = bool(config_main["regtools"])
	mode       = str(config_main["mode"])
	fasta      = config_main["fasta"]
	gtf        = config_main["gtf"]
	# ── Parse sample list ────────────────────────────────────────────────────
	samples = []
	looms = {}
	#if input_list != None:
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
			loom_found = [
				value.strip().strip('"')
				for value in parts
				if re.search(
					r"\.loom$",
					value.strip().strip('"'),
					flags=re.IGNORECASE,
				)
			]
			#loom_found = [s for s in parts if re.search(".loom$", s)]
			if loom_found:
				loom_file = str(loom_found[0].strip('"'))
				looms[sample.replace(".", "")] = str(loom_file)
			samples.append((sample.replace(".", ""), bam_file))
	# ── Per-sample execution scripts ─────────────────────────────────────────
	output_abs  = str(Path(output).resolve())   # already absolute; defensive copy
	master_exc  = str(Path(output_abs) / f"{Path(output_abs).name}_exc.sh")
	script_path = str(Path(__file__).resolve())
	sample_jobs = []
	with open(master_exc, "w") as master:
		for sample, bam_file in samples:
			bam_abs = str(Path(bam_file).resolve())
			loom_abs = looms.get(sample,None)
			if not os.path.isfile(bam_abs):
				print(f"[setup] WARNING: {bam_abs} not found — skipping {sample}")
				continue
			sample_dir = str(Path(output) / "Intermediate" / sample)
			os.makedirs(sample_dir, exist_ok=True)
			exc_script = str(Path(sample_dir) / f"{sample}_exc.py")
			with open(exc_script, "w") as exc:
				exc.write(_render_exc_script(
					script_path, sample, bam_abs, loom_abs,
					jxn_list, gene_regions, gene_gtfs,
					fasta, gtf, buffer, threads, output_abs,
					regtools, qc_step, velocyto, mode, config_func
				))
			# Write absolute path so the master script runs from any directory
			master.write(f"python3 {shlex.quote(exc_script)}\n")
			sample_jobs.append((sample, exc_script))
	print(f"[setup] Master execution script:\n\tsh {master_exc}\n")

	# ── Optional Slurm array script ───────────────────────────────────────────
	if sample_jobs:
		slurm_script = write_slurm_array_script(
			project_dir=output_abs,
			sample_jobs=sample_jobs,
			threads=threads,
		)
		print(f"[setup] Optional Slurm array job:\n\tsbatch {slurm_script}\n")
	else:
		print("[setup] WARNING: no valid sample runners were generated; "
			  "Slurm array script was not written")
	# ── Summary script ───────────────────────────────────────────────────────
	summary_script = str(Path(output_abs) / "junctScopeSummarize.py")
	with open(summary_script, "w") as summ:
		summ.write(_render_summary_script(output_abs, jxn_list, regtools, mode, qc_step))
	print(f"[setup] After all samples finish:\n\tpython3 {summary_script}\n")
	return jxn_list, samples, looms


def _render_exc_script(
	script_path, sample, bam_file, loom_file, jxn_list,
	gene_regions, gene_gtfs, fasta, gtf, buffer,
	threads, output, regtools, qc_step, velocyto, mode, config_func
) -> str:
	"""Render the Python source for a per-sample execution script."""
	merge_block = (
		f"""
from junctionScope import merge_sample_best_bcd_jxn
merge_sample_best_bcd_jxn(os.path.join({repr(output)}, "Intermediate", sample), loom_file = loom_file)
"""
		if mode == "sc"
		else ""
	)
	# if loom_file:
	# 	print(f"[107] loom_file: {loom_file!r}")
	# 	print(f"[107] loom exists: {Path(loom_file).exists() if loom_file else False}")
	# 	print(f"[107] loom size: {Path(loom_file).stat().st_size if loom_file and Path(loom_file).exists() else 'NA'}")
	return f"""#!/usr/bin/env python3
\"\"\"Auto-generated runner for sample {sample}.\"\"\"
import sys, os
sys.path.insert(0, {repr(str(Path(script_path).parent))})
from junctionScope import run_sample_junction
from junctionScope import sample_junc_qc
from junctionScope import count_sample_stats
from junctionScope import merge_sample_best_jxn
from junctionScope import merge_qc_within_sample
from junctionScope import *
from pathlib import Path
import scanpy as sc
import scvelo as scv
import cellrank as cr
import anndata as ad
import pyranges as pr

os.chdir({repr(str(Path(output).resolve()))})

sample   = {repr(sample)}
bam_file = {repr(bam_file)}
loom_file= {repr(loom_file)}
fasta    = {repr(fasta)}
buffer   = {buffer}
threads  = {threads}
regtools = {regtools}
qc_step  = {qc_step}
velocyto = {velocyto}
mode     = {repr(mode)}
config_func = {repr(config_func)}

jxn_list     = {repr(jxn_list)}
gene_regions = {repr(gene_regions)}
gene_gtfs    = {repr(gene_gtfs)}
gtf          = {repr(gtf)}

if velocyto:
	if not gtf:
		print(f"[{{sample}}] WARNING: gtf not provided — skipping velocyto")
	elif not validate_velocyto(verbose=True):
		print(f"[{{sample}}] WARNING: velocyto not available — skipping loom generation")
	else:
		print(f"[{{sample}}] Running velocyto ...")
		loom_file = bam_to_loom(
			bam_file=bam_file,
			gtf_ref=gtf,
			out_dir=os.path.join({repr(output)}, "Intermediate", sample),
			threads=threads
		)
		print(f"[{sample}] Velocyto loom: {loom_file}")

if qc_step == True:
	print(f"[{{sample}}] Sample QC Processing ...")
	qc_dict = count_sample_stats(
		input_file      = bam_file,
		mode            = mode,
		verbose         = True,
		threads         = threads,
		reference_fasta = fasta
	)

for jxn in jxn_list:
	gene = jxn["gene"]
	print(f"[{{sample}}] junction: {{jxn['junction_name']}}")
	jxn_res = run_sample_junction(
		sample      = sample,
		bam_file    = bam_file,
		jxn_entry   = jxn,
		gene_region = gene_regions.get(gene, None),
		fasta       = fasta,
		gtf         = gene_gtfs.get(gene, None),
		output_dir  = os.path.join({repr(output)}, "Intermediate", sample, jxn["junction_name"]),
		buffer      = buffer,
		threads     = threads,
		regtools    = regtools,
		mode        = mode,
		config_func = config_func
	)
	if qc_step == True:
		qc_file_name = '.'.join([sample, jxn["junction_name"], "jxn.qc.tsv"])
		qc_file = os.path.join({repr(output)}, "Intermediate", sample, jxn["junction_name"], qc_file_name)
		print(f"[{{sample}}] junction: {{jxn['junction_name']}}: Merge QC ")
		sample_junc_qc(jxn_res,qc_dict,
			mode = mode,
			output = qc_file)

merge_sample_best_jxn(os.path.join({repr(output)}, "Intermediate", sample), mode = mode, use_regtools = regtools)
if qc_step:
	merge_qc_within_sample(os.path.join({repr(output)}, "Intermediate"), sample, os.path.join({repr(output)}, "Intermediate", sample, "summaries"))

{merge_block}

print(f"[{{sample}}] all junctions complete.")
"""


def _render_summary_script(output: str, jxn_list: list, regtools: bool = False, mode: str = 'sc', qc_step: bool = False) -> str:
	"""Render the Python source for the cross-sample summary script."""
	# if regtools and validate_regtools():
	# 	if mode == 'sc':
	# 		header = "\t".join(OUTPUT_COLUMNS_sc)
	# 	else:
	# 		header = "\t".join(OUTPUT_COLUMNS_bulk)
	# else:
	# 	if mode == 'sc':
	# 		header = "\t".join(OUTPUT_COLUMNS_noReg_sc)
	# 	else:
	# 		header = "\t".join(OUTPUT_COLUMNS_noReg_bulk)
	header = "\t".join(get_output_columns(mode, regtools))
	output_abs = str(Path(output).resolve())
	qc_block = (
		f"""
	merge_qc_across_samples(os.path.join(OUTPUT_DIR, "Intermediate"), jxn_name, os.path.join(OUTPUT_DIR, "qc_summaries"))
"""
		if qc_step == True
		else ""
	)
	merge_block = (
		"""
# ── Per-barcode merge ────────────────────────────────────────────────────────
from junctionScope import merge_bcd_best
merge_fullJxn_bcd_best(OUTPUT_DIR)
"""
		if mode == "sc"
		else ""
	)
	return f"""#!/usr/bin/env python3
\"\"\"Collect per-sample per-junction results into a single summary table
and merge per-barcode best-junction files across junctions and samples.\"\"\"
import glob, os, sys
import numpy as np
from pathlib import Path

# Absolute paths baked in at setup time — script runs correctly from any directory
OUTPUT_DIR  = {repr(output_abs)}
header      = {repr(header)}
jxn_list    = {repr(jxn_list)}

sys.path.insert(0, {repr(str(Path(__file__).resolve().parent))})
from junctionScope import *

# ── Sample-level summary ────────────────────────────────────────────────────
for jxn in jxn_list:
	jxn_name = jxn["junction_name"]
	pattern = os.path.join(OUTPUT_DIR, "Intermediate", "*", jxn_name, "*.junctScope.txt")
	output_file = os.path.join(OUTPUT_DIR, Path(OUTPUT_DIR).name + "_" + jxn_name + "_SampleSummary.txt")
	with open(output_file, "w") as out:
		out.write(header + "\\n")
		for txt in sorted(glob.glob(pattern, recursive=True)):
			with open(txt) as fh:
				# Skip the header row written by result_df.to_csv
				_hdr = fh.readline()
				out.write(fh.read())
	{qc_block}
	print(f"Sample summary for {{jxn_name}} written to {{output_file}}")

{merge_block}
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
  # Setup project (write per-sample and optional Slurm scripts, do not run)
  python junctionScope.py -c project.conf --setup-only
  sbatch <project>/<project>_slurm_array.sh

  # Setup + run all samples sequentially
  python junctionScope.py -c project.conf

  # Run a single sample (e.g. from a Slurm array)
  python junctionScope.py -c project.conf --sample SampleA --bam /path/to/A.cram
		""",
	)
	parser.add_argument("-c", "--config", help="Path to .conf/.json file")
	parser.add_argument("-o", "--output", default = "junctionScope_Output", help="Project name for output folder")
	parser.add_argument("-j", "--jxn-table", help="Path to junction table")
	parser.add_argument("-i", "--input", help="Path bam/cram file list")
	parser.add_argument("-f", "--fasta", help="Path reference FASTA (required with --run-regtools if not in config)")
	parser.add_argument("-g", "--gtf", help="Path to reference GTF (required with --run-regtools if not in config)")
	parser.add_argument("-m","--mode", help="Input data sequencing type 'sc' for single-cell or 'bulk' for bulk RNASeq")
	parser.add_argument("-q","--qc", action="store_true", default = False, help="Perform QC step during sample processing (will add processing time)")
	parser.add_argument("-v","--velocyto", action="store_true", default = False, help="Generate loom file via velocyto (velocyto install required + will add processing time)")
	parser.add_argument("-l", "--loom", default = None, help="Path to loom file if running single sample")
	parser.add_argument("--setup-only", action="store_true", help="Write execution scripts but do not run samples")
	parser.add_argument("--run-regtools", action="store_true", default = False, help="Run regtools junction extraction and annotation steps")
	parser.add_argument("--sample", help="Run a single named sample")
	parser.add_argument("--bam", help="BAM/CRAM path (required with --sample)")
	args = parser.parse_args()

	# ── Build config_main, either from a config file or from CLI args alone ──
	if args.config:
		config = parse_config(args.config)
		config_main = config["main_config"]
		# Individual CLI flags override the config file when both are given
		if args.output is not None and args.output != "junctionScope_Output":
			config_main["proj_name"] = args.output
		if args.jxn_table:
			config_main["jxn_table"] = args.jxn_table
		if args.input:
			config_main["input"] = args.input
		if args.fasta:
			config_main["fasta"] = args.fasta
		if args.gtf:
			config_main["gtf"] = args.gtf
		if args.qc:
			config_main["qc_step"] = args.qc
		if args.velocyto:
			config_main["velocyto"] = args.velocyto
	else:
		# No config file — config_main is built entirely from CLI args.
		# jxn_table, input, and mode are mandatory in this path.
		if args.sample:
			missing = [
			name for name, val in
			[("--jxn-table", args.jxn_table), ("--mode", args.mode)]
			if not val
			]
		else:
			missing = [
				name for name, val in
				[("--jxn-table", args.jxn_table), ("--input", args.input), ("--mode", args.mode)]
				if not val
			]
		if missing:
			sys.exit(
				"[junctionScope] No --config provided — the following are required "
				f"on the command line: {', '.join(missing)}"
			)
		config_main = copy.deepcopy(config_default["main_config"])
		config_main["proj_name"]  = args.output
		config_main["jxn_table"]  = args.jxn_table
		config_main["input"]      = args.input
		config_main["mode"]       = args.mode
		config_main["qc_step"]    = args.qc
		config_main["velocyto"]   = args.velocyto
		if args.fasta:
			config_main["fasta"] = args.fasta
		if args.gtf:
			config_main["gtf"] = args.gtf
		config = {
			"main_config": config_main,
			"function_config": copy.deepcopy(config_default["function_config"]),
		}

	#config = parse_config(args.config)
	#config_main = config['main_config']

	config_map = {
		"sc": ["singlecell", "single-cell", "single cell","sc"],
		"bulk": ["rna", "rnaseq","rna-seq","rna seq","bulk","bulk rna"]
	}
	if args.mode is not None:
		(mode_cli,) = [k for k, v in config_map.items() if args.mode in v]
		config_main["mode"] = mode_cli
		config['main_config']["mode"] = mode_cli
	if config_main["mode"] not in config_map.keys():
		(mode_conf,) = [k for k, v in config_map.items() if config_main["mode"] in v]
		config_main["mode"] = mode_conf
		config['main_config']["mode"] = mode_conf
	mode = config_main["mode"]

	if config_main['regtools'] == False and args.run_regtools:
		config_main["regtools"] = args.run_regtools
	regtools = config_main["regtools"]

	# ── Shared helper: run all junctions for one sample ─────────────────────
	def _run_one(sample, bam_file, jxn_list, loom_file):
		bam_abs    = str(Path(bam_file).resolve())
		loom_abs    = str(Path(loom_file).resolve())
		fasta   = config_main["fasta"]
		buffer  = int(config_main["buffer"])
		threads = int(config_main["threads"])
		proj_name  = config_main["proj_name"]
		qc_step    = bool(config_main["qc_step"])
		velocyto    = bool(config_main["velocyto"])
		gene_regions = config_main.get("_gene_regions", {})
		gene_gtfs    = config_main.get("_gene_gtfs", {})
		if qc_step == True:
			print(f"[{sample}] Sample QC Processing ...")
			qc_dict = count_sample_stats(
				input_file      = bam_abs,
				mode            = mode,
				verbose         = True,
				threads         = threads,
				reference_fasta = fasta
			)
		# Bug D: validate velocyto and guard gtf before calling
		if velocyto:
			gtf_path = config_main.get("gtf", "")
			if not gtf_path:
				print(f"[{sample}] WARNING: gtf not provided — skipping velocyto")
			elif not validate_velocyto(verbose=True):
				print(f"[{sample}] WARNING: velocyto not available — skipping loom generation")
			else:
				print(f"[{sample}] Running velocyto ...")
				bam_to_loom(bam_abs, gtf_path,
					os.path.join(proj_name, "Intermediate", sample), threads)
		for jxn in jxn_list:
			sample_dir = os.path.join(proj_name, "Intermediate", sample, jxn['junction_name'])
			gene = jxn["gene"]
			print(f"[{sample}] junction: {jxn['junction_name']}")
			jxn_res = run_sample_junction(
				sample      = sample,
				bam_file    = bam_abs,
				jxn_entry   = jxn,
				gene_region = gene_regions.get(gene, None),
				fasta       = fasta,
				gtf         = gene_gtfs.get(gene, None),
				output_dir  = sample_dir,
				buffer      = buffer,
				threads     = threads,
				regtools    = regtools,
				mode        = mode,
				config_func = config["function_config"]
			)
			if qc_step == True:
				qc_file_name = '.'.join([sample, jxn["junction_name"], "jxn.qc.tsv"])
				qc_file = os.path.join(proj_name, "Intermediate", sample, jxn['junction_name'], qc_file_name)
				print(f"[{sample}] junction: {jxn['junction_name']}: Merge QC")
				sample_junc_qc(jxn_res, qc_dict,
					mode=mode, output=qc_file)
		# Bug H: call both merge functions (generated script has both)
		merge_sample_best_jxn(os.path.join(proj_name, "Intermediate", sample),
			mode=mode, use_regtools=regtools)
		# Bug K: bcd merge only for sc mode
		if mode == "sc":
			# print(f"[107] loom_file: {loom_abs!r}")
			# print(f"[107] loom exists: {Path(loom_abs).exists() if loom_abs else False}")
			# print(f"[107] loom size: {Path(loom_abs).stat().st_size if loom_abs and Path(loom_abs).exists() else 'NA'}")
			merge_sample_best_bcd_jxn(os.path.join(proj_name, "Intermediate", sample), loom_file = loom_abs)
		# Bug I: "summaries" was unquoted bare name
		if qc_step:
			merge_qc_within_sample(
				os.path.join(proj_name, "Intermediate"), sample,
				os.path.join(proj_name, "Intermediate", sample, "summaries"))



	# Single sample mode (--sample + --bam) ----------
	if args.sample:
		if not args.bam:
			sys.exit("--sample requires --bam")
		validate_config(config_main, require_input=False)
		jxn_list, gene_regions, gene_gtfs = _setup_references(config_main)
		_run_one(args.sample.replace(".", ""), args.bam, jxn_list, args.loom)
		return

	# Full setup ----------------
	validate_config(config_main, require_input=True)
	jxn_list, samples, looms = setup(config)

	if args.setup_only:
		return
	
	for sample, bam_file in samples:
		if os.path.isfile(bam_file):
			loom_abs = looms.get(sample,None)
			_run_one(sample, bam_file, jxn_list, loom_abs)
		else:
			print(f"[main] WARNING: {bam_file} not found — skipping {sample}")
	result = subprocess.run(
	    [sys.executable, os.path.join(config_main["proj_name"],'junctScopeSummarize.py')], capture_output=True, text=True
	)
	# Access the script's output and errors
	print("Output:", result.stdout)
	if result.stderr:
		print("Errors:", result.stderr)



if __name__ == "__main__":
	main()