#! /bin/env bash

# Example Use
# sh junctionScope_setup.sh -c juncitonScope_Proj.conf

usage() { echo "Usage: $0 [-c <config file>]" 1>&2; exit 1; }


# built and setup with (in order): 
# ml RegTools/0.5.2-foss-2021b
# ml SAMtools/1.9-foss-2018b
# ml OpenBLAS/0.3.18-GCC-11.2.0


# Defaults
BUFFER=200
THREADS=1
STRAND=0
NTSEQ="NULL"

while getopts ":c:" opt; do
  case ${opt} in
  c)
  	CONFIG=${OPTARG}
  	;;
  	\?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

if [ -f ${CONFIG} ]; then
	source ./${CONFIG}
fi

if [ -z "${INPUT}" ] || [ -z "${GENE}" ] || [ -z "${JXNCOORD}" ] || [ -z "${NTSEQ}" ] || [ -z "${FASTA}" ] || [ -z "${GTF}" ]; then
    usage
fi

mkdir -p ${OUTPUT}

# Subset GTF to chr for reduced read time
chr=$(echo ${JXNCOORD} | cut -d ':' -f 1)
grep "\"${GENE}\"" ${GTF} > ${OUTPUT}/temp_gene.gtf
GTF=$(pwd)/${OUTPUT}/temp_gene.gtf


# Get gene coordinate from GTF
GENE_start=$(grep "\"${GENE}\"" ${GTF} | cut -f4 | sort -n | head -n1)
GENE_end=$(grep "\"${GENE}\"" ${GTF} | cut -f5 | sort -rn | head -n1)
GENE_REGION="${chr}:${GENE_start}-${GENE_end}"

# Get directory name where src code is stored
code_dir=$(pwd)/$(dirname ${0})/

# add base start/base end to config file

if [ ${NTSEQ} == "NULL" ]; then
	NTSEQ=$(sh ${code_dir}getNTseq.sh -r ${JXNCOORD} -w 5 -f ${FASTA})
fi

# Get number of columns
inFileCols=$(awk -F'\t' '{print NF; exit}' ${INPUT})




# Set up workflow for each sample of from input file

while IFS=$'\t' read -a line; do
	if [ ${#line[@]} -eq 1 ]; then
		file=${line[0]}
		file="${file//[$'\r\n']/}" # Correct end of line issues
		file="${file//\"/}" # Remove any double quotations
		sample=$(basename ${file%.*}) # Extract sample name from file name
	else
		file=${line[1]}
		file="${file//[$'\r\n']/}"
		file="${file//\"/}"
		sample=${line[0]}
		sample="${sample//\"/}"
	fi
	if [ -f ${file} ]; then
		sample="${sample//./}"
		mkdir -p ${OUTPUT}/Intermediate/${sample}
			cat <<EOF > ${OUTPUT}/Intermediate/${sample}/${sample}_exc.sh
#! /bin/env bash

# Sequence Method target junction counts
sh ${code_dir}extractJxnRegion.sh -i ${file} -r ${JXNCOORD} -o Intermediate/${sample}/${sample}.region.sam -b ${BUFFER} -t ${THREADS}
sh ${code_dir}extractJxnSequence.sh -i Intermediate/${sample}/${sample}.region.sam -n ${NTSEQ} -o Intermediate/${sample}/${sample}.region.seq.sam
jxnCount_seq=\$(sh ${code_dir}countJxnSequence.sh -i Intermediate/${sample}/${sample}.region.seq.sam)

# Coordinate Method gene and target junction counts
sh ${code_dir}extractJxnRegtools.sh -i ${file} -r ${GENE_REGION} -b 1000 -o Intermediate/${sample}/${sample}.gene.bed
sh ${code_dir}annotateJxnRegtools.sh -i Intermediate/${sample}/${sample}.gene.bed -f ${FASTA} -g ${GTF} -o Intermediate/${sample}/${sample}.gene.bed.anno
jxnCount_coord=\$(sh ${code_dir}matchJxnRegtools.sh -i Intermediate/${sample}/${sample}.gene.bed.anno -G ${GENE} -r ${JXNCOORD})
sh ${code_dir}getJxnCountsGene.sh -i Intermediate/${sample}/${sample}.gene.bed.anno -G ${GENE} -o Intermediate/${sample}/${sample}.gene.jxn.counts

# Get max target junction count between seq and coord methods
jxnCount_max=\$(echo "\${jxnCount_seq} \${jxnCount_coord}" | tr ' ' '\n' | sort -rn | head -n 1)

# Get total sample junction counts
jxnCount_samp=\$(sh ${code_dir}getJxnCountsSample.sh -i ${file} -t ${THREADS})

# Write out sample results to file
echo -e "${sample}\t${GENE}\t${JXNCOORD}\t${GENE_REGION}\t${NTSEQ}\t\${jxnCount_seq}\t\${jxnCount_coord}\t\${jxnCount_max}\t\${jxnCount_samp}" | paste -d'\t' - Intermediate/${sample}/${sample}.gene.jxn.counts > Intermediate/${sample}/${sample}.junctScope.txt
EOF
			echo -e "sh Intermediate/${sample}/${sample}_exc.sh" >> ${OUTPUT}/${OUTPUT}_exc.sh
	fi
done < ${INPUT}


cat <<EOF > ${OUTPUT}/junctScopeSummarize.sh
#! /bin/env bash

# Go to project env
#cd \$(pwd)/${OUTPUT}/

# Set header
echo -e "sample\tgene\ttarget_junction_coord\ttarget_gene_coord\tnt_sequence\ttargetJxn_seq_count\ttargetJxn_coord_count\ttargetJxn_count_max\ttotal_sample_jxn_count\tgeneJxn_count_mean\tgeneJxn_count_median\tgeneJxn_count_sum\tgeneJxn_counts\tgeneJxn_transcripts" > ${OUTPUT}_output.txt

# Cat all sample outputs to single file
find . -name "*.junctScope.txt" -type f -exec cat {} + >> ${OUTPUT}_output.txt
EOF

