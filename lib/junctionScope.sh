#! /bin/env bash

# Example Use
# sh junctionScope.sh -i sample_inputfiles.txt -r chr1:1000-2000 -n ATCGGACATT -f ref.fa -g ref.gtf -o junctionScope_output -b 200 -s 0 -t 8

usage() { echo "Usage: $0 [-i <list file>] [-r <chr1:1000-2000>] [-n <nucleotide sequence>] [-f <FASTA reference>] [-g <GTF reference>] [-o <output directory>] [-b <int nucleotide region buffer>] [-s <0|1|2>] [-t <int threads>]" 1>&2; exit 1; }

# built with
# ml SAMtools/1.9-foss-2018b
# ml RegTools/0.5.2-foss-2021b

# Defaults
BUFFER=200
THREADS=1
STRAND=0

while getopts ":i:r:n:f:g:o:b:s:t:" opt; do
  case ${opt} in
  i)
  	INPUT=${OPTARG}
  	;;
  r)
  	REGION=${OPTARG}
  	;;
  n)
  	NTSEQ=${OPTARG}
  	;;
  f)
  	FASTA=${OPTARG}
  	;;
  g)
  	GTF=${OPTARG}
  	;;
  o)
  	OUTPUT=${OPTARG}
  	;;
  b)
  	BUFFER=${OPTARG}
  	;;
  s)
  	STRAND=${OPTARG}
  	;;
  t)
  	THREADS=${OPTARG}
  	;;
  	\?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ -z "${INPUT}" ] || [ -z "${REGION}" ] || [ -z "${NTSEQ}" ] || [ -z "${FASTA}" ] || [ -z "${GTF}" ]; then
    usage
fi

mkdir -p ${OUTPUT}
mkdir -p ${OUTPUT}/extractJxnRegion
mkdir -p ${OUTPUT}/extractJxnSequence
mkdir -p ${OUTPUT}/extractJxnRegtools
mkdir -p ${OUTPUT}/annotateJxnRegtools

while IFS=$'\t' read -r sample file; do
	file="${file//[$'\r\n']/}"
	# setup
	# samtools extract
	echo "sh extractJxnRegion.sh -i ${file} -r ${REGION} -o ${OUTPUT}/extractJxnRegion/${sample}.region.sam -b ${BUFFER} -t ${THREADS}" >> ${OUTPUT}/extractJxnRegion_batch.sh
	echo "sh extractJxnSequence.sh -i ${OUTPUT}/extractJxnRegion/${sample}.region.sam -n ${NTSEQ} -o ${OUTPUT}/extractJxnSequence/${sample}.region.seq.sam" >> ${OUTPUT}/extractJxnSequence_batch.sh
	# regtools extract
	echo "sh extractJxnRegtools.sh -i ${file} -r ${REGION} -s ${STRAND} -o ${OUTPUT}/extractJxnRegtools/${sample}.region.bed" >> ${OUTPUT}/extractJxnRegtools_batch.sh
	echo "sh annotateJxnRegtools.sh -i ${OUTPUT}/extractJxnRegtools/${sample}.region.bed -f ${FASTA} -g ${GTF} -o ${OUTPUT}/annotateJxnRegtools/${sample}.region.bed.anno" >> ${OUTPUT}/annotateJxnRegtools_batch.sh
done < ${INPUT}

echo "sh ${OUTPUT}/extractJxnRegion_batch.sh" >> junctionScope_exc.sh
echo "sh ${OUTPUT}/extractJxnSequence_batch.sh" >> junctionScope_exc.sh
echo "sh ${OUTPUT}/extractJxnRegtools_batch.sh" >> junctionScope_exc.sh
echo "sh ${OUTPUT}/annotateJxnRegtools_batch.sh" >> junctionScope_exc.sh
echo "sh countJxnSequence.sh -d ${OUTPUT}/extractJxnSequence/ -o ${OUTPUT}/output.JxnCount.txt" >> junctionScope_exc.sh
echo "sh matchJxnRegtools.sh -d ${OUTPUT}/annotateJxnRegtools/ -r ${REGION} -o ${OUTPUT}/outputJxnCount.regtools.txt" >> junctionScope_exc.sh

