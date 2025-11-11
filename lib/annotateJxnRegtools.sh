#! /bin/bash

# Example Use
# sh annotateJxnRegtools.sh -i input.region.bed -f ref.fa -g ref.gtf -r -o output.region.bed.anno


usage() { echo "Usage: $0 [-i <BED12 file>] [-f <FASTA reference>] [-g <GTF reference>] [-o <output text file>]" 1>&2; exit 1; }

while getopts ":i:f:g:o:" opt; do
  case ${opt} in
  i)
  	INPUT=${OPTARG}
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

if [ -z "${INPUT}" ] || [ -z "${FASTA}" ] || [ -z "${GTF}" ] || [ -z "${OUTPUT}" ]; then
    usage
fi


regtools junctions annotate -o ${OUTPUT} ${INPUT} ${FASTA} ${GTF}


