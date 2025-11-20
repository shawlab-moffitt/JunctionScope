#! /bin/bash

# Example Use
# sh extractJxnSequence.sh -i input.region.sam -n ATCGGACATT -o output.seq.sam

usage() { echo "Usage: $0 [-i <SAM file>] [-n <nucleotide sequence>] [-o <output SAM file>]" 1>&2; exit 1; }

while getopts ":i:n:o:" opt; do
  case ${opt} in
  i)
  	INPUT=${OPTARG}
  	;;
  n)
  	NTSEQ=${OPTARG}
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

if [ -z "${INPUT}" ] || [ -z "${NTSEQ}" ] || [ -z "${OUTPUT}" ]; then
    usage
fi

grep --no-filename ${NTSEQ} ${INPUT} > ${OUTPUT}





