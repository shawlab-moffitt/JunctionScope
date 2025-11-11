#! /bin/bash

# Example Use
# sh extractJxnRegtools.sh -i input.cram -r chr1:1000-2000 -s XS -o output.region.bed

# Default unstranded
STRAND=0
# (0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR)

usage() { echo "Usage: $0 [-i <CRAM/BAM/SAM file>] [-r <chr1:1000-2000>] [-s <0|1|2>] [-o <BED12 file>]" 1>&2; exit 1; }

while getopts ":i:r:s:o:" opt; do
  case ${opt} in
  i)
  	INPUT=${OPTARG}
  	;;
  r)
  	REGION=${OPTARG}
  	;;
  s)
  	STRAND=${OPTARG}
  	((STRAND == 0 || STRAND == 1 || STRAND == 2)) || usage
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

if [ -z "${INPUT}" ] || [ -z "${REGION}" ] || [ -z "${OUTPUT}" ]; then
    usage
fi

regtools junctions extract -r ${REGION} -s ${STRAND} -o ${OUTPUT} ${INPUT}


