#! /bin/bash


# Example Use
# sh extractJxnRegion.sh -i input.cram -r chr1:1000-2000 -o output.region.sam -b 200 -t 8

usage() { echo "Usage: $0 [-i <CRAM/BAM/SAM file>] [-r <chr1:1000-2000>] [-o <output SAM file>] [-b <int nucleotide region buffer>] [-t <int threads>] [-f <flag to filter out multimapped reads>]" 1>&2; exit 1; }

# Defaults
BUFFER=200
THREADS=1
FILTER=false

while getopts ":i:r:o:b:t:f" opt; do
  case ${opt} in
  i)
  	INPUT=${OPTARG}
  	;;
  r)
  	REGION=${OPTARG}
  	;;
  o)
  	OUTPUT=${OPTARG}
  	;;
  b)
  	BUFFER=${OPTARG}
  	;;
  t)
  	THREADS=${OPTARG}
  	;;
  f)
    FILTER=true
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

# Add buffer nt around region of interest
if [[ ${BUFFER} -gt 0 ]]
then
	chr=$(echo ${REGION} | cut -d ':' -f 1)
	pos=$(echo ${REGION} | cut -d ':' -f 2,3)
	str=$(($(echo ${pos} | cut -d '-' -f 1) - ${BUFFER}))
	end=$(($(echo ${pos} | cut -d '-' -f 2) + ${BUFFER}))
	REGION=${chr}:${str}-${end}
fi

if [ ${FILTER} ]; then
  samtools view -F 0x900 -@ ${THREADS} ${INPUT} ${REGION} > ${OUTPUT}
else
  samtools view -@ ${THREADS} ${INPUT} ${REGION} > ${OUTPUT}
fi


