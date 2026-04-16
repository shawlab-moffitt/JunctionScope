#! /bin/bash


# Example Use
# sh extractJxnRegtools.sh -i input.cram -r chr1:1000-2000 -b 200 -o output.region.bed

# Default unstranded
#STRAND=0
# (0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR)

#usage() { echo "Usage: $0 [-i <CRAM/BAM file>] [-r <chr1:1000-2000>] [-s <0|1|2>] [-b <int nucleotide region buffer>] [-o <BED12 file>]" 1>&2; exit 1; }
usage() { echo "Usage: $0 [-i <CRAM/BAM file>] [-r <chr1:1000-2000>] [-b <int nucleotide region buffer>] [-o <BED12 file>]" 1>&2; exit 1; }


BUFFER=200

while getopts ":i:r:b:o:" opt; do
#while getopts ":i:r:s:b:o:" opt; do
  case ${opt} in
  i)
  	INPUT=${OPTARG}
  	;;
  r)
  	REGION=${OPTARG}
  	;;
  #s)
  #	STRAND=${OPTARG}
  #	((STRAND == 0 || STRAND == 1 || STRAND == 2)) || usage
  #	;;
  b)
    BUFFER=${OPTARG}
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

# Add buffer nt around region of interest
if [[ ${BUFFER} -gt 0 ]]
then
  chr=$(echo ${REGION} | cut -d ':' -f 1)
  pos=$(echo ${REGION} | cut -d ':' -f 2,3)
  str=$(($(echo ${pos} | cut -d '-' -f 1) - ${BUFFER}))
  end=$(($(echo ${pos} | cut -d '-' -f 2) + ${BUFFER}))
  REGION=${chr}:${str}-${end}
fi

# Check regtools version to catch strand argument input change
#reg_ver=$(regtools --version 2>&1 | awk '{print $2}' | sed -n 3p)
#if [[ "${reg_ver:0:1}" == "1" ]]; then
#  if [[ ${STRAND} == 0 ]]; then
#    ${STRAND}="XS"
#  elif [[ ${STRAND} == 1 ]]; then
#    ${STRAND}="RF"
#  elif [[ ${STRAND} == 2 ]]; then
#    ${STRAND}="FR"
#  fi
#fi

STRAND0=0
STRAND1=1
STRAND2=2
reg_ver=$(regtools --version 2>&1 | awk '{print $2}' | sed -n 3p)
if [[ "${reg_ver:0:1}" == "1" ]]; then
  ${STRAND0}="XS"
  ${STRAND1}="RF"
  ${STRAND2}="FR"
fi

#regtools junctions extract -r ${REGION} -s ${STRAND} -o ${OUTPUT} ${INPUT}

regtools junctions extract -r ${REGION} -s ${STRAND0} -o ${OUTPUT} ${INPUT}
regtools junctions extract -r ${REGION} -s ${STRAND1} ${INPUT} | cat >> ${OUTPUT}
regtools junctions extract -r ${REGION} -s ${STRAND2} ${INPUT} | cat >> ${OUTPUT}

