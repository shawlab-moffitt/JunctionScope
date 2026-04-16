#! /bin/bash


# Example Use
# sh getJxnCountsSample.sh -i input.cram -t 8 -o output.txt

usage() { echo "Usage: $0 [-i <CRAM/BAM/SAM file>] [-t <int threads>] [-o <output TXT file>]" 1>&2; exit 1; }

# Defaults
THREADS=1
OUTPUT="NULL"

while getopts ":i:t:o:" opt; do
  case ${opt} in
  i)
    INPUT=${OPTARG}
    ;;
  t)
    THREADS=${OPTARG}
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


if [ -z "${INPUT}" ]; then
    usage
fi

if [[ ${OUTPUT} == "NULL" ]]; then
  samtools view -F 0x900 -@ ${THREADS} ${INPUT} | cut -f 6 | grep -o "N" | wc -l | awk '{print $1}'
else
  samtools view -F 0x900 -@ ${THREADS} ${INPUT} | cut -f 6 | grep -o "N" | wc -l | awk '{print $1}'> ${OUTPUT}
fi





