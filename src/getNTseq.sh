#! /bin/env bash



# Example Use
# sh getNTseq.sh -r chr1:1000-2000 -w 5 -f ref.fa

# NTSEQ=$(sh ${code_dir}getNTseq.sh -r ${JXNCOORD} -w 5 -f ${FASTA})

usage() { echo "Usage: $0 [-r <chr1:1000-2000>] [-f <FASTA reference>] [-w <int flanking window>] [-i <flag for reverse compliment output>] [-t <int threads>]" 1>&2; exit 1; }

# Defaults
REVCOMP=0
WINDOW=5

while getopts ":r:f:w:i" opt; do
  case ${opt} in
  r)
  	REGION=${OPTARG}
  	;;
  f)
  	FASTA=${OPTARG}
  	;;
  w)
    WINDOW=${OPTARG}
    ;;
  i)
  	REVCOMP=1
  	;;
  	\?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

if [ -z "${REGION}" ] || [ -z "${FASTA}" ]; then
    usage
fi

if [[ ${WINDOW} -gt 0 ]]
then
  chr=$(echo ${REGION} | cut -d ':' -f 1)
  pos=$(echo ${REGION} | cut -d ':' -f 2,3)
  str=$(echo ${pos} | cut -d '-' -f 1)
  end=$(echo ${pos} | cut -d '-' -f 2)
  left_start=$((str - WINDOW))
  left_end=$((str - 1))
  right_start=$((end + 1))
  right_end=$((end + WINDOW))
  REGIONs=${chr}:${left_start}-${left_end}
  REGIONe=${chr}:${right_start}-${right_end}
fi


tmpfile=$(mktemp)
trap 'rm -f "${tmpfile}"' EXIT
echo ${REGIONs} > ${tmpfile}
echo ${REGIONe} >> ${tmpfile}
if [ ${REVCOMP} == 1 ]; then
  samtools faidx ${FASTA} -r ${tmpfile} -i | sed -n '2p;4p' | paste -sd ''
else
  samtools faidx ${FASTA} -r ${tmpfile} | sed -n '2p;4p' | paste -sd ''
fi

