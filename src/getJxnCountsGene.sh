#! /bin/bash


# Example Use
# sh getJxnCountsGene.sh -i input.gene.bed.anno -o output.txt

usage() { echo "Usage: $0 [-i <bed anno file>] [-G <GENE SYMBOL>] [-o <TXT file>]" 1>&2; exit 1; }



while getopts ":i:G:o:" opt; do
  case ${opt} in
  i)
  	INPUT=${OPTARG}
  	;;
  G)
    GENE=${OPTARG}
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

if [ -z "${INPUT}" ] || [ -z "${OUTPUT}" ] || [ -z "${GENE}" ]; then
    usage
fi


# https://stackoverflow.com/a/41762669
median() {
  arr=($(printf '%d\n' "${@}" | sort -n))
  nel=${#arr[@]}
  if (( $nel % 2 == 1 )); then     # Odd number of elements
    val="${arr[ $(($nel/2)) ]}"
  else                             # Even number of elements
    (( j=nel/2 ))
    (( k=j-1 ))
    (( val=(${arr[j]} + ${arr[k]})/2 ))
  fi
  echo $val
}

# Get sum and array of counts
declare -a counts
declare -a transcriptIDs
sum=0
while IFS=$'\t' read -a line; do
  kd=${line[11]}
  ka=${line[12]}
  gn=${line[14]}
  if [[ ${kd} == "1" ]] && [[ ${ka} == "1" ]] && [[ ${gn} == ${GENE} ]]; then
    counts+=(${line[4]})
    transcriptIDs+=(${line[15]})
    sum=$((sum+${line[4]}))
  fi
done < ${INPUT}

# Get mean count
num_elem=(${#counts[@]})
if [ ${num_elem} -gt 0 ]; then
  mean=$((${sum}/${num_elem}))
else
  mean=${counts[0]}
fi

# Get median count
if [ ${num_elem} -gt 0 ]; then
  med=$(median ${counts[@]})
else
  med=${counts[0]}
fi

# Collapse count vector
counts_coll=$(printf "%s," ${counts[@]})
counts_coll=${counts_coll%,}

# Collapse vector of transcript IDs
transcriptIDs_coll=$(printf "%s," ${transcriptIDs[@]})
transcriptIDs_coll=${transcriptIDs_coll%,}

echo -e "${mean}\t${med}\t${sum}\t${counts_coll}\t${transcriptIDs_coll}" > ${OUTPUT}


