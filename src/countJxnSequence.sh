#! bin/bash/


# Example Use
# sh countJxnSequence.sh -d region.seq.sam -o region.seq.txt

unset -v INPUT OUTPUT

usage() { echo "Usage: $0 [-i <SAM file>] [-o <output text file>]" 1>&2; exit 1; }

OUTPUT="NULL"

while getopts ":i:o:" opt; do
  case ${opt} in
  i)
  	INPUT=${OPTARG}
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
  wc -l ${INPUT} | awk '{print $1}'
else
  wc -l ${INPUT} | awk '{print $1}' > ${OUTPUT}
fi


#echo -e "count_samtools\tSAM_file" > ${OUTPUT}
#wc -l ${DIR}/* | sed '$ d' | sed 's/^[[:space:]]*//' | awk '{print $1}' >> ${OUTPUT}
# list files and line count | remove last line | reorder columns with tab delim