#! bin/bash/

# Example Use
# sh countJxnSequence.sh -d /.region.seq.sam/file/directory/ -o outputJxnCount.region.seq.txt


usage() { echo "Usage: $0 [-i <SAM file directory>] [-o <output text file>]" 1>&2; exit 1; }


while getopts ":d:o:" opt; do
  case ${opt} in
  d)
  	DIR=${OPTARG}
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

if [ -z "${DIR}" ] || [ -z "${OUTPUT}" ]; then
    usage
fi

wc -l ${DIR}/* > ${OUTPUT}
sed -i '$ d' ${OUTPUT}


