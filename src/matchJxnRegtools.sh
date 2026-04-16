#! /bin/env bash


# Example Use
# sh matchJxnRegtools.sh -d region.bed.anno -r chr1:1000-2000 -o JxnCount.regtools.txt


usage() { echo "Usage: $0 [-i <annotated BED file>] [-G <GENE SYMBOL>] [-r <chr1:1000-2000>] [-f <flag to report region offset>] [-o <output text file>]" 1>&2; exit 1; }

OUTPUT="NULL"
OFFSET=0

while getopts ":i:G:r:o:f" opt; do
  case ${opt} in
  i)
    INPUT=${OPTARG}
    ;;
  G)
    GENE=${OPTARG}
    ;;
  r)
    REGION=${OPTARG}
    ;;
  f)
    OFFSET=1
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

if [ -z "${INPUT}" ] || [ -z "${REGION}" ] || [ -z "${GENE}" ]; then
    usage
fi

chrom=$(echo "${REGION}" | cut -d: -f1)
start=$(echo "${REGION}" | cut -d: -f2 | cut -d- -f1)
end=$(echo "${REGION}" | cut -d: -f2 | cut -d- -f2)

start_min=$((start - 1))
start_max=$((start + 1))
end_min=$((end - 1))
end_max=$((end + 1))


if [ ${OFFSET} == 1 ]; then
  countOut=$(awk -v chrom="$chrom" -v gene="$GENE" \
    -v sm="$start_min" -v s="$start" -v sp="$start_max" \
    -v em="$end_min" -v e="$end" -v ep="$end_max" \
    '
    $1==chrom && $2>=sm && $2<=sp && $3>=em && $3<=ep && $15==gene {
      start_offset = ($2==s) ? "0" : (($2==sm) ? "-1" : "+1")
      end_offset   = ($3==e) ? "0" : (($3==em) ? "-1" : "+1")
      print $5"\t"start_offset"\t"end_offset
      found = 1
  }
  END {
      if (!found) 
          print "\tNA\tNA\tNA"
  }
  ' ${INPUT})
else
  countOut=$(awk -v chrom="$chrom" -v gene="$GENE" \
    -v sm="$start_min" -v s="$start" -v sp="$start_max" \
    -v em="$end_min" -v e="$end" -v ep="$end_max" \
  '
  $1==chrom && $2>=sm && $2<=sp && $3>=em && $3<=ep && $15==gene {
      start_offset = ($2==s) ? "0" : (($2==sm) ? "-1" : "+1")
      end_offset   = ($3==e) ? "0" : (($3==em) ? "-1" : "+1")
      print $5
      found = 1
  }
  END {
      if (!found) 
          print "NA"
  }
  ' ${INPUT})
fi



if [[ ${OUTPUT} != "NULL" ]]; then
  echo -e ${countOut} > ${OUTPUT}
else
  echo -e ${countOut}
fi
