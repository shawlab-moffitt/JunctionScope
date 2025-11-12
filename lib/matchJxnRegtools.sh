#! /bin/env bash

# Example Use
# sh matchJxnRegtools.sh -d /region.bed.anno/file/directory/ -r chr1:1000-2000 -o outputJxnCount.regtools.region.anno.txt


usage() { echo "Usage: $0 [-d <annotated BED file directory>] [-r <chr1:1000-2000>] [-o <output text file>]" 1>&2; exit 1; }

while getopts ":d:r:o:" opt; do
  case ${opt} in
  d)
    DIR=${OPTARG}
    ;;
  r)
    REGION=${OPTARG}
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

if [ -z "${DIR}" ] || [ -z "${REGION}" ] || [ -z "${OUTPUT}" ]; then
    usage
fi

chrom=$(echo "${REGION}" | cut -d: -f1)
start=$(echo "${REGION}" | cut -d: -f2 | cut -d- -f1)
end=$(echo "${REGION}" | cut -d: -f2 | cut -d- -f2)

start_min=$((start - 1))
start_max=$((start + 1))
end_min=$((end - 1))
end_max=$((end + 1))

echo -e "regtools_file\tregtools_score\tstart_offset\tend_offset" > "${OUTPUT}"

for f in "${DIR}"/*; do
    match_count=0

    awk -v chrom="$chrom" \
        -v sm="$start_min" -v s="$start" -v sp="$start_max" \
        -v em="$end_min" -v e="$end" -v ep="$end_max" \
        -v fname="$(basename "$f")" \
    '
    $1==chrom && $2>=sm && $2<=sp && $3>=em && $3<=ep {
        start_offset = ($2==s) ? "0" : (($2==sm) ? "-1" : "+1")
        end_offset   = ($3==e) ? "0" : (($3==em) ? "-1" : "+1")
        print fname"\t"$5"\t"start_offset"\t"end_offset
        found = 1
    }
    END {
        if (!found) 
            print fname"\t0\tNA\tNA"
    }
    ' "$f" >> "${OUTPUT}"

done
