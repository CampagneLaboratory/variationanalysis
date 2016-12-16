#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

ALIGNMENTS="$*"


if [ -z "${SBI_GENOME+set}" ]; then
    SBI_GENOME="/data/genomes/Homo_sapiens.ucsc.hg19"
    echo "SBI_GENOME set to ${SBI_GENOME}. Change the variable to influence the genome used (must be indexed with goby build-sequence-cache)."
fi
if [ -z  "${SBI_NUM_THREADS+set}" ]; then
    SBI_NUM_THREADS="2"
    echo "SBI_NUM_THREADS set to ${SBI_NUM_THREADS}. Change the variable to influence the number of parallel jobs."
fi
echo "variables: ${SBI_GENOME} ${SBI_NUM_THREADS}"

goby ${memory_requirement} suggest-position-slices ${ALIGNMENTS} --number-of-slices 200 -o slices.tsv
grep -v targetIdStart slices.tsv >slices
echo " discover-sequence-variants -n 1 -t 1 --genome  ${SBI_GENOME} --format  SEQUENCE_BASE_INFORMATION  ${ALIGNMENTS} --max-coverage-per-site 10000" >command.txt

cut -f3,6 slices  | awk 'BEGIN{count=1} {print "-s "$1" -e " $2" -o out-part-"(count++)}' >boundaries
parallel -j2 --plus  --progress goby ${memory_requirement}  `cat command.txt`  :::: boundaries



concat.sh ${memory_requirement} out-part-*.sbi -o out-concat.sbi
