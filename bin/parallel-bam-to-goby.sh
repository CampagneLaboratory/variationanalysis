#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

assertGobyInstalled
assertParallelInstalled

ALIGNMENTS="$*"
if [ "$#" -eq 1 ]; then
   case ${ALIGNMENTS} in *.bam) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .bam`;; esac
   case ${ALIGNMENTS} in *.sam) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .sam`;; esac
   case ${ALIGNMENTS} in *.cram) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .cram`;; esac
else
    OUTPUT_BASENAME="out-concat"
fi
echo "Will write Goby alignment to ${OUTPUT_BASENAME}"

if [ -z "${SBI_GENOME+set}" ]; then
    SBI_GENOME="/data/genomes/Homo_sapiens.ucsc.hg19"
    echo "SBI_GENOME set to ${SBI_GENOME}. Change the variable to influence the genome used (must be indexed with goby build-sequence-cache)."
fi
if [ -z  "${SBI_NUM_THREADS+set}" ]; then
    SBI_NUM_THREADS="2"
    echo "SBI_NUM_THREADS set to ${SBI_NUM_THREADS}. Change the variable to influence the number of parallel jobs."
fi
echo "variables: ${SBI_GENOME} ${SBI_NUM_THREADS}"

goby ${memory_requirement} suggest-position-slices ${ALIGNMENTS} --number-of-slices 60 -o slices.tsv
grep -v targetIdStart slices.tsv >slices
echo " concatenate-alignments --genome  ${SBI_GENOME}  ${ALIGNMENTS} " >command.txt

cut -f3,6 slices  | awk 'BEGIN{count=1} {print "-s "$1" -e " $2" -o out-part-"(count++)}' >boundaries
parallel -j2 --plus  --progress goby ${memory_requirement}  `cat command.txt`  :::: boundaries

goby ${memory_requirement} concatenate-alignments out-part-*.sbi -o ${OUTPUT_BASENAME}
