#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

assertGobyInstalled

ALIGNMENTS="$*"
if [ "$#" -eq 1 ]; then
   case ${ALIGNMENTS} in *.bam) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .bam`;; esac
   case ${ALIGNMENTS} in *.sam) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .sam`;; esac
   case ${ALIGNMENTS} in *.cram) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .cram`;; esac
else
    OUTPUT_BASENAME="out-concat"
fi
echo "Will write Goby alignment to ${OUTPUT_BASENAME}"

if [ -z "${INCLUDE_INDELS+set}" ]; then
    INCLUDE_INDELS="false"
    echo "INCLUDE_INDELS set to ${INCLUDE_INDELS}. Change the variable to influence whether indels are included in the SBI file."
fi
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
echo " discover-sequence-variants -n 1 -t 1 --genome  ${SBI_GENOME} --format  SEQUENCE_BASE_INFORMATION  ${ALIGNMENTS} \
    --call-indels  ${INCLUDE_INDELS} --processor realign_near_indels \
    --max-coverage-per-site 10000" >command.txt

cut -f3,6 slices  | awk 'BEGIN{count=1} {print "-s "$1" -e " $2" -o out-part-"(count++)}' >boundaries
parallel -j2 --plus  --progress goby ${memory_requirement}  `cat command.txt`  :::: boundaries

concat.sh ${memory_requirement} -i out-part-*.sbi -o ${OUTPUT_BASENAME}.sbi