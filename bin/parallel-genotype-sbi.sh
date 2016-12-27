#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

assertGobyInstalled
assertParallelInstalled

ALIGNMENTS="$*"
if [ -z "${OUTPUT_BASENAME+set}" ]; then

    if [ "$#" -eq 1 ]; then
       case ${ALIGNMENTS} in *.bam) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .bam`;; esac
       case ${ALIGNMENTS} in *.sam) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .sam`;; esac
       case ${ALIGNMENTS} in *.cram) OUTPUT_BASENAME=`basename ${ALIGNMENTS} .cram`;; esac
    else
        OUTPUT_BASENAME="out-concat"
    fi

    echo "OUTPUT_BASENAME set to ${OUTPUT_BASENAME}. Change the variable to influence the output basename."
fi

echo "Will write Goby alignment to ${OUTPUT_BASENAME}"
if [ -z "${SBI_GENOTYPE_VARMAP+set}" ]; then
    echo "Set "SBI_GENOTYPE_VARMAP to the filename of a varmap to annotate true genotypes while you generate the .sbi file."
    VARMAP_OPTION=" "
else
   if [ -z "${REF_SAMPLING_RATE+set}" ]; then
        REF_SAMPLING_RATE="0.1"
        echo "Define REF_SAMPLING_RATE to change the default sampling rate. Set to ${REF_SAMPLING_RATE}"
    fi
    VARMAP_OPTION=" -x SequenceBaseInformationOutputFormat:random-seed=3764 -x SequenceBaseInformationOutputFormat:sampling-rate=${REF_SAMPLING_RATE} -x SequenceBaseInformationOutputFormat:true-genotype-map=${SBI_GENOTYPE_VARMAP} "
fi


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
echo " discover-sequence-variants -n 0 -t 1 --genome  ${SBI_GENOME} --format  SEQUENCE_BASE_INFORMATION  ${ALIGNMENTS} \
    --call-indels  ${INCLUDE_INDELS} --processor realign_near_indels \
    --max-coverage-per-site 10000 -x HTSJDKReaderImpl:force-sorted=true ${VARMAP_OPTION} " >command.txt

cut -f3,6 slices  | awk 'BEGIN{count=1} {print "-s "$1" -e " $2" -o out-part-"(count++)}' >boundaries
parallel -j${SBI_NUM_THREADS} --plus  --progress goby ${memory_requirement}  `cat command.txt`  :::: boundaries

concat.sh ${memory_requirement} -i out-part-*.sbi -o ${OUTPUT_BASENAME}