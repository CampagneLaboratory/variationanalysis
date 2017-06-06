#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

assertGobyInstalled
assertParallelInstalled
assertBctftoolsInstalled
. ./configure.sh

# Call Genotypes across an alignment, using a DL variationanalysis model and Goby3:
# Usage: parallel-calls.sh 10g model alignment.entries
MODEL_PATH=$1
shift
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

if [ -z "${GOBY_NUM_SLICES+set}" ]; then
    GOBY_NUM_SLICES="50"
    echo "GOBY_NUM_SLICES set to ${GOBY_NUM_SLICES}. Change the variable to influence whether how many slices to use when processing an alignment."
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

if [ -z "${REALIGN_AROUND_INDELS+set}" ]; then
    echo "Set REALIGN_AROUND_INDELS to true to enable realignment around indels."
    REALIGNMENT_OPTION=" "
else
   if [ "${REALIGN_AROUND_INDELS}" == "true" ]; then
       echo "REALIGN_AROUND_INDELS set to ${REALIGN_AROUND_INDELS}. ENABLED realignment around indels."
       REALIGNMENT_OPTION="--processor realign_near_indels"
   else
       echo "REALIGN_AROUND_INDELS set to ${REALIGN_AROUND_INDELS}. DISABLED realignment around indels."
       REALIGNMENT_OPTION=" "
   fi
fi


echo "variables: ${SBI_GENOME} ${SBI_NUM_THREADS}"

set -x
goby 8g suggest-position-slices ${ALIGNMENTS} --modulo 1000 --restrict-per-chromosome --number-of-slices ${GOBY_NUM_SLICES} -o slices.tsv
grep -v targetIdStart slices.tsv >slices

if [  -z "${LIMIT_TO_CHROMOSOME+set}" ]; then
  echo "Set LIMIT_TO_CHROMOSOME to call only one chromosome";
 else
  # eliminate any slice that does not match the chromosome:
  grep "${LIMIT_TO_CHROMOSOME}," slices >filtered_slices
  mv filtered_slices slices
fi

echo " discover-sequence-variants -n 0 -t 1 --genome  ${SBI_GENOME} --format  GENOTYPES  ${ALIGNMENTS} \
    --call-indels  ${INCLUDE_INDELS} ${REALIGNMENT_OPTION} \
    --max-coverage-per-site 1000 -x HTSJDKReaderImpl:force-sorted=true \
    -x GenotypesOutputFormat:model-path=${MODEL_PATH} " >command.txt

# keep a log of the commands that were used to generate this dataset:
cp command.txt command-`date +%h_%d_%H_%M`.txt

cut -f3,6 slices  | awk 'BEGIN{count=1} {print "-s "$1" -e " $2" -o calls-part-"(count++)".vcf" }' >boundaries
parallel --bar --eta -j${SBI_NUM_THREADS} --plus  --progress goby ${memory_requirement} `cat command.txt`  :::: boundaries

# keep a log of the commands that were used to generate this dataset:
cp command.txt command-`date +%h_%d_%H_%M`.txt

echo "Concatenating VCF files."
goby ${memory_requirement} fdr calls-part-*.vcf -o ${OUTPUT_BASENAME}.vcf
bgzip ${OUTPUT_BASENAME}.vcf
bcftools index ${OUTPUT_BASENAME}.vcf.gz
