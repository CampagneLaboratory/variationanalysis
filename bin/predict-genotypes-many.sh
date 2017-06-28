#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh
. configure.sh

assertGobyInstalled
assertParallelInstalled
if [ $# -lt 3 ]; then
 echo "usage: predict-genotypes-many MODEL_PATH MODEL_LABEL [out-part-*.sbi]"
 echo "usage: e.g., predict-genotypes-many /home/campagne/data/genotypes/jun_19/CNG_NA12878_model/1498053463106 bestscore out-part-*.sbi"

 exit 1;
fi

MODEL_PATH=$1
MODEL_LABEL=$2
shift
shift
SBI_INPUTS="$*"
MODEL=`basename ${MODEL_PATH}`

echo "usage: result will be written to sorted-${MODEL}-${MODEL_LABEL}.vcf.gz"
if [ -z  "${SBI_NUM_THREADS+set}" ]; then
    SBI_NUM_THREADS="2"
    echo "SBI_NUM_THREADS set to ${SBI_NUM_THREADS}. Change the variable to influence the number of parallel jobs."
fi

if [ -z "${MINI_BATCH_SIZE+set}" ]; then
    MINI_BATCH_SIZE="2048"
    echo "MINI_BATCH_SIZE set to ${MINI_BATCH_SIZE}. Change the variable to switch the mini-batch-size."
fi

echo " predict-genotypes.sh 10g -m ${MODEL_PATH} -l ${MODEL_LABEL}  --no-cache --mini-batch-size ${MINI_BATCH_SIZE} \
-f --format VCF -i   " >command.txt

# keep a log of the commands that were used to generate this dataset:
cp command.txt command-`date +%h_%d_%H_%M`.txt

ls -1 ${SBI_INPUTS} >files
parallel --bar --eta -j${SBI_NUM_THREADS} --plus  --progress `cat command.txt`  :::: files

cat ${MODEL}-${MODEL_LABEL}-*.vcf | vcf-sort >sorted-${MODEL}-${MODEL_LABEL}.vcf
bgzip -f sorted-${MODEL}-${MODEL_LABEL}.vcf
tabix -f sorted-${MODEL}-${MODEL_LABEL}.vcf.gz