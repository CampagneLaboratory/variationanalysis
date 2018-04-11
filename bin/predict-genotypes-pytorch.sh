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
OUTPUT_VEC=$3
INPUT_VEC_PATH=$4
MODEL=`basename ${MODEL_PATH}`

echo "usage: result will be written to sorted-${MODEL}-${MODEL_LABEL}.vcf.gz"

if [ -z "${MINI_BATCH_SIZE+set}" ]; then
    MINI_BATCH_SIZE="2048"
    echo "MINI_BATCH_SIZE set to ${MINI_BATCH_SIZE}. Change the variable to switch the mini-batch-size."
fi

MODEL_TIME=`basename ${MODEL_PATH}`
echo "Running predict-dataset to create output vector file from trained PyTorch model..."
predict-dataset.sh --model-path ${MODEL_PATH} --model-label ${MODEL_LABEL} \
   --problem "genotyping:${INPUT_VEC_PATH}" ${PREDICT_MAX_RECORDS} \
   --dataset ${DATASET_NAME} --mini-batch-size ${MINI_BATCH_SIZE} \
   --output "${OUTPUT_VEC}" --checkpoint-key ${CHECKPOINT_KEY} ${VEC_INPUTS}
