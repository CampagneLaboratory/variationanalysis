#!/usr/bin/env bash
FEATURE_MAPPER=$1
GPU=$2

if [ -e configure.sh ]; then
 echo "Loading configure.sh"
 source configure.sh
fi

if [ "$#" -ne 2 ]; then
   echo "Argument missing. You must provide a feature mapper classname to use in the iteration."
   exit 1;
fi
if [ -z "${LEARNING_RATE+set}" ]; then
    LEARNING_RATE="50"
    echo "LEARNING_RATE set to ${LEARNING_RATE}. Change the variable to switch the learning rate."
fi
if [ -z "${MINI_BATCH_SIZE+set}" ]; then
    MINI_BATCH_SIZE="128"
    echo "MINI_BATCH_SIZE set to ${MINI_BATCH_SIZE}. Change the variable to switch the mini-batch-size."
fi
if [ -z "${DATASET+set}" ]; then
    DATASET="NA12878-random-labeled-"
    echo "DATASET set to ${DATASET}. Change the variable to switch dataset."
fi
if [ -z "${GPU+set}" ]; then
    GPU="3"
    echo "GPU set to ${GPU}. Change the integer to switch to another GPU card."
fi
if [ -z "${VAL_SUFFIX+set}" ]; then
    VAL_SUFFIX="validation"
    echo "VAL_SUFFIX set to ${VAL_SUFFIX}. Change the variable to switch the validation suffix."
fi
if [ ! -e "${DATASET}${VAL_SUFFIX}.sbi" ]; then
    echo "The validation set was not found: ${DATASET}${VAL_SUFFIX}.sbi  "
       exit 1;
fi
if [ ! -e "${DATASET}train.sbi" ]; then
    echo "The training set was not found: ${DATASET}train.sbi  "
       exit 1;
fi
if [ ! -e "${DATASET}test.sbi" ]; then
        echo "The test set was not found: ${DATASET}test.sbi  "
           exit 1;
fi
echo "Iteration for FEATURE_MAPPER=${FEATURE_MAPPER}"

export FORCE_PLATFORM=native
#rm ${DATASET}train*.cf ${DATASET}${VAL_SUFFIX}*cf
train-genotype.sh 10g -t ${DATASET}train.sbi -v ${DATASET}${VAL_SUFFIX}.sbi \
 --mini-batch-size ${MINI_BATCH_SIZE}  -r ${LEARNING_RATE} --feature-mapper ${FEATURE_MAPPER} -x 10000 --build-cache-then-stop

export FORCE_PLATFORM=cuda
MODEL_DIR=`train-genotype.sh 10g -t ${DATASET}train.sbi -v ${DATASET}${VAL_SUFFIX}.sbi \
  --mini-batch-size ${MINI_BATCH_SIZE} -r ${LEARNING_RATE} --feature-mapper ${FEATURE_MAPPER} -x 10000 \
  --random-seed 90129 \
  --early-stopping-num-epochs 10 --gpu-device ${GPU} \
  | tee output-${RANDOM}.log |grep "model directory:"|cut -d " " -f 3`

`ls -1tr models|tail -1`
predict-genotypes.sh 10g -m ${MODEL_DIR} -l bestAUC+F1 -f -i ${DATASET}test.sbi --num-variants-expected 5205
