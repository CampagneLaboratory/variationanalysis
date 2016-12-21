FEATURE_MAPPER=$1
if [ "$#" -ne 1 ]; then
   echo "Argument missing. You must provide a feature mapper classname to use in the iteration."
   exit 1;
fi

if [ -z "${DATASET+set}" ]; then
    DATASET="NA12878-random-labeled-"
    echo "DATASET set to ${DATASET}. Change the variable to switch dataset."
fi
if [ -z "${VAL_SUFFIX+set}" ]; then
    VAL_SUFFIX="validation"
    echo "VAL_SUFFIX set to ${VAL_SUFFIX}. Change the variable to switch the validation suffix."
fi
if [ -n -e ${DATASET}${VAL_SUFFIX}.sbi ]; then
    echo "The validation set was not found: ${DATASET}${VAL_SUFFIX}.sbi  "
       exit 1;
fi
if [ -n -e ${DATASET}train.sbi ]; then
    echo "The training set was not found: ${DATASET}train.sbi  "
       exit 1;
fi
if [ -n -e ${DATASET}test.sbi ]; then
        echo "The test set was not found: ${DATASET}test.sbi  "
           exit 1;
fi
echo "Iteration for FEATURE_MAPPER=${FEATURE_MAPPER}"
export FORCE_PLATFORM=native
train-genotype.sh 10g -t ${DATASET}training.sbi -v ${DATASET}validation.sbi \
 --mini-batch-size 2048 -r 5 --feature-mapper ${FEATURE_MAPPER} --build-cache-then-stop
export FORCE_PLATFORM=cuda
train-genotype.sh 10g -t ${DATASET}training.sbi -v ${DATASET}${VAL_SUFFIX}.sbi \
  --mini-batch-size 2048 -r 5 --feature-mapper ${FEATURE_MAPPER} --early-stopping-num-epochs 1
MODEL_TIMESTAMP=`ls -1tr models|tail -1`
predict-genotypes.sh 10g -m models/${MODEL_TIMESTAMP} -l bestscore -f -i ${DATASET}test.sbi
