#!/usr/bin/env bash

USAGE_STR=$(cat <<-END
    Usage: monitor-model-training.sh -m model-directory -c checkpoint-key -p model-prefix -t test-sbi -d dataset [-n]
    Monitor the progress of a model being trained by PyTorch
    -m should point to the directory of the trained PyTorch model being used with a config.properties within.
    -c should be the checkpoint key used and -p should be the model prefix (usually best or latest)
    -t should be a path to the test SBI file, and -d should be the name of the dataset used for the sbi.
    -n disables writing progress out to a file
END
)

. `dirname "${BASH_SOURCE[0]}"`/common.sh

if [ $# -eq 0 ]; then
    echo "${USAGE_STR}"
    exit 0
fi

function assertEvaluateInstalled {
    evaluate-genotypes-vec.sh -h >/dev/null 2>&1 || { echo >&2 "This script requires evaluate-genotypes-vec.sh from Variation to be in your path. Aborting. Check to make sure it is in your path before running, then try again."; exit 1; }
}

assertEvaluateInstalled

function assertLogInstalled {
    log-evaluate.sh -h >/dev/null 2>&1 || { echo >&2 "This script requires log-evaluate.sh from GenotypeTensors to be in your path. Aborting. Check to make sure it is in your path before running, then try again."; exit 1; }
}

assertLogInstalled

MODEL_DIR=""
CHECKPOINT_KEY=""
MODEL_PREFIX=""
DATASET_SBI=""
DATASET_NAME=""
NONVERBOSE=0

while getopts ":hm:c:p:t:d:n" opt; do
    case "${opt}" in
        h)
            echo "${USAGE_STR}"
            exit 0
            ;;
        m)
            MODEL_DIR=$OPTARG
            ;;
        c)
            CHECKPOINT_KEY=$OPTARG
            ;;
        p)
            MODEL_PREFIX=$OPTARG
            ;;
        t)
            DATASET_SBI=$OPTARG
            ;;
        d)
            DATASET_NAME=$OPTARG
            ;;
        n)
            NONVERBOSE=1
            ;;
        \?)
            echo "Invalid option: -${OPTARG}" 1>&2
            exit 1;
            ;;
        :)
            echo "Invalid Option: -$OPTARG requires an argument" 1>&2
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))

ERROR_STR=$(cat <<-END
   You must specify all options: -m, -c, -p, -t, and -d. If you are unsure about usage, run with no arguments or with -h.
END
)

if [ -z "${MODEL_DIR}" ] || [ -z "${CHECKPOINT_KEY}" ] || [ -z "${MODEL_PREFIX}" ] || [ -z "${DATASET_SBI}" ] || [ -z "${DATASET_NAME}" ]; then
    echo "${ERROR_STR}"
    exit 1
fi

FULL_MODEL_PATH=$(python - <<EOF
import os
print(os.path.join("${MODEL_DIR}", "models", "pytorch_{}_{}.t7".format("${CHECKPOINT_KEY}", "${MODEL_PREFIX}")))
EOF
)

if [ ! -e "${FULL_MODEL_PATH}" ]; then
    echo "The model was not found: ${FULL_MODEL_PATH}"
    exit 1;
fi

echo "Logging changes on ${FULL_MODEL_PATH}..."

PREV_MODEL_MD5=""

LOG_PROGRESS_PATH="${CHECKPOINT_KEY}_progress.log"
if [ ${NONVERBOSE} -ne 0 ]; then
    LOG_PROGRESS_PATH="/dev/null"
fi

# Compare current hash against previous; if different, run evaluate and log results from evaluate to log file
while sleep 1
do
    CURR_MODEL_MD5=$(md5sum "${FULL_MODEL_PATH}")
    if [[ "${CURR_MODEL_MD5}" != "${PREV_MODEL_MD5}" ]] ; then
        echo "=========Change detected on ${FULL_MODEL_PATH}========="
        RANDOM_OUTPUT_SUFFIX="${RANDOM}"
        NOW=$(date)
        echo "========= NEW EVALUATION RUN ON ${NOW} FOR MODEL ${FULL_MODEL_PATH} ON DATASET ${DATASET_SBI} =========" >> ${LOG_PROGRESS_PATH}
        echo "Running evaluate on ${FULL_MODEL_PATH} with dataset ${DATASET_SBI}..."
        evaluate-genotypes-vec.sh -m "${MODEL_DIR}" -c "${CHECKPOINT_KEY}" -p "${MODEL_PREFIX}" -t "${DATASET_SBI}" -d "${DATASET_NAME}" -r "${RANDOM_OUTPUT_SUFFIX}" >>"${LOG_PROGRESS_PATH}" 2>&1
        dieIfError "Unable to evaluate genotypes with the given parameters"
        echo "Logging results of evaluate on ${FULL_MODEL_PATH} with dataset ${DATASET_SBI}..."
        log-evaluate.sh --dataset ${DATASET_NAME} --model-path "${MODEL_DIR}" --checkpoint-key "${CHECKPOINT_KEY}" --model-label "${MODEL_PREFIX}" --vcf-path "output-${RANDOM_OUTPUT_SUFFIX}" --output-path "${CHECKPOINT_KEY}_log.tsv" >>"${LOG_PROGRESS_PATH}" 2>&1
        dieIfError "Unable to log results of evaluation with the given parameters"
        PREV_MODEL_MD5="${CURR_MODEL_MD5}"
        rm -r "output-${RANDOM_OUTPUT_SUFFIX}"
        echo "Evaluation and logging for next epoch of ${FULL_MODEL_PATH} complete"
    fi
done