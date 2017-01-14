#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/common.sh
FEATURE_MAPPER=$1
GPU=$2
PREVIOUS_MODEL_PATH=$3
PREVIOUS_MODEL_NAME=$4

if [ -e configure.sh ]; then
 echo "Loading configure.sh"
 source configure.sh
fi

cat << EOF | cat>configure-downsampling.sh
export TRAIN_SUFFIX="${TRAIN_SUFFIX}-ds"
export VAL_SUFFIX="${VAL_SUFFIX}-ds"
export TRAINING_OPTIONS="${TRAINING_OPTIONS} --previous-model-path ${PREVIOUS_MODEL_PATH} --previous-model-name ${PREVIOUS_MODEL_NAME} "
EOF

echo "Continuing training with downsampled training and validation sets.."
iterate-genotype.sh ${FEATURE_MAPPER} ${GPU}