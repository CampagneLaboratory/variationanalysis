#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

assertGobyInstalled
assertParallelInstalled
assertBctftoolsInstalled
. ./configure.sh

# Call somatic variations across one or more alignments, using DL variationanalysis models and Goby3:
# Usage: parallel-somatic-calls.sh 10g somatic_model genotype_model alignment.entries [alignment.entries]+
SOMATIC_MODEL_PATH=$1
GENOTYPE_MODEL_PATH=$2
shift
shift
ALIGNMENTS="$*"

export GOBY_DSV_OPTIONS=" -x SomaticVariationOutputFormat:model-path=${SOMATIC_MODEL_PATH} "
export GOBY_GENOTYPE_FORMAT="SOMATIC_VARIATIONS"
parallel-calls.sh ${memory_requirement} ${GENOTYPE_MODEL_PATH}
