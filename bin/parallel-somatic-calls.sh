#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

assertGobyInstalled
assertParallelInstalled
assertBctftoolsInstalled
. ./configure.sh

# Call somatic variations across one or more alignments, using DL variationanalysis models and Goby3:
# Usage: parallel-somatic-calls.sh 10g covariates.tsv somatic_model genotype_model alignment.entries [alignment.entries]+
GOBY_DSV_COVARIATES="$1"
SOMATIC_MODEL_PATH="$2"
GENOTYPE_MODEL_PATH="$3"
shift
shift
shift
ALIGNMENTS="$*"


if [ -z "${FOCUS_VARMAP+set}" ]; then
    FOCUS_VARMAP_OPTION=" "
    echo "FOCUS_VARMAP not set.  Set the variable to path of varmap to restrict output to only sites in the varmap."
else
    FOCUS_VARMAP_OPTION=" -x SomaticVariationOutputFormat:focus-on-varmap=${FOCUS_VARMAP} "
fi

export GOBY_DSV_OTHER_OPTIONS=" --covariates ${GOBY_DSV_COVARIATES} ${FOCUS_VARMAP_OPTION} -x SomaticVariationOutputFormat:model-path=${SOMATIC_MODEL_PATH} "
export GOBY_GENOTYPE_FORMAT="SOMATIC_VARIATIONS"
parallel-calls.sh ${memory_requirement} ${GENOTYPE_MODEL_PATH} ${ALIGNMENTS}