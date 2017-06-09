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
export GOBY_DSV_OTHER_OPTIONS=" --covariates ${GOBY_DSV_COVARIATES} -x SomaticVariationOutputFormat:model-path=${SOMATIC_MODEL_PATH} "
export GOBY_GENOTYPE_FORMAT="SOMATIC_VARIATIONS"
parallel-calls.sh ${memory_requirement} ${GENOTYPE_MODEL_PATH} ${ALIGNMENTS}