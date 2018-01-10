#!/usr/bin/env bash

DISTRIBUTION=$1
shift
set -x

if [ -f ${DISTRIBUTION}/gatk ]; then
    bin/parallel-gatk4-realign-filtered.sh $1/gatk "$@"
elif [ -f ${DISTRIBUTION}/GenomeAnalysisTK.jar ]; then
    bin/parallel-gatk3-realign-filtered.sh $1/GenomeAnalysisTK.jar "$@"
else
    echo "Unable to execute GATK distribution at ${DISTRIBUTION}"
fi
