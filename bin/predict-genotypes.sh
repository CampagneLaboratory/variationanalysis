#!/usr/bin/env bash
export FORCE_PLATFORM="native"
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

java -Xmx${memory_requirement} -cp ${GDLVA_JAR} -Dlogback.configurationFile=${SLF4J_CONFIG} \
    org.campagnelab.dl.genotype.tools.PredictG ${other_parameters}