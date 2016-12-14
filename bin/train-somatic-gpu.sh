#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

java -Djava.io.tmpdir=${TMPDIR}  -Xmx${memory_requirement} -cp ${DLVA_JAR} -Dlogback.configurationFile=${SLF4J_CONFIG}   \
    org.campagnelab.dl.somatic.learning.TrainSomaticModelOnGPU ${other_parameters}