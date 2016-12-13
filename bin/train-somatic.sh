#!/usr/bin/env bash
export FORCE_PLATFORM=cuda
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

java -Djava.io.tmpdir=${TMPDIR} -Xmx${memory_requirement} -cp ${DLVA_JAR} -Dlogback.configurationFile=${SLF4J_CONFIG}   \
    org.campagnelab.dl.somatic.learning.TrainModelS ${other_parameters}
