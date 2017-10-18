#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

java -Djava.io.tmpdir=${TMPDIR} -Xmx${memory_requirement} -cp ${GDLVA_JAR} -Dlogback.configurationFile=${SLF4J_CONFIG}   \
    org.campagnelab.dl.genotype.learning.TrainModelGS ${other_parameters}
