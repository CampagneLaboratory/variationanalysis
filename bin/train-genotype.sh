#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

java -Djava.io.tmpdir=${TMPDIR} -Xmx${memory_requirement} -cp ${GDLVA_JAR} -Dlogback.configurationFile=${SLF4J_CONFIG} \
    -Dframework.parallelWrapper.prefetchBuffer=3 \
     -Dframework.parallelWrapper.averagingFrequency=3 \
     -Dframework.parallelWrapper.numWorkers=3   \
    org.campagnelab.dl.genotype.learning.TrainModelG ${other_parameters}
