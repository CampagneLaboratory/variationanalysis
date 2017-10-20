#!/usr/bin/env bash
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

java -Xmx${memory_requirement} -cp ${GDLVA_JAR} -Dlogback.configurationFile=${SLF4J_CONFIG} \
     -Dframework.parallelWrapper.prefetchBuffer=32 -Dframework.parallelWrapper.averagingFrequency=7 -Dframework.parallelWrapper.numWorkers=8   \
      -Djava.io.tmpdir=${TMPDIR} \
    org.campagnelab.dl.genotype.tools.TrainModelGS ${other_parameters}
