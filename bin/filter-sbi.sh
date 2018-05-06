#!/usr/bin/env bash
FORCE_PLATFORM=native
. `dirname "${BASH_SOURCE[0]}"`/setup.sh

java -ea -Xmx${memory_requirement} -cp ${GDLVA_JAR} -Dlogback.configurationFile=${SLF4J_CONFIG} \
    org.campagnelab.dl.genotype.tools.FilterSBI ${other_parameters}