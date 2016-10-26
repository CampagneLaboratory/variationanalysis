#!/usr/bin/env bash
memory_requirement=$1
shift
other_parameters=$*
set +x
WORKING_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if [[ $OSTYPE == "cygwin" ]]; then
    WORKING_DIR=`cygpath -m "${WORKING_DIR}"`
fi

export DLVA_HOME=${WORKING_DIR}
DLVA_JAR=${DLVA_HOME}/model-training/target/model-training-1.0.3-SNAPSHOT-bin.jar
SLF4J_CONFIG=${DLVA_HOME}/config/goby-logback.xml

java -Xmx${memory_requirement} -cp ${DLVA_JAR} -Dlogback.configurationFile=${SLF4J_CONFIG}  ${DLVA_JAR} \
    org.campagnelab.dl.varanalysis.tools.Randomize ${other_parameters}