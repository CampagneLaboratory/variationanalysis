#!/usr/bin/env bash
memory_requirement=$1
shift
other_parameters=$*

export DISTRIBUTION_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
if [[ $OSTYPE == "cygwin" ]]; then
    DISTRIBUTION_DIR=`cygpath -m "${DISTRIBUTION_DIR}"`
fi
export VERSION=`cat ${DISTRIBUTION_DIR}/VERSION.txt`
export DLVA_HOME=${DISTRIBUTION_DIR}
export DLVA_JAR=${DLVA_HOME}/gpus/target/gpus-${VERSION}.jar:${DLVA_HOME}/somatic/target/somatic-${VERSION}-bin.jar
export SLF4J_CONFIG=${DLVA_HOME}/config/logback.xml