#!/usr/bin/env bash
memory_requirement=$1
shift
other_parameters=$*

WORKING_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
if [[ $OSTYPE == "cygwin" ]]; then
    WORKING_DIR=`cygpath -m "${WORKING_DIR}"`
fi
VERSION=`cat ../VERSION.txt`
export DLVA_HOME=${WORKING_DIR}
DLVA_JAR=${DLVA_HOME}/gpus/target/gpus-${VERSION}.jar:${DLVA_HOME}/somatic/target/somatic-${VERSION}-bin.jar
SLF4J_CONFIG=${DLVA_HOME}/config/logback.xml