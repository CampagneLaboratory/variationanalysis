#!/usr/bin/env bash -x
memory_requirement=$1
shift
other_parameters=$*

WORKING_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
if [[ $OSTYPE == "cygwin" ]]; then
    WORKING_DIR=`cygpath -m "${WORKING_DIR}"`
fi

export DLVA_HOME=${WORKING_DIR}
DLVA_JAR=${DLVA_HOME}/model-training/target/model-training-1.0.3-SNAPSHOT-bin.jar
SLF4J_CONFIG=${DLVA_HOME}/model-training/config/logback.xml