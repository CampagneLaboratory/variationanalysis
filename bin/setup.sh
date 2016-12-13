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
if [ -n "${FORCE_PLATFORM+set}" ]; then
    LATEST_PLATFORM=`cat ${DISTRIBUTION_DIR}/PROFILE.txt`
    if [ "$FORCE_PLATFORM" == "$LATEST_PLATFORM" ]; then
           export DLVA_JAR=${DLVA_HOME}/gpus/target/gpus-${VERSION}.jar:${DLVA_HOME}/somatic/target/somatic-${VERSION}-bin-${FORCE_PLATFORM}.jar
    fi
fi

if [ -n "${DLVA_JAR+set}" ]; then
    export DLVA_JAR=${DLVA_HOME}/gpus/target/gpus-${VERSION}.jar:${DLVA_HOME}/somatic/target/somatic-${VERSION}-bin-native.jar
fi
export SLF4J_CONFIG=${DLVA_HOME}/config/logback.xml