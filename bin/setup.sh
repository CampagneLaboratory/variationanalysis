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
EXECUTION_PLATFORM=`tail -1 ${DISTRIBUTION_DIR}/PROFILE.txt`

if [ -n "${FORCE_PLATFORM+set}" ]; then
    EXECUTION_PLATFORM="${FORCE_PLATFORM}"
fi
export DLVA_JAR=${DLVA_HOME}/gpus/target/gpus-${VERSION}.jar:${DLVA_HOME}/somatic/target/somatic-${VERSION}-bin-${EXECUTION_PLATFORM}.jar
export GDLVA_JAR=${DLVA_HOME}/gpus/target/gpus-${VERSION}.jar:${DLVA_HOME}/genotype/target/genotype-${VERSION}-bin-${EXECUTION_PLATFORM}.jar

export SLF4J_CONFIG=${DLVA_HOME}/config/logback.xml

function assertGobyInstalled {
    goby 1g version >/dev/null 2>&1 || { echo >&2 "This script requires goby but it's not installed. Aborting. Install Goby and add the distribution folder to your path, then try again."; exit 1; }
}