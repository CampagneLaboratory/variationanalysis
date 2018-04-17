# Common functions used in various scripts.

function dieIfError {
    if [ ! $? == 0 ]; then
     echo "An error was encountered ($1)"
     exit;
    fi
}

function dieAndExitWithErrorIfError {
    if [ ! $? == 0 ]; then
     echo "An error was encountered ($1)"
     exit 1;
    fi
}
function assertGobyInstalled {
    goby 1g version >/dev/null 2>&1 || { echo >&2 "This script requires goby but it's not installed. Aborting. Install Goby and add the distribution folder to your path, then try again."; exit 1; }
}

function assertParallelInstalled {
    echo donothing |parallel echo >/dev/null 2>&1 || { echo >&2 "This script requires GNU parallel, but it's not installed. Aborting. Install GNU parallel and try again."; exit 1; }
}

function assertBctftoolsInstalled {
  bcftools >/dev/null 2>&1 || if [ $? !=1 ]; then  { echo >&2 "This script requires bcftools but it's not installed. Aborting. Install bcftools and add it to your path, then try again."; exit 1; }
  fi
}

function assertTabixInstalled {
  tabis >/dev/null 2>&1 ||{ echo >&2 "This script requires tabix but it's not installed. Aborting. Install tabix and add it to your path, then try again."; exit 1; }
}

function loadConfigure {

    if [ -e configure.sh ]; then
     echo "Loading configure.sh"
     source configure.sh
    fi

}


export DISTRIBUTION_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
if [[ $OSTYPE == "cygwin" ]]; then
    DISTRIBUTION_DIR=`cygpath -m "${DISTRIBUTION_DIR}"`
fi
export VERSION=`cat ${DISTRIBUTION_DIR}/VERSION.txt`
export DLVA_HOME=${DISTRIBUTION_DIR}

function resetPlatform {
    export EXECUTION_PLATFORM=`tail -1 ${DISTRIBUTION_DIR}/PROFILE.txt`
    if [ -n "${FORCE_PLATFORM+set}" ]; then
         EXECUTION_PLATFORM="${FORCE_PLATFORM}"
    fi
    echo "Running with EXECUTION_PLATFORM=${EXECUTION_PLATFORM}"
    export DLVA_JAR=${DLVA_HOME}/gpus/target/gpus-${VERSION}.jar:${DLVA_HOME}/somatic/target/somatic-${VERSION}-bin-${EXECUTION_PLATFORM}.jar
    export GDLVA_JAR=${DLVA_HOME}/gpus/target/gpus-${VERSION}.jar:${DLVA_HOME}/genotype/target/genotype-${VERSION}-bin-${EXECUTION_PLATFORM}.jar
    export SLF4J_CONFIG=${DLVA_HOME}/config/logback.xml
}
resetPlatform
