#!/usr/bin/env bash
if [ $# -lt 1 ]; then
  echo "The first argument must be memory (i.e., 4g)"
  exit 1
fi
memory_requirement=$1
shift
other_parameters=$*

export DISTRIBUTION_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. && pwd )"
if [[ $OSTYPE == "cygwin" ]]; then
    DISTRIBUTION_DIR=`cygpath -m "${DISTRIBUTION_DIR}"`
fi

. ${DISTRIBUTION_DIR}/bin/common.sh
