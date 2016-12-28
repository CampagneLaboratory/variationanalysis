# Common functions used in various scripts.

function dieIfError {
    if [ ! $? == 0 ]; then
     echo "An error was encountered ($1)"
     exit;
    fi
}
function assertGobyInstalled {
    goby 1g version >/dev/null 2>&1 || { echo >&2 "This script requires goby but it's not installed. Aborting. Install Goby and add the distribution folder to your path, then try again."; exit 1; }
}

function assertParallelInstalled {
    echo donothing |parallel echo >/dev/null 2>&1 || { echo >&2 "This script requires GNU parallel, but it's not installed. Aborting. Install GNU parallel and try again."; exit 1; }
}

function loadConfigure {

    if [ -e configure.sh ]; then
     echo "Loading configure.sh"
     source configure.sh
    fi

}