


## SOURCES

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source $SCRIPT_DIR/DECoN.function.sh


## PARAM

THREADS=$1
[ "$THREADS" == "" ] && THREADS=1


INSTALL $THREADS
