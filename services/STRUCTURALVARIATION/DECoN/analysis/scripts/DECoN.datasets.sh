

## SOURCES

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source $SCRIPT_DIR/DECoN.function.sh


## PARAM

DATASETS_FOLDER=$1
[ "$DATASETS_FOLDER" == "" ] && DATASETS_FOLDER=/STARK/data/DECoN/DATASETS

THREADS=$2
[ "$THREADS" == "" ] && THREADS=1

VERBOSE=1


## Datasets list by folders
if true; then

    # DATASETS
    if true; then
        (($VERBOSE)) && echo "#[INFO] Analyses by DATASETS [$DATASETS_FOLDER]"
        for INPUT_FOLDER in $(find $DATASETS_FOLDER/ -mindepth 1 -maxdepth 1 -type d); do 
            time DECON_DATASET $INPUT_FOLDER
        done;
    fi;

fi;

