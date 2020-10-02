


## SOURCES

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source $SCRIPT_DIR/DECoN.function.sh


## PARAM

REPOSITORY_FOLDER=$1
[ "$REPOSITORY_FOLDER" == "" ] && REPOSITORY_FOLDER=/STARK/data/DECoN/Repository

THREADS=$2
[ "$THREADS" == "" ] && THREADS=1

VERBOSE=1


## datasets as run in repository
if true; then

    VERBOSE=1

    # RUNS
    if true; then
        (($VERBOSE)) && echo "#[INFO] Analyses by RUN"
        for INPUT_FOLDER in $(find $REPOSITORY_FOLDER/ -mindepth 3 -maxdepth 3 -type d); do 
            time DECON_RUN $INPUT_FOLDER
        done;
    fi;

    # PROJECTS
    if true; then
        (($VERBOSE)) && echo "#[INFO] Analyses by PROJECTS"
        for INPUT_FOLDER in $(find $REPOSITORY_FOLDER/ -mindepth 2 -maxdepth 2 -type d); do
            time DECON_PROJECT $INPUT_FOLDER
        done;
    fi;

    # GROUP
    if true; then
        (($VERBOSE)) && echo "#[INFO] Analyses by GROUP"
       for INPUT_FOLDER in $(find $REPOSITORY_FOLDER/ -mindepth 1 -maxdepth 1 -type d); do
            time DECON_GROUP $INPUT_FOLDER
        done;
    fi;

    # REPOSITORY
    if true; then
        (($VERBOSE)) && echo "#[INFO] Analyses by REPOSITORY"
        for INPUT_FOLDER in $(find $REPOSITORY_FOLDER/ -mindepth 0 -maxdepth 0 -type d); do
            time DECON_REPOSITORY $INPUT_FOLDER
        done;
    fi;


fi;
