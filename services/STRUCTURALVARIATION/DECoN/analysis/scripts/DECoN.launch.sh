#!/bin/bash

## SOURCES

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source $SCRIPT_DIR/DECoN.function.sh


## PARAM

INPUT_FOLDER=$1
[ "$INPUT_FOLDER" == "" ] && echo "#[ERROR] NO input folder" && exit 0

DECoN_main_results=$2
[ "$DECoN_main_results" == "" ] && DECoN_main_results=/STARK/data/DECoN/results && echo "#[WARNING] NO ooutput folder. Default /STARK/data/DECoN/results"

ANALYSIS_TYPE=$3
[ "$ANALYSIS_TYPE" == "" ] && ANALYSIS_TYPE=RUN && echo "#[WARNING] NO analysis type. Default RUN"

VERBOSE=$4

ANALYSIS_FUNCTION="DECON_"$ANALYSIS_TYPE

$ANALYSIS_FUNCTION $INPUT_FOLDER



## datasets as run in repository
if false; then


    # RUNS
    if [ "$ANALYSIS_TYPE" == "RUN" ]; then
        (($VERBOSE)) && echo "#[INFO] Analyses by RUN"
        for INPUT_FOLDER in $(find $REPOSITORY_FOLDER/ -mindepth 3 -maxdepth 3 -type d); do 
            time DECON_RUN $INPUT_FOLDER
        done;
    fi;

    # PROJECTS
    if [ "$ANALYSIS_TYPE" == "PROJECT" ]; then
        (($VERBOSE)) && echo "#[INFO] Analyses by PROJECTS"
        for INPUT_FOLDER in $(find $REPOSITORY_FOLDER/ -mindepth 2 -maxdepth 2 -type d); do
            time DECON_PROJECT $INPUT_FOLDER
        done;
    fi;

    # GROUP
    if [ "$ANALYSIS_TYPE" == "GROUP" ]; then
        (($VERBOSE)) && echo "#[INFO] Analyses by GROUP"
       for INPUT_FOLDER in $(find $REPOSITORY_FOLDER/ -mindepth 1 -maxdepth 1 -type d); do
            time DECON_GROUP $INPUT_FOLDER
        done;
    fi;

    # REPOSITORY
    if [ "$ANALYSIS_TYPE" == "REPOSITORY" ]; then
        (($VERBOSE)) && echo "#[INFO] Analyses by REPOSITORY"
        for INPUT_FOLDER in $(find $REPOSITORY_FOLDER/ -mindepth 0 -maxdepth 0 -type d); do
            time DECON_REPOSITORY $INPUT_FOLDER
        done;
    fi;

    # DATASETS
    if [ "$ANALYSIS_TYPE" == "DATASE" ]; then
        (($VERBOSE)) && echo "#[INFO] Analyses by DATASETS [$DATASETS_FOLDER]"
        for INPUT_FOLDER in $(find $DATASETS_FOLDER/ -mindepth 1 -maxdepth 1 -type d); do 
            time DECON_DATASET $INPUT_FOLDER
        done;
    fi;


fi;
