#!/bin/bash

# Parameters
FOLDER=$1
[ "$FOLDER" == "" ] && FOLDER="."
PATTERNS=$2
[ "$PATTERNS" == "" ] && PATTERNS="*"
LEVEL_MIN=$3
[ "$LEVEL_MIN" == "" ] && LEVEL_MIN="0"
LEVEL_MAX=$4
[ "$LEVEL_MAX" == "" ] && LEVEL_MAX="0"
OUTPUT=$5
[ "$OUTPUT" == "" ] && OUTPUT=$FOLDER"/index.idx"
OUTPUT_TMP=$6
[ "$OUTPUT_TMP" == "" ] && OUTPUT_TMP=$OUTPUT".tmp"

# Command
(cd $FOLDER; find $PATTERNS -mindepth $LEVEL_MIN -maxdepth $LEVEL_MAX | sort -ru > $OUTPUT_TMP; cp -f $OUTPUT_TMP $OUTPUT; rm -f $OUTPUT_TMP)
