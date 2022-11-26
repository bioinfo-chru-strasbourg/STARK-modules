#!/bin/bash

# Parameters
FOLDER=$1
[ "$FOLDER" == "" ] && FOLDER="."
PATTERNS=$2
[ "$PATTERNS" == "" ] && PATTERNS="*"
FILES_PATTERNS=$3
[ "$FILES_PATTERNS" == "" ] && FILES_PATTERNS=""
LEVEL_MIN=$4
[ "$LEVEL_MIN" == "" ] && LEVEL_MIN="0"
LEVEL_MAX=$5
[ "$LEVEL_MAX" == "" ] && LEVEL_MAX="0"
OUTPUT=$6
[ "$OUTPUT" == "" ] && OUTPUT=$FOLDER"/index.idx"
OUTPUT_TMP=$7
[ "$OUTPUT_TMP" == "" ] && OUTPUT_TMP=$OUTPUT".tmp"

# Command
echo "cd $FOLDER; find $PATTERNS -mindepth $LEVEL_MIN -maxdepth $LEVEL_MAX $FILES_PATTERNS | sort -ru | xargs ls -t"
#echo "(cd $FOLDER; find $PATTERNS -mindepth $LEVEL_MIN -maxdepth $LEVEL_MAX $FILES_PATTERNS | sort -ru | xargs ls -t > $OUTPUT_TMP; cp -f $OUTPUT_TMP $OUTPUT; rm -f $OUTPUT_TMP)"
(cd $FOLDER; find $PATTERNS -mindepth $LEVEL_MIN -maxdepth $LEVEL_MAX $FILES_PATTERNS | sort -ru | xargs ls -t > $OUTPUT_TMP; cp -f $OUTPUT_TMP $OUTPUT; rm -f $OUTPUT_TMP)
#(cd $FOLDER; find $PATTERNS -mindepth $LEVEL_MIN -maxdepth $LEVEL_MAX $FILES_PATTERNS -print0 | xargs -0 sort -ru | xargs ls -t > $OUTPUT_TMP; cp -f $OUTPUT_TMP $OUTPUT; rm -f $OUTPUT_TMP)
#(cd $FOLDER; find $PATTERNS -mindepth $LEVEL_MIN -maxdepth $LEVEL_MAX $FILES_PATTERNS -exec ls -l {} "+" | sort -ru > $OUTPUT_TMP; cp -f $OUTPUT_TMP $OUTPUT; rm -f $OUTPUT_TMP)
