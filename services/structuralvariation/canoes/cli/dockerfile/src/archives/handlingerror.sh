#!/bin/bash
#################################
##
## Handling CANOES errors
##
#################################

#NB this two path DONOT be the last var in the conf.json !!!!!!!!!!!!

#get Run path in conf.json which regroup run analysis configuration
#echo "WARNING Handling canoes issues and copy logfiles in samples"
#while read -r line; do
#	#echo "#[INFO] $line"
#	if grep -q "runPath"* <<< "$line"; then
#		IFS=':' read -ra PARSE <<< "$line"
#		for VAR in "${PARSE[@]}"; do
#			if grep -q "runPath" <<< "$line"; then
#				RUN_TMP=$VAR
#				RUN=${RUN_TMP:1:-2}
#			fi;
#		done;
#	elif grep -q "repositoryPath"* <<< "$line"; then
#		IFS=':' read -ra PARSE <<< "$line"
#		for VAR in "${PARSE[@]}"; do
#			if grep -q "repositoryPath" <<< "$line"; then
#				REPO_TMP=$VAR
#				REPO=${REPO_TMP:1:-2}
#			fi;
#		done;
#	elif grep -q "outputPath"* <<< "$line"; then
#		IFS=':' read -ra PARSE <<< "$line"
#		for VAR in "${PARSE[@]}"; do
#			if grep -q "outputPath" <<< "$line"; then
#				OUTPUT_TMP=$VAR
#				OUTPUT=${OUTPUT_TMP:1:-2}
#			fi;
#		done;
#	fi;
#
#done < "conf.json"

#Copy logs in folder if errors in snakemake script
echo "#[INFO] CANOES repository $1"
echo "#[INFO] CANOES output $2"
#echo "#[INFO] RUN folder $3"

REPO=$1
OUTPUT=$2

if [ ! -z "$(ls -A $OUTPUT)" ]; then
	if [[ $OUTPUT == "/STARK/output"* ]]; then
		echo "#[INFO] Nothing todo output and repository are the same folder"
	else
		for SAMPLE in $(find $OUTPUT -maxdepth 1 -mindepth 1 -type d ! -name "CANOES" ! -name "logs" ! -name ".*"); do
			#if sample is in repository
			if [ -d "$REPO/$(basename $SAMPLE)" ]; then
				echo "#[INFO] Sample handling error $SAMPLE"
				if [[ -d "$SAMPLE" ]] && [[ -d "$REPO/$(basename $SAMPLE)/CANOES" ]]; then
					echo "WARNING $REPO/$(basename $SAMPLE)/CANOES already exists !"
					mkdir -p "$REPO/$(basename $SAMPLE)/CANOES/logs" "$REPO/$(basename $SAMPLE)/CANOES/logs/snakemake" && rsync -ar --force "$OUTPUT/logs" "$REPO/$(basename $SAMPLE)/CANOES/"
					rsync -ar --force "$OUTPUT/.snakemake/log" "$REPO/$(basename $SAMPLE)/CANOES/logs/snakemake/"
				elif [[ -d "$SAMPLE" ]] && [[ ! -d "$REPO/$(basename $SAMPLE)/CANOES" ]]; then
					echo "#[INFO] Copy logs for sample $(basename $SAMPLE)"
					mkdir -p "$REPO/$(basename $SAMPLE)/CANOES" "$REPO/$(basename $SAMPLE)/CANOES/logs/snakemake" && rsync -ar "$OUTPUT/logs" "$REPO/$(basename $SAMPLE)/CANOES/"
					rsync -ar --force "$OUTPUT/.snakemake/log" "$REPO/$(basename $SAMPLE)/CANOES/logs/snakemake/"
				else
					echo "ERROR $SAMPLE does not exists !"
				fi;
			else
				echo "#[INFO] skip $SAMPLE, probably comming from other run in project"
			fi;
		done;
	fi;
else
	echo "WARNING nothing to do, no results generated"
fi;
