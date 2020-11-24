#!/bin/bash

# Usage 
function usage {
	echo "# USAGE: $(basename $0) --input=<RUN_INPUT>,<RUN_INPUT>,... --nbDaysBack=<DAYS_BACK> --genome=<GENOME_FOLDER> ";
	echo "# -i/--input 		INPUT directory where runs are";
	echo "# -t/--nbDaysBack	Analyze run OLDER than X days";
	echo "# -g/--genome		GENOME folder";
	echo "# -h/--help 		HELP option";
	echo "#";
}

####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "i:t:g:h" --long "input:,nbDaysBack:,genome:,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
	:
fi;
PARAM=$@

eval set -- "$ARGS"
while true
do
	case "$1" in
		-i|--input)
			GROUP_INPUT="$2"
			GROUP_INPUT=$(echo $GROUP_INPUT | tr "," " ")
			shift 2 
			;;
		-t|--nbDaysBack)
			DAYS="$2"
			shift 2
			;;
		-g|--genome)
			GENOME="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 0
			;;
		--) shift
			break 
			;;
		*) 	echo "# Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done

if [[ -z $GROUP_INPUT ]]; then
echo "Error : No group found." && echo "" && usage && exit 1;
fi

#Verified if STARK analysis is complete, and if it is, recover the run name
COMPLETE=$(find ${GROUP_INPUT}/*/STARKComplete.txt -mtime -$DAYS 2>/dev/null)
if [[ ! -z $COMPLETE ]]; then
	RUN_NAME=$(cat `echo $COMPLETE` | awk -v OFS="\n" '{print $3}')
fi

#Search if the run hasn't been analysed yet by microservices, and if so, find the corresponding SampleSheet
#In the SampleSheet, search for the microservices on the "Description" line of the Header
for RUN in $RUN_NAME; do
	echo $RUN
	if [[ -z `grep ${RUN} /home1/BAS/grentziv/Listener/Analyzed_runs.txt` ]]; then
		SAMPLESHEET=$(find ${GROUP_INPUT}/${RUN} -iname "*SampleSheet.csv" | paste -sd, - | awk -F "," '{print $1}')
		MICROSERVICE=$(grep "Description" $SAMPLESHEET | awk -F "," '{if($1=="Description"){$1=""; print $0}}')
		BED=$(find ${GROUP_INPUT}/${RUN} -iname "*.bed" | paste -sd, - | awk -F "," '{print $1}')
	fi
	if echo $MICROSERVICE | grep -q "CANOES"; then
		docker run -dti --name=service-canoes -v ${GROUP_INPUT}/${RUN}:${GROUP_INPUT}/${RUN} -v ${GENOME}:${GENOME} canoes:1.1 /bin/bash /app/bin/canoes run -r ${GROUP_INPUT}/${RUN} -s ${SAMPLESHEET} -l ${BED} -e A,F,M -g ${GENOME}/hg19.fa -o ${GROUP_INPUT}/${RUN}/CANOES
		while [ "$( docker container inspect -f '{{.State.Status}}' service-canoes )" == "running" ]; do
			sleep 5m
		done
		docker rm -f service-canoes
	fi
	echo $RUN >> /home1/BAS/grentziv/Listener/Analyzed_runs.txt
done



