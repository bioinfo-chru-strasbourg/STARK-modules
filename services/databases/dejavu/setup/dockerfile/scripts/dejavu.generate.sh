#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARK_DEJAVU"
SCRIPT_DESCRIPTION="STARK DEJAVU ANNOVAR databases generation"
SCRIPT_RELEASE="0.12.1"
SCRIPT_DATE="06/01/2021"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS"
SCRIPT_LICENCE="GNU-AGPL"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-05/09/2017: Script creation\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.1b-07/09/2017: Add generation of ANNOVAR generic database.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.9.2b-02/11/2018: Use BCFTOOLS instead of VCFTOOLS.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.10.0-12/08/2020: Many changes.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.11.0-29/09/2020: Add STARK module json files.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.12.0-17/12/2020: Add sample filter parameters.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.12.1-06/01/2021: Add group project filter, annotation, calculation and nomen fields parameters, rebuld search goup/project folders.\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)

source $ENV_CONFIG 1>/dev/null 2>/dev/null

ENV_TOOLS=$(find -L $SCRIPT_DIR/.. -name toold.app)

source $ENV_TOOLS 1>/dev/null 2>/dev/null


# Header
function header () {
	echo "#######################################";
	echo "# $SCRIPT_NAME [$SCRIPT_RELEASE-$SCRIPT_DATE]";
	echo "# $SCRIPT_DESCRIPTION ";
	echo "# $SCRIPT_AUTHOR @ $SCRIPT_COPYRIGHT © $SCRIPT_LICENCE";
	echo "#######################################";
}

# Release
function release () {
	echo "# RELEASE NOTES:";
	echo -e $RELEASE_NOTES
}

# Usage
function usage {
	echo "# USAGE: $(basename $0) [options...]";
	echo "# --application=<STRING|FILE>                   APP name or APP file configuration of the APPLICATION.";
	echo "#                                               Must be in the STARK APPS folder if relative path";
	echo "#                                               Default: 'default.app'";
	echo "# --app_folder|--application_folder             STARK Application folder";
	echo "#                                               Used to detect STARK Repository folders";
	echo "#                                               Default: default STARK Application folder";
	echo "# -r|--repo_folder|--repository_folder          STARK Repository folders";
	echo "#                                               List of STARK Repository folders containing group and project results";
	echo "#                                               Format: 'folder1[,folder2...]'";
	echo "#                                               Default: default STARK Repository folder";
	echo "# --dejavu_folder                               Output DejaVu folder";
	echo "#                                               Default: '.'";
	echo "# --dejavu_release                              Output DejaVu release";
	echo "#                                               Default: 'date +%Y%m%d-%H%M%S'";
	echo "# --dejavu_annotation                           Output VCF DejaVu annotation";
	echo "#                                               Default: 'HOWARD_ANNOTATION_REPORT' STARK parameter";
	echo "#                                               Example: 'ALL,snpeff' for all annotations";
	echo "#                                               Tip: use 'none' for no annotation, 'ALL,snpeff,snpeff_hgvs' for ALL annotation";
	echo "# --dejavu_calculation                          Output VCF DejaVu calculation";
	echo "#                                               Default: 'HOWARD_CALCULATION_REPORT' STARK parameter";
	echo "#                                               Example: 'VAF,VAF_STATS,DP_STATS,VARTYPE,NOMEN'";
	echo "#                                               Tip: use 'none' for no calculation";
	echo "# --dejavu_nomen_fields                         Output VCF DejaVu NOMEN field";
	echo "#                                               Default: 'HOWARD_NOMEN_FIELDS' STARK parameter";
	echo "#                                               Example: 'hgvs', 'snpeff_hgvs'";
	echo "# --sample_exclude                              Exclude sample pattern (regexp)";
	echo "#                                               Format: '<group>/<project>/<sample_pattern>[,<group>/<project>/<run>/<sample_pattern>]'";
	echo "#                                               Example: 'GENOME/GERMLINE/.*/.*CORIEL.*' to exclude all *CORIEL* samples";
	echo "#                                               Default: ''";
	echo "# --sample_exclude_file                         Exclude sample pattern (regexp) within a file";
	echo "#                                               Format: same as --sample_exclude parameter";
	echo "#                                               Default: <STARK_FOLDER_CONFIG>/dejavu/sample_exclude.conf'";
	echo "# --group_project_list                          Include only group/project list (shell like)";
	echo "#                                               These folders must be well structured as group/project within repository folders";
	echo "#                                               Format: '<group>/<project>[,<group>/<project>]'";
	echo "#                                               Example: 'GENOME/*' to include only run from group GENOME and all project";
	echo "#                                               Default: '' (empty), filter will be used";
	echo "# --group_project_list_file                     Include only group/project/run list (shell like) within a file";
	echo "#                                               Format: same as --group_project_list parameter";
	echo "#                                               Default: <STARK_FOLDER_CONFIG>/dejavu/group_project_list.conf'";
	echo "# --group_project_filter                        Include only group/project/run pattern (shell like)";
	echo "#                                               These patterns will be used to detect well structured folders, as group/project, within repository folders";
	echo "#                                               Format: '<group>/<project>/<run>[,<group>/<project>/<run>]'";
	echo "#                                               Example: 'GENOME/*/19*' to include only run of year 2019 from group GENOME and all project";
	echo "#                                               Default: '*/*/*', all groups, projects and runs";
	echo "# --group_project_filter_file                   Include only group/project/run pattern (shell like) within a file";
	echo "#                                               Format: same as --group_project_filter parameter";
	echo "#                                               Default: <STARK_FOLDER_CONFIG>/dejavu/group_project_filter.conf'";
	echo "# --tmp                                         Temporary folder";
	echo "#                                               Default: default STARK Temporary folder";
	echo "# --bcftools                                    BCFTools application binary";
	echo "#                                               Default: default STARK configuration or 'bcftools'";
	echo "# --tabix                                       TABix application binary";
	echo "#                                               Default: default STARK configuration or 'tabix'";
	echo "# --bgzip                                       BGZip application binary";
	echo "#                                               Default: default STARK configuration or 'bgzip'";
	echo "# --annovar                                     ANNOVAR application binary folder";
	echo "#                                               Default: default STARK configuration or ''";
	echo "# --verbose                                     VERBOSE";
	echo "# --debug                                       DEBUG";
	echo "# --release                                     RELEASE";
	echo "# --help                                        HELP";
	echo "#";
}

# header
header;


####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "e:r:vdnh" --long "env:,app:,application:,app_folder:,application_folder:,repo_folder:,repository_folder:,dejavu_folder:,dejavu_release:,dejavu_annotation:,dejavu_calculation:,dejavu_nomen_fields:,sample_exclude:,sample_exclude_file:,group_project_list:,group_project_list_file:,group_project_filter:,group_project_filter_file:,tmp:,bcftools:,tabix:,bgzip:,annovar:,verbose,debug,release,help" -- "$@" 2> /dev/null)

eval set -- "$ARGS"
while true
do
	case "$1" in
		-e|--env|--app|--application)
			APP="$2"
			shift 2
			;;
		--app_folder|--application_folder)
			APP_FOLDER="$2";
			shift 2
			;;
		-r|--repo_folder|--repository_folder)
			REPO_FOLDER=$(echo "$2" | tr "," " ");
			shift 2
			;;
		--dejavu_folder)
			DEJAVU_FOLDER="$2";
			shift 2
			;;
		--dejavu_release)
			RELEASE="$2";
			shift 2
			;;
		--dejavu_annotation)
			DEJAVU_ANNOTATION="$2";
			shift 2
			;;
		--dejavu_calculation)
			DEJAVU_CALCULATION="$2";
			shift 2
			;;
		--dejavu_nomen_fields)
			DEJAVU_NOMEN_FIELDS="$2";
			shift 2
			;;
		--sample_exclude)
			SAMPLE_EXCLUDE=$(echo "$2" | tr "," " ");
			shift 2
			;;
		--sample_exclude_file)
			SAMPLE_EXCLUDE_FILE=$2;
			shift 2
			;;
		--group_project_list)
			GROUP_PROJECT_LIST=$(echo "$2" | tr "," " ");
			shift 2
			;;
		--group_project_list_file)
			GROUP_PROJECT_LIST_FILE=$2;
			shift 2
			;;
		--group_project_filter)
			GROUP_PROJECT_FILTER=$(echo "$2" | tr "," " ");
			shift 2
			;;
		--group_project_filter_file)
			GROUP_PROJECT_FILTER_FILE=$2;
			shift 2
			;;
		--tmp)
			TMP_FOLDER_TMP="$2";
			shift 2
			;;
		--bcftools)
			BCFTOOLS="$2";
			shift 2
			;;
		--tabix)
			TABIX="$2";
			shift 2
			;;
		--bgzip)
			BGZIP="$2";
			shift 2
			;;
		--annovar)
			ANNOVAR="$2";
			shift 2
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-d|--debug)
			VERBOSE=1
			DEBUG=1
			shift 1
			;;
		-n|--release)
			release;
			exit 0
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


## PARAMETERS
##############




# ENV
#########

#echo "APP=$APP"; exit;
(($VERBOSE)) && [ ! -z "$APP" ] && echo "#[INFO] Search Application '$APP'"

ENV=$(find_app "$APP" "$STARK_FOLDER_APPS")
source_app "$APP" "$STARK_FOLDER_APPS" 1

export ENV
export APP

(($VERBOSE)) && [ ! -z "$APP" ] && [ ! -z "$ENV" ] && echo "#[INFO] Application '$APP' found ('$ENV')"
(($VERBOSE)) && [ ! -z "$APP" ] && [ -z "$ENV" ] && echo "#[INFO] Application '$APP' NOT found"



# APP FOLDER
if [ ! -z "$APP_FOLDER" ] && [ -d "$APP_FOLDER" ]; then
    STARK_FOLDER_APPS=$APP_FOLDER
fi;


# REPO_FOLDER
REPO_FOLDER=$(echo $REPO_FOLDER | tr "," " ")

# FOLDER DEJAVU
if [ -z "$DEJAVU_FOLDER" ]; then
    DEJAVU_FOLDER="."
fi;

# BCFTOOLS
if [ -z "$BCFTOOLS" ]; then
    BCFTOOLS="bcftools"
fi;

# TABIX
if [ -z "$TABIX" ]; then
    TABIX="tabix"
fi;

# BGZIP
if [ -z "$BGZIP" ]; then
    BGZIP="bgzip"
fi;

# VCFSTATS
if [ -z "$VCFSTATS" ]; then
	VCFSTATS=$(whereis VcfStats.jar | awk -F" " '{print $2}')
fi;


# DEJAVU ANNOTATION
#DEJAVU_ANNOTATION=$HOWARD_ANNOTATION_REPORT
if [ -z "$DEJAVU_ANNOTATION" ]; then
    DEJAVU_ANNOTATION=$HOWARD_ANNOTATION_REPORT
fi;

# DEJAVU CALCULATION
if [ -z "$DEJAVU_CALCULATION" ]; then
    DEJAVU_CALCULATION=$HOWARD_CALCULATION_REPORT
fi;

# DEJAVU NOMEN FIELDS
if [ -z "$DEJAVU_NOMEN_FIELDS" ]; then
    DEJAVU_NOMEN_FIELDS=$HOWARD_NOMEN_FIELDS
fi;




if [ -z "$RELEASE" ]; then
	RELEASE=$(date +%Y%m%d-%H%M%S)
fi;



# DEJAVU
DEJAVU=$DEJAVU_FOLDER
mkdir -p $DEJAVU/$RELEASE

DEJAVU_SUBFOLDER_ANNOVAR=annovar
DEJAVU_FOLDER_ANNOVAR=$DEJAVU/$RELEASE/$DEJAVU_SUBFOLDER_ANNOVAR
mkdir -p $DEJAVU_FOLDER_ANNOVAR

DEJAVU_SUBFOLDER_VCF=vcf
DEJAVU_FOLDER_VCF=$DEJAVU/$RELEASE/$DEJAVU_SUBFOLDER_VCF
mkdir -p $DEJAVU_FOLDER_VCF

DEJAVU_SUBFOLDER_LOG=log
DEJAVU_FOLDER_LOG=$DEJAVU/$RELEASE/$DEJAVU_SUBFOLDER_LOG
mkdir -p $DEJAVU_FOLDER_LOG

DEJAVU_SUBFOLDER_STATS=stats
#DEJAVU_FOLDER_STATS=$DEJAVU/$RELEASE/$DEJAVU_SUBFOLDER_STATS
#mkdir -p $DEJAVU_FOLDER_STATS

# SAMPLE_EXCLUDE_FILE
if [ ! -e $SAMPLE_EXCLUDE_FILE ] || [ "$SAMPLE_EXCLUDE_FILE" == "" ]; then
	SAMPLE_EXCLUDE_FILE="$STARK_FOLDER_CONFIG/dejavu/sample_exclude.conf"
fi
# Load sample exclude patterns
if [ -e $SAMPLE_EXCLUDE_FILE ] && [ "$SAMPLE_EXCLUDE_FILE" != "" ]; then
	SAMPLE_EXCLUDE=$SAMPLE_EXCLUDE" "$(cat $SAMPLE_EXCLUDE_FILE | grep -v "^[ \t]*#" | tr "\n" " ")
fi;



# GROUP_PROJECT_LIST
if [ ! -e $GROUP_PROJECT_LIST_FILE ] || [ "$GROUP_PROJECT_LIST_FILE" == "" ]; then
	GROUP_PROJECT_LIST_FILE="$STARK_FOLDER_CONFIG/dejavu/group_project_list.conf"
fi
# Load GROUP_PROJECT_LIST patterns
if [ -e $GROUP_PROJECT_LIST_FILE ] && [ "$GROUP_PROJECT_LIST_FILE" != "" ]; then
	GROUP_PROJECT_LIST=$GROUP_PROJECT_LIST" "$(cat $GROUP_PROJECT_LIST_FILE | grep -v "^[ \t]*#" | tr "\n" " ")
fi;

if [ "$GROUP_PROJECT_LIST" == "" ] || [ "$GROUP_PROJECT_LIST" == " " ]; then
	GROUP_PROJECT_LIST=""
fi;


# GROUP_PROJECT_FILTER
if [ ! -e $GROUP_PROJECT_FILTER_FILE ] || [ "$GROUP_PROJECT_FILTER_FILE" == "" ]; then
	GROUP_PROJECT_FILTER_FILE="$STARK_FOLDER_CONFIG/dejavu/group_project_filter.conf"
fi
# Load GROUP_PROJECT_FILTER patterns
if [ -e $GROUP_PROJECT_FILTER_FILE ] && [ "$GROUP_PROJECT_FILTER_FILE" != "" ]; then
	GROUP_PROJECT_FILTER=$GROUP_PROJECT_FILTER" "$(cat $GROUP_PROJECT_FILTER_FILE | grep -v "^[ \t]*#" | tr "\n" " ")
fi;

if [ "$GROUP_PROJECT_FILTER" == "" ] || [ "$GROUP_PROJECT_FILTER" == " " ]; then
	GROUP_PROJECT_FILTER="*/*/*"
fi;

# DEJAVU FOLDER
#mkdir -p $DEJAVU_FOLDER_LOG
# mkdir -p $DEJAVU/$RELEASE/annovar
# mkdir -p $DEJAVU/$RELEASE/vcf
# mkdir -p $DEJAVU/$RELEASE/log


# ASSEMBLY PREFIX

ASSEMBLY_PREFIX_DEFAULT="hg19"


# VCF pattern

VCF_PATTERN=".final.vcf.gz"


# INFO SUFFIX (e.g. release): $PREFIX_dejavu.$GROUP.$PROJECT$SUFFIX.txt
#SUFFIX=".$RELEASE"
SUFFIX=""

# TMP
if [ "$TMP_FOLDER_TMP" == "" ]; then TMP_FOLDER_TMP=/tmp; fi;
TMP=$TMP_FOLDER_TMP/dejavu_$RANDOM$RANDOM$RANDOM$RANDOM
mkdir -p $TMP


# MK
MK=$DEJAVU_FOLDER_LOG/$RELEASE.mk
> $MK

# LOG
LOG=$DEJAVU_FOLDER_LOG/$RELEASE.log
> $LOG


#### INFO
(($VERBOSE)) && echo "#[INFO] RELEASE: $RELEASE"
(($VERBOSE)) && echo "#[INFO] STARK FOLDER APPLICATIONS: $STARK_FOLDER_APPS"
(($VERBOSE)) && echo "#[INFO] DEJAVU FOLDER: $DEJAVU_FOLDER"

(($VERBOSE)) && echo "#[INFO] REPOSITORY FOLDER: "
(($VERBOSE)) && for RF in $REPO_FOLDER; do echo "#[INFO]    "$RF; done

(($VERBOSE)) && echo "#[INFO] GROUP/PROJECT LIST FILE:"
(($VERBOSE)) && echo "#[INFO]    $GROUP_PROJECT_LIST_FILE"
(($VERBOSE)) && echo "#[INFO] GROUP/PROJECT LIST:"
(($VERBOSE)) && for GPL in $GROUP_PROJECT_LIST; do echo "#[INFO]    "$GPL; done

(($VERBOSE)) && echo "#[INFO] GROUP/PROJECT FILTER FILE:"
(($VERBOSE)) && echo "#[INFO]    $GROUP_PROJECT_FILTER_FILE"
(($VERBOSE)) && echo "#[INFO] GROUP/PROJECT FILTER:"
(($VERBOSE)) && for GPF in $GROUP_PROJECT_FILTER; do echo "#[INFO]    "$GPF; done

(($VERBOSE)) && echo "#[INFO] SAMPLE EXCLUDE FILE:"
(($VERBOSE)) && echo "#[INFO]    $SAMPLE_EXCLUDE_FILE"
(($VERBOSE)) && echo "#[INFO] SAMPLE EXCLUDE:"
(($VERBOSE)) && for SE in $SAMPLE_EXCLUDE; do echo "#[INFO]    "$SE; done

(($VERBOSE)) && echo "#[INFO] DEJAVU ANNOTATION: $DEJAVU_ANNOTATION"
(($VERBOSE)) && echo "#[INFO] DEJAVU CALCULATION: $DEJAVU_CALCULATION"
(($VERBOSE)) && echo "#[INFO] DEJAVU NOMEN FIELDS: $DEJAVU_NOMEN_FIELDS"

(($VERBOSE)) && echo "#[INFO] THREADS: $THREADS"


(($DEBUG)) && echo "#[INFO] TMP: $TMP"
(($DEBUG)) && echo "#[INFO] BCFTOOLS: $BCFTOOLS"
(($DEBUG)) && echo "#[INFO] BGZIP: $BGZIP"
(($DEBUG)) && echo "#[INFO] TABIX: $TABIX"
(($DEBUG)) && echo "#[INFO] VCFSTATS: $VCFSTATS"



### Find Group folders
########################


GP_FOLDER_LIST=""


# FOLDER=$1
# [ "$FOLDER" == "" ] && FOLDER="."
# PATTERNS=$2
# [ "$PATTERNS" == "" ] && PATTERNS="*"
# FILES_PATTERNS=$3
# [ "$FILES_PATTERNS" == "" ] && FILES_PATTERNS=""
# LEVEL_MIN=$4
# [ "$LEVEL_MIN" == "" ] && LEVEL_MIN="0"
# LEVEL_MAX=$5
# [ "$LEVEL_MAX" == "" ] && LEVEL_MAX="0"
# OUTPUT=$6
# [ "$OUTPUT" == "" ] && OUTPUT=$FOLDER"/index.idx"
# OUTPUT_TMP=$7
# [ "$OUTPUT_TMP" == "" ] && OUTPUT_TMP=$OUTPUT".tmp"

# (cd $FOLDER; find $PATTERNS -mindepth $LEVEL_MIN -maxdepth $LEVEL_MAX $FILES_PATTERNS | sort -ru | xargs ls -t > $OUTPUT_TMP; cp -f $OUTPUT_TMP $OUTPUT; rm -f $OUTPUT_TMP)

(($VERBOSE)) && echo "#"
#(($VERBOSE)) && echo "#[INFO] DEJAVU database repository/group/project detection"
echo "#[INFO] DEJAVU database repository/group/project detection"



if [ ! -z "$REPO_FOLDER" ]; then

	#GP_FOLDER_LIST=$GROUP_PROJECT_LIST

	# Group Project List
	if [ ! -z "$GROUP_PROJECT_LIST" ]; then
		(($VERBOSE)) && echo "#[INFO] DEJAVU database repository/group/project detection from GROUP/PROJECT list"
		GP_FOLDER_LIST="";
		for RF in $REPO_FOLDER; do
			for GPR_LIST in $GROUP_PROJECT_LIST; do
				GP_FOLDER_LIST="$GP_FOLDER_LIST	$RF/$GPR_LIST"
			done;
		done
	fi;

	#echo "GP_FOLDER_LIST=$GP_FOLDER_LIST";

	if [ -z "$GP_FOLDER_LIST" ]; then

		(($VERBOSE)) && echo "#[INFO] DEJAVU database repository/group/project detection from GROUP/PROJECT filter"

		# Index method
		if true; then

			# Index file
			INDEX=$TMP/index.idx

			# TODO
			GPR_FILTERS=$GROUP_PROJECT_FILTER		
			#GPR_FILTERS="*/*/* */*/19*"


			# Repository filter
			REPOSITORIES_FILTER="";
			for RF in $REPO_FOLDER; do
				for GPR_FILTER in $GPR_FILTERS; do
					REPOSITORIES_FILTER="$REPOSITORIES_FILTER $RF/$GPR_FILTER"
				done;
			done

			#(($VERBOSE)) && echo "#[INFO] DEJAVU database repository/group/project filter"
			#(($VERBOSE)) && echo "#[INFO]    $REPOSITORIES_FILTER"

			# Command param
			LEVEL_MIN=2
			LEVEL_MAX=2
			FILES_PATTERNS=" -name *$VCF_PATTERN "
			# Command
			#CMD="find $REPOSITORIES_FILTER -mindepth $LEVEL_MIN -maxdepth $LEVEL_MAX $FILES_PATTERNS | sort -u > $INDEX"
			#eval $CMD

			if [ -z "$GP_FOLDER_LIST" ]; then
				find $REPOSITORIES_FILTER -mindepth $LEVEL_MIN -maxdepth $LEVEL_MAX $FILES_PATTERNS | sort -u > $INDEX
				GP_FOLDER_LIST=$(cat $INDEX | grep -v " " | xargs dirname | xargs dirname  | xargs dirname | sort -u)
			else
				GP_FOLDER_LIST=""
			fi;

		# ls method
		else
			if [ ! -z "$REPO_FOLDER" ]; then
				#GP_FOLDER_LIST=$(find -L $REPO_FOLDER -maxdepth 2 -mindepth 2 -type d 2>/dev/null)
				GP_FOLDER_LIST=$(ls $(ls $REPO_FOLDER/*/*/*/STARKCopyComplete.txt 2>/dev/null | xargs dirname | xargs dirname | sort -u | sed 's#$#/*/*/*'$VCF_PATTERN'#g') 2>/dev/null | xargs dirname | xargs dirname | xargs dirname | sort -u)
			fi;
		fi;

	fi;

fi;


if [ -z "$GP_FOLDER_LIST" ]; then

	(($VERBOSE)) && echo "#[INFO] DEJAVU database repository/group/project detection from applications"

    for ENV_DEF in $(find -L $STARK_FOLDER_APPS -name '*.app' -type f | sed s#$STARK_FOLDER_APPS/## | sort -f -t'/' -k2.3 -k2.2 -k2.1) $(find -L $STARK_FOLDER_APPS -name '*.plugapp' -type f | sed s#$STARK_FOLDER_APPS/## | sort -f -t'/' -k2.3 -k2.2 -k2.1); do
    
    	# APP INFO
        APP_FOLDER_ARCHIVES=$(source_app "$ENV_DEF"  2>/dev/null; echo $FOLDER_ARCHIVES)
        APP_FOLDER_REPOSITORY=$(source_app "$ENV_DEF"  2>/dev/null; echo $FOLDER_REPOSITORY)
        APP_GROUP=$(source_app "$ENV_DEF"  2>/dev/null; echo $APP_GROUP)
        APP_PROJECT=$(source_app "$ENV_DEF"  2>/dev/null; echo $APP_PROJECT)

        # REPOSITORY FOLDER
        if [ ! -z "$REPO_FOLDER" ] && [ -d "$REPO_FOLDER" ] && [ -d "$REPO_FOLDER/$APP_GROUP/$APP_PROJECT" ]; then
        	APP_FOLDER_ARCHIVES=$REPO_FOLDER
        fi;

        # GROUP FOLDER
        if [ -z "$APP_GROUP" ]; then
            APP_GROUP="UNKNOWN"
        fi;

        # PROJECT FOLDER
        if [ -z "$APP_PROJECT" ]; then
            APP_PROJECT="UNKNOWN"
        fi;
        
        # Add repository group project folder 
        if [ ! -z "$APP_FOLDER_ARCHIVES" ] && [ ! -z "$APP_GROUP" ] && [ ! -z "$APP_PROJECT" ] && [ -d "$APP_FOLDER_ARCHIVES/$APP_GROUP/$APP_PROJECT" ]; then
	        (($VERBOSE)) && echo "#[INFO] DEJAVU database '$APP_GROUP/$APP_PROJECT' repository '$APP_FOLDER_ARCHIVES' found"
	        GP_FOLDER_LIST="$GP_FOLDER_LIST\n$APP_FOLDER_ARCHIVES/$APP_GROUP/$APP_PROJECT"
	    fi;
	    if [ ! -z "$APP_FOLDER_REPOSITORY" ] && [ ! -z "$APP_GROUP" ] && [ ! -z "$APP_PROJECT" ] && [ -d "$APP_FOLDER_REPOSITORY/$APP_GROUP/$APP_PROJECT" ]; then
	        (($VERBOSE)) && echo "#[INFO] DEJAVU database '$APP_GROUP/$APP_PROJECT' repository '$APP_FOLDER_REPOSITORY' found"
	        GP_FOLDER_LIST="$GP_FOLDER_LIST\n$APP_FOLDER_REPOSITORY/$APP_GROUP/$APP_PROJECT"
	    fi;
        
        
    done;

else

	for GP_FOLDER in $GP_FOLDER_LIST; do
		if [ -d "$GP" ]; then
			REPO=$(dirname $(dirname "$GP_FOLDER"))
			GROUP=$(basename $(dirname "$GP_FOLDER"))
			PROJECT=$(basename "$GP_FOLDER")

	        (($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' repository '$REPO' "
	        #GP_FOLDER_LIST="$GP_FOLDER_LIST\n$APP_FOLDER_ARCHIVES/$APP_GROUP/$APP_PROJECT"
	    fi;
	done;

fi;


# Repository found

#GP_FOLDER_LIST_UNIQ=$(echo -e $GP_FOLDER_LIST | sort -u)
GP_FOLDER_LIST_UNIQ=$(ls -d $GP_FOLDER_LIST 2>/dev/null | grep -v " " | sort -u)
#GP_FOLDER_LIST_UNIQ="${GP_FOLDER_LIST_UNIQ// /\\ }"
GP_FOLDER_LIST_UNIQ_COUNT=$(echo $GP_FOLDER_LIST_UNIQ | wc -w)

# echo "GP_FOLDER_LIST=$GP_FOLDER_LIST";
# echo "GP_FOLDER_LIST_UNIQ=$GP_FOLDER_LIST_UNIQ";
# echo -e $GP_FOLDER_LIST | sed "s/ /\\ /g"


#(($VERBOSE)) && for RF in echo -e $GP_FOLDER_LIST; do echo "#[INFOO]    "$RF; done



(($VERBOSE)) && echo "#[INFO] DEJAVU database repository/group/project found [$GP_FOLDER_LIST_UNIQ_COUNT]:"
#(($VERBOSE)) && echo $GP_FOLDER_LIST_UNIQ
(($VERBOSE)) && for RF in $GP_FOLDER_LIST_UNIQ; do echo "#[INFO]    "$RF; done



### DEJAVU database copy file
##############################

(($VERBOSE)) && echo "#"
#(($VERBOSE)) && echo "#[INFO] DEJAVU database file copy"
echo "#[INFO] DEJAVU database file copy"

(($DEBUG)) && echo "#[INFO] DEJAVU database file pattern '$VCF_PATTERN'"



for GP_FOLDER in $GP_FOLDER_LIST_UNIQ; do

	REPO=$(dirname $(dirname $GP_FOLDER))
	GROUP=$(basename $(dirname $GP_FOLDER))
	PROJECT=$(basename $GP_FOLDER)

#echo "$REPO/$GROUP/$PROJECT"
#continue

	# Filter
	SAMPLE_EXCLUDE_PARAM_GREP=""
	for SAMPLE_FILTER in $SAMPLE_EXCLUDE; do

		SAMPLE_FILTER_GROUP=$(echo $SAMPLE_FILTER | awk -F/ '{print $1}')
		SAMPLE_FILTER_PROJECT=$(echo $SAMPLE_FILTER | awk -F/ '{print $2}')
		SAMPLE_FILTER_RUN=$(echo $SAMPLE_FILTER | awk -F/ '{print $3}')
		SAMPLE_FILTER_SAMPLE=$(echo $SAMPLE_FILTER | awk -F/ '{print $4}')

		if [[ $GROUP =~ $SAMPLE_FILTER_GROUP ]] && [[ $PROJECT =~ $SAMPLE_FILTER_PROJECT ]]; then
			if [ "$SAMPLE_EXCLUDE_PARAM_GREP" == "" ]; then
				SEP=""
			else
				SEP="|"
			fi;
			SAMPLE_EXCLUDE_PARAM_GREP="$SAMPLE_EXCLUDE_PARAM_GREP$SEP$REPO/$GROUP/$PROJECT/$SAMPLE_FILTER_RUN/$SAMPLE_FILTER_SAMPLE"
		fi;

	done;

	# No filter
	if [ "$SAMPLE_EXCLUDE_PARAM_GREP" == "" ]; then
		SAMPLE_EXCLUDE_PARAM_GREP="ALLSAMPLEARESELECTED"
	fi;


	# NB VARIANT
	NB_VCF=$(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN 2>/dev/null | grep -vE $SAMPLE_EXCLUDE_PARAM_GREP 2>/dev/null | wc -l)

	(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' repository '$REPO' $NB_VCF VCF found"

	# If at least 1 vcf
	if [ $NB_VCF -gt 0 ]; then
		
		> $MK.$GROUP.$PROJECT.log
		> $MK.$GROUP.$PROJECT.err

		# TMP folder creation
		mkdir -p $TMP/$GROUP/$PROJECT
		#cp -f $(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN) $TMP/$GROUP/$PROJECT/ 2>/dev/null
		cp -f $(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN | grep -vE $SAMPLE_EXCLUDE_PARAM_GREP) $TMP/$GROUP/$PROJECT/ 2>/dev/null
		NB_VCF_FOUND=$(ls $TMP/$GROUP/$PROJECT/*.vcf.gz | wc -w)
		(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' $NB_VCF_FOUND VCF considered files (some files/samples may be found multiple times)"
		
	fi;

done;

### Database generation process

(($VERBOSE)) && echo "#"
#(($VERBOSE)) && echo "#[INFO] DEJAVU database generation process"
echo "#[INFO] DEJAVU database generation process"

for GP_FOLDER in $GP_FOLDER_LIST_UNIQ; do

	REPO=$(dirname $(dirname "$GP_FOLDER"))
	GROUP=$(basename $(dirname "$GP_FOLDER"))
	PROJECT=$(basename "$GP_FOLDER")

	NB_VCF=$(ls -l $TMP/$GROUP/$PROJECT/* 2>/dev/null | wc -l);
	#echo "NBVCF: $NB_VCF"; exit 0;

	#(($VERBOSE)) && echo "#"
	#(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' process"

	if (($NB_VCF)); then

		if [ ! -s "$DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done" ]; then

			# MK files
			> $MK
			echo "%.vcf.gz.tbi: %.vcf.gz
				$TABIX $<

			" > $MK

			VCFGZ_LIST=""
			#for VCF in $(ls $TMP/$GROUP/$PROJECT/*); do
			VCFGZ_NB=0
			for VCF in $TMP/$GROUP/$PROJECT/*; do
				SAMPLE_NAME=$(basename $VCF | cut -d. -f1)
				#echo "VCF: "$VCF
				if [ -s $VCF ] && (($(grep ^# -cv $VCF))); then
					echo "$VCF.simple.vcf.gz: $VCF
						mkdir $<.sort.
						if $BCFTOOLS annotate -x FILTER,QUAL,ID,INFO $< 1>/dev/null 2>/dev/null; then \
							$BCFTOOLS annotate -x FILTER,QUAL,ID,INFO $< | $BCFTOOLS norm -m -any -c s --fasta-ref $GENOMES/current/$ASSEMBLY.fa | $BCFTOOLS norm --rm-dup=exact | $BCFTOOLS +fixploidy  -- -f 2 | $BCFTOOLS +setGT  -- -t . -n 0 | $BCFTOOLS sort -T $<.sort. -o \$@ -O z 2>/dev/null; \
						else \
							echo '##fileformat=VCFv4.1' > \$@.tmp; \
							echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$SAMPLE_NAME' >> \$@.tmp; \
							$BGZIP -c \$@.tmp > \$@; \
							rm -f \$@.tmp; \
						fi;
						$TABIX \$@;
						rm -f $<.sort.*
					" >> $MK
					VCFGZ_LIST="$VCFGZ_LIST $VCF.simple.vcf.gz"
					((VCFGZ_NB++))
					#| sed s/ID=PL,Number=G/ID=PL,Number=./gi ,^FORMAT/GT
				fi;
			done

			#echo $VCFGZ_LIST > $TMP/$GROUP/$PROJECT/VCF_LIST
			
			# Minimum VCF
			echo "$TMP/$GROUP/$PROJECT/dejavu.simple.vcf: $VCFGZ_LIST" >> $MK
			if [ $VCFGZ_NB -gt 1 ]; then
				echo "	$BCFTOOLS merge --force-samples $TMP/$GROUP/$PROJECT/*.simple.vcf.gz | $BCFTOOLS norm -m -any | $BCFTOOLS norm --rm-dup=exact | $BCFTOOLS +setGT  -- -t . -n 0 | $BCFTOOLS +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS > \$@;" >> $MK
			else
				echo "	$BCFTOOLS norm -m -any $VCFGZ_LIST | $BCFTOOLS +setGT  -- -t . -n 0 | $BCFTOOLS +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS > \$@;" >> $MK
			fi;

			# BCFTOOLS stats
			echo "$TMP/$GROUP/$PROJECT/dejavu.stats.bcftools: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
				$BCFTOOLS stats $TMP/$GROUP/$PROJECT/dejavu.simple.vcf > $TMP/$GROUP/$PROJECT/dejavu.stats.bcftools
				$BCFTOOLS plugin counts $TMP/$GROUP/$PROJECT/dejavu.simple.vcf > $TMP/$GROUP/$PROJECT/dejavu.stats.bcftools.counts
			" >> $MK

			# VCFSTATS stats
			if [ "$VCFSTATS" != "" ]; then
				echo "$TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats: $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi
					mkdir -p $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats
					java -jar $VCFSTATS --inputFile $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz --outputDir $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats --referenceFile $GENOMES/current/$ASSEMBLY.fa \$\$(bgzip -dc $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz |  perl -ne 'print \"\$\$1\n\" if /##INFO=<ID=(.*?),/' | awk '{print \"--infoTag \"\$\$1\":All\"}')
				" >> $MK
			fi;

			# DejaVu
			echo "$TMP/$GROUP/$PROJECT/dejavu.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
				cp $TMP/$GROUP/$PROJECT/dejavu.simple.vcf $TMP/$GROUP/$PROJECT/dejavu.vcf
			" >> $MK

			# Annotation
			if [ "$DEJAVU_ANNOTATION" != "" ]; then
				echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
					+$HOWARD --input=$< --output=\$@.tmp.vcf --config=$HOWARD_CONFIG --config_annotation=$HOWARD_CONFIG_ANNOTATION --annotation=$DEJAVU_ANNOTATION --calculation=$DEJAVU_CALCULATION --nomen_fields=$DEJAVU_NOMEN_FIELDS --annovar_folder=$ANNOVAR --annovar_databases=$ANNOVAR_DATABASES --snpeff_jar=$SNPEFF --snpeff_databases=$SNPEFF_DATABASES --multithreading --threads=$THREADS --snpeff_threads=$THREADS --split=1000 --tmp=$TMP_FOLDER_TMP --env=$CONFIG_TOOLS
					mkdir $<.sort.
					$BCFTOOLS sort -T $<.sort. \$@.tmp.vcf | $BCFTOOLS annotate -x INFO/AN,INFO/AC,INFO/AF  | $BCFTOOLS +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS > \$@;
					#rm \$@.tmp.vcf
					rm \$@.tmp*
					rm -f $<.sort.*
				" >> $MK
				# --multithreading --threads=$THREADS
				# --norm_options='--multiallelics=-any,--rm-dup=exact'

				# echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf: $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz  $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi
				# 	if ((\$\$($BCFTOOLS view $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz -h | grep '##INFO=<ID=ANN,' -c))); then \
				# 		$BCFTOOLS view -O v $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz | sed -e 's/\([;=[:space:]]\)ANN\([,;=[:space:]]\)/\1EFF\2/' | bcftools view -O z -o \$@.tmp.eff.vcf.gz; \
				# 		$TABIX \$@.tmp.eff.vcf.gz; \
				# 		$BCFTOOLS annotate -a \$@.tmp.eff.vcf.gz -c EFF $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz -o \$@; \
				# 		rm \$@.tmp*; \
				# 	else \
				# 		$BCFTOOLS view $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz > $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf; \
				# 	fi;
				# " >> $MK

				# java -Xmx4G -jar /STARK/tools/snpeff/current/bin/snpEff.jar -i vcf -classic -formatEff  -o vcf hg19 annotated.vcf > annotated.eff4.vcf

			else
				echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
					cp $TMP/$GROUP/$PROJECT/dejavu.simple.vcf $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf
				" >> $MK
			fi;


			# snpEff
			if [ "$SNPEFF" != "" ] && [ -e $SNPEFF ]; then
				echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf: $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz  $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi
						$JAVA -Xmx4G -jar $SNPEFF -i vcf -classic -formatEff -o vcf $ASSEMBLY $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz -dataDir $SNPEFF_DATABASES > $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf
					" >> $MK
			fi;

			# echo "$TMP/$GROUP/$PROJECT/VCF_LIST: $VCFGZ_LIST
			# ls $^ > $TMP/$GROUP/$PROJECT/VCF_LIST

			# " >> $MK

			echo "%.vcf.gz: %.vcf
				$BGZIP -c $< > \$@

			" >> $MK

			echo "%.tsv.gz: %.tsv
				$BGZIP -c $< > \$@

			" >> $MK


			echo "%.gz.tbi: %.gz
				$TABIX $<

			" >> $MK

			echo "%.tsv: %.vcf
				+$HOWARD --input=$< --output=\$@ --env=$CONFIG_TOOLS --translation=TSV

			" >> $MK


			#echo "awk -F'\t' '{AF=\$\$6/($NB_VCF_FOUND*2)} {print \$\$1\"\\t\"\$\$2\"\\t\"\$\$3\"\\t\"\$\$4\"\\t\"AF}'"
			#exit 0

			#echo "$DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.txt: $TMP/$GROUP/$PROJECT/dejavu.vcf
			echo "$TMP/$GROUP/$PROJECT/dejavu.percent: $TMP/$GROUP/$PROJECT/dejavu.vcf
				$BCFTOOLS query -f'%CHROM\t%POS\t%REF\t%ALT\t%AF\t%NS\t%AN\t%AC\t%AC_Hom\t%AC_Het\n' $< > \$@

			" >> $MK
			#$BCFTOOLS query -f'%CHROM\t%POS\t%REF\t%ALT\t%AF\t%NS\t%AN\t%AC\t%AC_Hom\t%AF_Het\n' $< > \$@
			echo "$TMP/$GROUP/$PROJECT/dejavu.annovar: $TMP/$GROUP/$PROJECT/dejavu.vcf
				perl $ANNOVAR/convert2annovar.pl --format vcf4old  --allallele --outfile \$@.tmp $<
				cat \$@.tmp | cut -f1-5 > \$@
				rm -f \$@.tmp

			" >> $MK

			echo "$TMP/$GROUP/$PROJECT/dejavu.annovar.percent: $TMP/$GROUP/$PROJECT/dejavu.annovar $TMP/$GROUP/$PROJECT/dejavu.percent
				paste $^ | cut -f1-5,10-15 > \$@

			" >> $MK
			
			echo "$TMP/$GROUP/$PROJECT/dejavu.txt: $TMP/$GROUP/$PROJECT/dejavu.annovar.percent
				#perl $ANNOVAR/index_annovar.sh $< --outfile \$@
				perl $SCRIPT_DIR/index_annovar.pl $< --outfile \$@

			" >> $MK

			#cat $MK
			#(($VERBOSE)) && echo "#[INFO] DEJAVU database generation process..."
			#make -j $THREADS -f $MK $TMP/$GROUP/$PROJECT/dejavu.txt $TMP/$GROUP/$PROJECT/dejavu.vcf $TMP/$GROUP/$PROJECT/dejavu.tsv $TMP/$GROUP/$PROJECT/dejavu.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.tsv $TMP/$GROUP/$PROJECT/dejavu.stats.bcftools $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz.tbi 1>>$MK.$GROUP.$PROJECT.log 2>>$MK.$GROUP.$PROJECT.err
			make -j $THREADS -f $MK $TMP/$GROUP/$PROJECT/dejavu.txt $TMP/$GROUP/$PROJECT/dejavu.vcf $TMP/$GROUP/$PROJECT/dejavu.tsv.gz $TMP/$GROUP/$PROJECT/dejavu.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.tsv.gz $TMP/$GROUP/$PROJECT/dejavu.stats.bcftools $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz.tbi 1>>$MK.$GROUP.$PROJECT.log 2>>$MK.$GROUP.$PROJECT.err
			# grep "\*\*\*" $MK.err
			(($DEBUG)) && grep "\*\*\*" $MK.$GROUP.$PROJECT.log -B30
			(($DEBUG)) && grep "\*\*\*" $MK.$GROUP.$PROJECT.err -B30
			
			# $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi 

			# echo
			if (($(cat $MK.$GROUP.$PROJECT.log $MK.$GROUP.$PROJECT.err | grep "\*\*\*" -c))); then
				echo "#[ERROR] File '$DEJAVU/$RELEASE/dejavu.$GROUP.$PROJECT.txt' generation..."
				(($DEBUG)) && cat $MK.$GROUP.$PROJECT.log $MK.$GROUP.$PROJECT.err | grep "\*\*\*" -B 80
				#exit 1
			else
			
				# STATS
				> $TMP/$GROUP/$PROJECT/dejavu.stats.txt

				# SAMPLES
				NB_SAMPLES=$(echo $(grep "^#CHROM" $TMP/$GROUP/$PROJECT/dejavu.vcf | wc -w)" - 9" | bc)
				NB_VARIANTS=$(grep -cv "^#" $TMP/$GROUP/$PROJECT/dejavu.vcf)
				
				echo "### Number of samples: $NB_SAMPLES" >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
				echo "### Number of variants: $NB_VARIANTS" >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
				
				# VARIANTS
				echo "### variants frequency" >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
				# for p in $(seq 0 100); do
				#for p in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do awk -F"\t" -v p=$p '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print "# "SUM+0" variants out of "TOT" ("PERC"%) found more than "(p*100)"% in the set"}' $TMP/$GROUP/$PROJECT/dejavu.txt; done >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
				for p in $(seq 1 100); do awk -F"\t" -v p=$( echo "scale=2; $p/100" | bc) '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print "# "SUM+0" variants out of "TOT" ("PERC"%) found more than "(p*100)"% in the set"}' $TMP/$GROUP/$PROJECT/dejavu.txt; done >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
				#cat $TMP/$GROUP/$PROJECT/dejavu.stats.txt

				# TSV
				echo "$NB_SAMPLES" >> $TMP/$GROUP/$PROJECT/dejavu.stats.nb_samples
				echo "$NB_VARIANTS" >> $TMP/$GROUP/$PROJECT/dejavu.stats.nb_variants

				echo -e "#NB_variant\tTOT_variant\tPercent\tmore_than" > $TMP/$GROUP/$PROJECT/dejavu.stats.tsv
				#for p in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do awk -F"\t" -v p=$p '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print SUM+0"\t"TOT"\t"PERC"\t"(p*100)}' $TMP/$GROUP/$PROJECT/dejavu.txt; done > $TMP/$GROUP/$PROJECT/dejavu.stats.frequency.tsv
				for p in $(seq 1 100); do awk -F"\t" -v p=$( echo "scale=2; $p/100" | bc) '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print SUM+0"\t"TOT"\t"PERC"\t"(p*100)}' $TMP/$GROUP/$PROJECT/dejavu.txt; done >> $TMP/$GROUP/$PROJECT/dejavu.stats.tsv
				#cat $TMP/$GROUP/$PROJECT/dejavu.stats.tsv



				# find assembly prefix
				ASSEMBLY_PREFIX_VCF=$(grep "^##reference=" $TMP/$GROUP/$PROJECT/dejavu.vcf | cut -d= -f2 | xargs basename | sed 's/.fa$//' | sed 's/.fasta$//')
				if [ "$ASSEMBLY_PREFIX_VCF" == "" ]; then
					ASSEMBLY_PREFIX=$ASSEMBLY_PREFIX_DEFAULT
				else 
					ASSEMBLY_PREFIX=$ASSEMBLY_PREFIX_VCF
				fi;
				ASSEMBLY_PREFIX_ANNOVAR=$ASSEMBLY_PREFIX"_"
				
				# SUFFIX
				SUFFIX_ANNOVAR=$SUFFIX
				
				PATTERN=$ASSEMBLY_PREFIX_ANNOVAR"dejavu."$GROUP.$PROJECT$SUFFIX_ANNOVAR
				mkdir -p $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/

				cp $TMP/$GROUP/$PROJECT/dejavu.txt $DEJAVU_FOLDER_ANNOVAR/$PATTERN.txt
				cp $TMP/$GROUP/$PROJECT/dejavu.txt.idx $DEJAVU_FOLDER_ANNOVAR/$PATTERN.txt.idx

				#cp $TMP/$GROUP/$PROJECT/dejavu.vcf $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/minimal.vcf
				cp $TMP/$GROUP/$PROJECT/dejavu.tsv.gz $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/minimal.tsv.gz
				cp $TMP/$GROUP/$PROJECT/dejavu.vcf.gz $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/minimal.vcf.gz
				cp $TMP/$GROUP/$PROJECT/dejavu.vcf.gz.tbi $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/minimal.vcf.gz.tbi
				#cp $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.vcf
				cp $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.vcf.gz
				cp $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.vcf.gz.tbi
				cp $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.eff.vcf.gz
				cp $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz.tbi $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.eff.vcf.gz.tbi
				cp $TMP/$GROUP/$PROJECT/dejavu.annotated.tsv.gz $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/annotated.tsv.gz

				mkdir -p $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/$DEJAVU_SUBFOLDER_STATS
				cp -R $TMP/$GROUP/$PROJECT/dejavu.stats* $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/$DEJAVU_SUBFOLDER_STATS/

				(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' generated with $NB_SAMPLES samples and $NB_VARIANTS variants"
				

				# end
				#(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' generated for release $RELEASE"
				echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' generated for release $RELEASE" > $DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done
				echo "#[INFO] release=$RELEASE" >> $DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done
				echo "#[INFO] samples=$NB_SAMPLES" >> $DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done
				echo "#[INFO] variants=$NB_VARIANTS" >> $DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done
				(($VERBOSE)) && echo "$RELEASE" > $DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.release
				

				# STARK module json

				# Database definition
				if [ ! -s $DEJAVU/STARK.database ]; then
					echo '

						{
							"code": "dejavu",
							"name": "DejaVu",
							"fullname": "STARK DejaVu databases",
							"website": "",
							"description": "STARK DejaVu databases is a compilation of all samples variants for each group/project, useful to calculate population frequencies"
						}
				
					' > $DEJAVU/STARK.database
				fi;

				# CLeaning
				rm -rf $TMP/$GROUP/$PROJECT


			fi;
		
		fi;

	else

		(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' without VCF files"

	fi;

done


# STARK.database and STARK.database.release
echo '

		{
			"release": "'$RELEASE'",
			"date": "'$RELEASE'",
			"files": [ "release" ],
			"assembly": [ "'$ASSEMBLY'" ],
			"download": {
				"methode": "DejaVu Databases generation script ['$SCRIPT_RELEASE'-'$SCRIPT_DATE']",
				"date": "'.$(date).'"
			}
		}

' > $DEJAVU/$RELEASE/STARK.database.release


# Latest symlink
rm -f $DEJAVU/latest
ln -s $RELEASE/ $DEJAVU/latest



(($VERBOSE)) && echo "#"
(($VERBOSE)) && echo "#[INFO] DEJAVU database release '$RELEASE' done."
(($VERBOSE)) && echo "$RELEASE" > $DEJAVU/$RELEASE/release

rm -rf $TMP

exit 0;