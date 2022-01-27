#!/bin/bash
#################################
##
## NGS environment
##
#################################

SCRIPT_NAME="STARK_DEJAVU"
SCRIPT_DESCRIPTION="STARK DEJAVU ANNOVAR databases generation"
SCRIPT_RELEASE="0.12.3"
SCRIPT_DATE="27/01/2022"
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
RELEASE_NOTES=$RELEASE_NOTES"# 0.12.1-11/01/2021: Add group project filter, annotation, calculation and nomen fields parameters, rebuld search goup/project folders.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.12.2-29/10/2021: Fix error in input vcf.\n";
RELEASE_NOTES=$RELEASE_NOTES"# 0.12.3-27/01/2022: Fix error dejavu generation.\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Configuration
ENV_CONFIG=$(find -L $SCRIPT_DIR/.. -name config.app)

source $ENV_CONFIG 1>/dev/null 2>/dev/null

ENV_TOOLS=$(find -L $SCRIPT_DIR/.. -name tools.app)

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
	echo "# --app_folder|--application_folder=<FOLDER>    STARK Application folder";
	echo "#                                               Used to detect STARK Repository folders";
	echo "#                                               Default: default STARK Application folder";
	echo "# -r|--repo_folder|--repository_folder=<LIST>   STARK Repository folders";
	echo "#                                               List of STARK Repository folders containing group and project results";
	echo "#                                               Format: 'folder1[,folder2...]'";
	echo "#                                               Default: default STARK Repository folder";
	echo "# --dejavu_folder=<FOLDER>                      Output DejaVu folder";
	echo "#                                               Default: '.'";
	echo "# --dejavu_release=<STRING>                     Output DejaVu release";
	echo "#                                               Default: 'date +%Y%m%d-%H%M%S'";
	# echo "# --dejavu_release_latest                       Define DejaVu release as latest (symlink, previous removed)";
	# echo "#                                               Default: no";
	echo "# --dejavu_release_symlink=<LIST>               Create symlink to DejaVu release (symlink, previous removed)";
	echo "#                                               Format: 'symlink1[,symlink2]'";
	echo "#                                               Example: 'latest,current'";
	echo "#                                               Default: ''";
	echo "# --dejavu_previous_release=<STRING>            Previous DejaVu release";
	echo "#                                               In order to detect changes (new vcf files), usually 'latest' release ";
	echo "#                                               Default: '' (full generation of DejaVu)";
	echo "# --dejavu_previous_copy_mode=<STRING>          Previous DejaVu release copy mode";
	echo "#                                               Previous DejaVu will be copied (presistence) or linked (low disk space)";
	echo "#                                               Options: 'symlink', 'copy'";
	echo "#                                               Default: 'symlink'";
	echo "# --dejavu_annotation=<LIST>                    Output VCF DejaVu annotation";
	echo "#                                               Default: 'HOWARD_ANNOTATION_REPORT' STARK parameter";
	echo "#                                               Example: 'ALL,snpeff' for all annotations";
	echo "#                                               Tip: use 'none' for no annotation, 'ALL,snpeff,snpeff_hgvs' for ALL annotation";
	echo "# --dejavu_calculation=<LIST>                   Output VCF DejaVu calculation";
	echo "#                                               Default: 'HOWARD_CALCULATION_REPORT' STARK parameter";
	echo "#                                               Example: 'VAF,VAF_STATS,DP_STATS,VARTYPE,NOMEN'";
	echo "#                                               Tip: use 'none' for no calculation";
	echo "# --dejavu_nomen_fields=<LIST>                  Output VCF DejaVu NOMEN field";
	echo "#                                               Default: 'HOWARD_NOMEN_FIELDS' STARK parameter";
	echo "#                                               Example: 'hgvs', 'snpeff_hgvs'";
	echo "# --dejavu_vcfstats=<STRING>                    Output VCFStats";
	echo "#                                               Default: no VCFStats output";
	echo "#                                               Example: 'hgvs', 'snpeff_hgvs'";
	echo "# --sample_exclude=<LIST>                       Exclude sample pattern (regexp)";
	echo "#                                               Format: '<group>/<project>/<sample_pattern>[,<group>/<project>/<run>/<sample_pattern>]'";
	echo "#                                               Example: 'GENOME/GERMLINE/.*/.*CORIEL.*' to exclude all *CORIEL* samples";
	echo "#                                               Default: ''";
	echo "# --sample_exclude_file=<FILE>                  Exclude sample pattern (regexp) within a file";
	echo "#                                               Format: same as --sample_exclude parameter";
	echo "#                                               Default: <STARK_FOLDER_CONFIG>/dejavu/sample_exclude.conf'";
	echo "# --group_project_list=<LIST>                   Include only group/project list (shell like)";
	echo "#                                               These folders must be well structured as group/project within repository folders";
	echo "#                                               Format: '<group>/<project>[,<group>/<project>]'";
	echo "#                                               Example: 'GENOME/*' to include only run from group GENOME and all project";
	echo "#                                               Default: '' (empty), filter will be used";
	echo "# --group_project_list_file=<FILE>              Include only group/project/run list (shell like) within a file";
	echo "#                                               Format: same as --group_project_list parameter";
	echo "#                                               Default: <STARK_FOLDER_CONFIG>/dejavu/group_project_list.conf'";
	echo "# --group_project_filter=<LIST>                 Include only group/project/run pattern (shell like)";
	echo "#                                               These patterns will be used to detect well structured folders, as group/project, within repository folders";
	echo "#                                               Format: '<group>/<project>/<run>[,<group>/<project>/<run>]'";
	echo "#                                               Example: 'GENOME/*/19*' to include only run of year 2019 from group GENOME and all project";
	echo "#                                               Default: '*/*/*', all groups, projects and runs";
	echo "# --group_project_filter_file=<FILE>            Include only group/project/run pattern (shell like) within a file";
	echo "#                                               Format: same as --group_project_filter parameter";
	echo "#                                               Default: <STARK_FOLDER_CONFIG>/dejavu/group_project_filter.conf'";
	echo "# --tmp=<FOLDER>                                Temporary folder";
	echo "#                                               Default: default STARK Temporary folder";
	echo "# --bcftools=<FILE>                             BCFTools application binary";
	echo "#                                               Default: default STARK configuration or 'bcftools'";
	echo "# --tabix=<FILE>                                TABix application binary";
	echo "#                                               Default: default STARK configuration or 'tabix'";
	echo "# --bgzip=<FILE>                                BGZip application binary";
	echo "#                                               Default: default STARK configuration or 'bgzip'";
	echo "# --annovar=<FILE>                              ANNOVAR application binary folder";
	echo "#                                               Default: default STARK configuration or ''";
	echo "# --vcfstats=<FILE>                             VCFStats application jar";
	echo "#                                               Default: default STARK configuration or '' or detected with wheris command";
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
ARGS=$(getopt -o "e:r:vdnh" --long "env:,app:,application:,app_folder:,application_folder:,repo_folder:,repository_folder:,dejavu_folder:,dejavu_release:,dejavu_release_latest,dejavu_release_symlink:,dejavu_previous_release:,dejavu_previous_copy_mode:,dejavu_annotation:,dejavu_calculation:,dejavu_nomen_fields:,dejavu_vcfstats,sample_exclude:,sample_exclude_file:,group_project_list:,group_project_list_file:,group_project_filter:,group_project_filter_file:,tmp:,bcftools:,tabix:,bgzip:,annovar:,vcfstats:,verbose,debug,release,help" -- "$@" 2> /dev/null)

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
		--dejavu_release_latest)
			RELEASE_LATEST=1
			shift 1
			;;
		--dejavu_release_symlink)
			RELEASE_SYMLINK=$(echo "$2" | tr "," " ");
			shift 2
			;;
		--dejavu_previous_release)
			PREVIOUS_RELEASE="$2";
			shift 2
			;;
		--dejavu_previous_copy_mode)
			PREVIOUS_COPY_MODE="$2";
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
		--dejavu_vcfstats)
			DEJAVU_VCFSTATS=1
			shift 1
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
		--vcfstats)
			VCFSTATS="$2";
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
REPO_FOLDER=$(echo $REPO_FOLDER | tr "," " " | tr " " "\n" | sort -u)


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


if [ -z "$PREVIOUS_COPY_MODE" ]; then
	PREVIOUS_COPY_MODE="symlink"
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
#(($VERBOSE)) && [ "$GROUP_PROJECT_FILTER" != "*/*/*" ] && for GPF in $(echo "$GROUP_PROJECT_FILTER" | tr " " "\n"); do echo "#[INFO]    "$GPF; done
#(($VERBOSE)) && for GPF in $GROUP_PROJECT_FILTER; do echo "#[INFO]    $GPF"; done
#($VERBOSE)) && for GPF in "$GROUP_PROJECT_FILTER"; do echo "$GPF" | sed "s/[^ ]/#[INFO]    $GPF/gi"; done
(($VERBOSE)) && echo "#[INFO]    $GROUP_PROJECT_FILTER"

(($VERBOSE)) && echo "#[INFO] SAMPLE EXCLUDE FILE:"
(($VERBOSE)) && echo "#[INFO]    $SAMPLE_EXCLUDE_FILE"
(($VERBOSE)) && echo "#[INFO] SAMPLE EXCLUDE:"
(($VERBOSE)) && for SE in $SAMPLE_EXCLUDE; do echo "#[INFO]    "$SE; done

(($VERBOSE)) && echo "#[INFO] DEJAVU ANNOTATION: $DEJAVU_ANNOTATION"
(($VERBOSE)) && echo "#[INFO] DEJAVU CALCULATION: $DEJAVU_CALCULATION"
(($VERBOSE)) && echo "#[INFO] DEJAVU NOMEN FIELDS: $DEJAVU_NOMEN_FIELDS"
(($VERBOSE)) && echo "#[INFO] DEJAVU VCFSTATS: "$( (($DEJAVU_VCFSTATS)) && echo "Yes" || echo "No")

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
				GP_FOLDER_LIST="$GP_FOLDER_LIST	$RF/$GPR_LIST/"
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
GP_FOLDER_LIST_UNIQ=$(ls -d $(echo -e $GP_FOLDER_LIST) 2>/dev/null | grep -v " " | sort -u)
#GP_FOLDER_LIST_UNIQ="${GP_FOLDER_LIST_UNIQ// /\\ }"
GP_FOLDER_LIST_UNIQ_COUNT=$(echo $GP_FOLDER_LIST_UNIQ | wc -w)


(($VERBOSE)) && echo "#[INFO] DEJAVU database repository/group/project found [$GP_FOLDER_LIST_UNIQ_COUNT]:"
#(($VERBOSE)) && echo $GP_FOLDER_LIST_UNIQ
(($VERBOSE)) && for RF in $GP_FOLDER_LIST_UNIQ; do echo "#[INFO]    "$RF; done




### DEJAVU database copy file
##############################

#(($VERBOSE)) && echo "#"
#(($VERBOSE)) && echo "#[INFO] DEJAVU database file copy"
#echo "#[INFO] DEJAVU database file copy"




if false; then

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
	#NB_VCF=$(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN 2>/dev/null | grep -vE $SAMPLE_EXCLUDE_PARAM_GREP 2>/dev/null | wc -l)
	VCF_LIST=$(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN 2>/dev/null | grep -vE $SAMPLE_EXCLUDE_PARAM_GREP 2>/dev/null)
	NB_VCF=$(echo $VCF_LIST | wc -w)

	(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' repository '$REPO' $NB_VCF VCF found"

	# If at least 1 vcf
	if [ $NB_VCF -gt 0 ]; then
		
		> $MK.$GROUP.$PROJECT.log
		> $MK.$GROUP.$PROJECT.err

		# TMP folder creation
		mkdir -p $TMP/$GROUP/$PROJECT
		#cp -f $(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN) $TMP/$GROUP/$PROJECT/ 2>/dev/null
		#cp -f $(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN | grep -vE $SAMPLE_EXCLUDE_PARAM_GREP) $TMP/$GROUP/$PROJECT/ 2>/dev/null
		cp -f $VCF_LIST $TMP/$GROUP/$PROJECT/ 2>/dev/null
		NB_VCF_FOUND=$(ls $TMP/$GROUP/$PROJECT/*.vcf.gz 2>/dev/null | wc -w)
		(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' $NB_VCF_FOUND VCF considered files (some files/samples may be found multiple times)"
		
	fi;

done;

fi;

### Database generation process

#(($VERBOSE)) && echo "#"
#(($VERBOSE)) && echo "#[INFO] DEJAVU database generation process"
#echo "#[INFO] DEJAVU database generation process"

#echo "GP_FOLDER_LIST_UNIQ=$GP_FOLDER_LIST_UNIQ"


GP_LIST=""
for GP_FOLDER in $GP_FOLDER_LIST_UNIQ; do
	#echo $GP_FOLDER
	GROUP=$(basename $(dirname "$GP_FOLDER"))
	PROJECT=$(basename "$GP_FOLDER")
	#echo "GROUP=$GROUP PROJECT=$PROJECT"
	GP_LIST="$GP_LIST $GROUP/$PROJECT"
done;

GP_LIST_UNIQ=$(echo $GP_LIST | tr " " "\n" | sort -u)
#GP_FOLDER_LIST_UNIQ="${GP_FOLDER_LIST_UNIQ// /\\ }"
GP_LIST_UNIQ_COUNT=$(echo $GP_LIST_UNIQ | wc -w)


(($VERBOSE)) && echo "#"
(($VERBOSE)) && echo "#[INFO] DEJAVU database group/project found [$GP_LIST_UNIQ_COUNT]:"
#(($VERBOSE)) && echo $GP_FOLDER_LIST_UNIQ
(($VERBOSE)) && for GP in $GP_LIST_UNIQ; do echo "#[INFO]    "$GP; done


(($VERBOSE)) && echo "#"
(($VERBOSE)) && echo "#[INFO] DEJAVU database VCF file pattern '$VCF_PATTERN'"


# (($VERBOSE)) && echo "#"
# (($VERBOSE)) && echo "#[INFO] DEJAVU database 'GROUP/PROJECT' process"


#for GP_FOLDER in $GP_FOLDER_LIST_UNIQ; do
#for GP_FOLDER in $(find -L $TMP -mindepth 2 -maxdepth 2 -type d); do
for GP_FOLDER in $GP_LIST_UNIQ; do


	(($VERBOSE)) && echo "#"

	#REPO=$(dirname $(dirname "$GP_FOLDER"))
	GROUP=$(basename $(dirname "$GP_FOLDER"))
	PROJECT=$(basename "$GP_FOLDER")

	(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' process..."

	#(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' file copy..."


	for GP_FOLDER in $GP_FOLDER_LIST_UNIQ; do

		REPO=$(dirname $(dirname $GP_FOLDER))
		GROUP_FOLDER=$(basename $(dirname $GP_FOLDER))
		PROJECT_FOLDER=$(basename $GP_FOLDER)

		if [ "$GROUP" == "$GROUP_FOLDER" ] && [ "$PROJECT" == "$PROJECT_FOLDER" ]; then

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
			#NB_VCF=$(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN 2>/dev/null | grep -vE $SAMPLE_EXCLUDE_PARAM_GREP 2>/dev/null | wc -l)
			VCF_LIST=$(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN 2>/dev/null | grep -vE $SAMPLE_EXCLUDE_PARAM_GREP 2>/dev/null)
			NB_VCF=$(echo $VCF_LIST | wc -w)

			(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' $NB_VCF VCF files found in repository '$REPO'"

			# If at least 1 vcf
			if [ $NB_VCF -gt 0 ]; then
				
				> $MK.$GROUP.$PROJECT.log
				> $MK.$GROUP.$PROJECT.err

				# TMP folder creation
				mkdir -p $TMP/$GROUP/$PROJECT
				#cp -f $(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN) $TMP/$GROUP/$PROJECT/ 2>/dev/null
				#cp -f $(find -L $GP_FOLDER/*/*/ -maxdepth 1 -name '*'$VCF_PATTERN -a ! -name '*.*-*'$VCF_PATTERN | grep -vE $SAMPLE_EXCLUDE_PARAM_GREP) $TMP/$GROUP/$PROJECT/ 2>/dev/null
				#cp -f $VCF_LIST $TMP/$GROUP/$PROJECT/ 2>/dev/null
				cp -sup $VCF_LIST $TMP/$GROUP/$PROJECT/ 2>/dev/null
				
			fi;

		fi;

	done;



	NB_VCF_FOUND=$(ls $TMP/$GROUP/$PROJECT/*.vcf.gz 2>/dev/null | wc -w)
	(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' $NB_VCF_FOUND VCF files copied (some files/samples may be found multiple times)"


	#(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' generation process..."

	#NB_VCF=$(ls -l $TMP/$GROUP/$PROJECT/* 2>/dev/null | wc -l);
	NB_VCF=$NB_VCF_FOUND;
	#echo "NBVCF: $NB_VCF"; exit 0;

	#(($VERBOSE)) && echo "#"
	#(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' process"

	if (($NB_VCF)); then

		# DEV

		# CHECKSUM
		CHECKSUM_FILE=$TMP/$GROUP/$PROJECT/dejavu.checksum
		stat -c '%n %y %s' $TMP/$GROUP/$PROJECT/*.vcf.gz | while read line; do
			basename "$line"
		done | sha1sum | cut -d\  -f1 > $CHECKSUM_FILE
		CHECKSUM=$(cat $CHECKSUM_FILE)

		# stat -c '%n %y %s' $TMP/$GROUP/$PROJECT/*.vcf.gz | while read line; do
		# 	basename "$line"
		# done 
		
		CHECKSUM_DONE=0

		if [ "$PREVIOUS_RELEASE" != "" ] && [ -e $DEJAVU/$PREVIOUS_RELEASE/vcf/$GROUP/$PROJECT/checksum ]; then

			PREVIOUS_CHECKSUM_FILE=$DEJAVU/$PREVIOUS_RELEASE/vcf/$GROUP/$PROJECT/checksum
			CHECKSUM_PREVIOUS=$(cat $PREVIOUS_CHECKSUM_FILE)

			#stat -c "%y %s %n" $TMP/$GROUP/$PROJECT/*.vcf.gz 
			#stat -c "%y %s %n" $TMP/$GROUP/$PROJECT/*.vcf.gz | sha1sum

			#stat -c "%y %s %n" $TMP/$GROUP/$PROJECT/*.vcf.gz > $TMP/$GROUP/$PROJECT/checksum
			#echo $TMP/$GROUP/$PROJECT/checksum
			#cat $TMP/$GROUP/$PROJECT/checksum | sha1sum

			#CHECKSUM=$(stat -c "%y %s %n" $TMP/$GROUP/$PROJECT/*.vcf.gz | sha1sum)
			(($DEBUG)) && echo "[INFO] PREVIOUS_RELEASE=$PREVIOUS_RELEASE"
			(($DEBUG)) && echo "[INFO] PREVIOUS_COPY_MODE=$PREVIOUS_COPY_MODE"
			(($DEBUG)) && echo "[INFO] CHECKSUM=$CHECKSUM"
			(($DEBUG)) && echo "[INFO] CHECKSUM_PREVIOUS=$CHECKSUM_PREVIOUS"
			(($DEBUG)) && echo "[INFO] PREVIOUS_COPY_MODE=$PREVIOUS_COPY_MODE"

			if [ "$CHECKSUM" == "$CHECKSUM_PREVIOUS" ]; then
				(($DEBUG)) && echo "[INFO] same checksum"

				SYMLINK_ERR=0
				if [ "$PREVIOUS_COPY_MODE" == "symlink" ]; then
					mkdir -p $DEJAVU/$RELEASE/vcf/$GROUP
					if ! ln -s ../../../$(realpath $DEJAVU/$PREVIOUS_RELEASE/vcf/$GROUP/$PROJECT | rev | awk -F/ '{print $4}' | rev )/vcf/$GROUP/$PROJECT/ $DEJAVU/$RELEASE/vcf/$GROUP/$PROJECT; then
						SYMLINK_ERR=1
						PREVIOUS_COPY_MODE="copy"
					fi;
					
				fi;

				if [ "$PREVIOUS_COPY_MODE" == "copy" ]; then
					mkdir -p $DEJAVU/$RELEASE/vcf/$GROUP/$PROJECT
					cp -rp $DEJAVU/$PREVIOUS_RELEASE/vcf/$GROUP/$PROJECT/* $DEJAVU/$RELEASE/vcf/$GROUP/$PROJECT/
				fi;

				CHECKSUM_DONE=1
				
				(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' from previous release '$PREVIOUS_RELEASE' ('$PREVIOUS_COPY_MODE' copy mode)"

			fi;


			#RELEASE_PREVIOUS=latest
			# Softlink
			#ln -s ../../../$(realpath $DEJAVU/$RELEASE_PREVIOUS/vcf/$GROUP/$PROJECT | rev | awk -F/ '{print $4}' | rev )/vcf/$GROUP/$PROJECT $DEJAVU/$RELEASE/vcf/$GROUP/$PROJECT
			# Copy
			#cp -r $DEJAVU/$RELEASE_PREVIOUS/vcf/$GROUP/$PROJECT/* $DEJAVU/$RELEASE/vcf/$GROUP/$PROJECT/

		fi;


		#continue;




		if [ ! -s "$DEJAVU_FOLDER_LOG/dejavu.$GROUP.$PROJECT.done" ] && ! (($CHECKSUM_DONE)); then

			# MK files
			> $MK

			# TABIX
			echo "%.vcf.gz.tbi: %.vcf.gz
				$TABIX $<

			" > $MK

			# EMPTY
			echo "%.empty.vcf:
				# Header
				echo '##fileformat=VCFv4.1' > \$@;
				# Head first line
				echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO' >> \$@;

			" > $MK

			# VCF fix and simplify
			VCFGZ_LIST=""
			#for VCF in $(ls $TMP/$GROUP/$PROJECT/*); do
			VCFGZ_NB=0
			for VCF in $TMP/$GROUP/$PROJECT/*.vcf.gz; do
				SAMPLE_NAME=$(basename $VCF | cut -d. -f1)
				#echo "VCF: "$VCF
				if [ -s $VCF ] && (($(grep ^# -cv $VCF))); then
					echo "$VCF.simple.vcf.gz: $VCF $VCF.empty.vcf
						mkdir $<.sort.
						# if $JAVA -jar $PICARD FixVcfHeader -I $< -O $<.tmp.fixed.vcf; then \
						# 	echo '#[INFO] VCF well-formed for $< (FixVcfHeader)' ; \
						# else \
						# 	echo '#[ERROR] VCF not well-formed for $< (FixVcfHeader)' ; \
						# 	cp $VCF.empty.vcf $<.tmp.fixed.vcf; \
						# fi;
						ln -s $< $<.tmp.fixed0.vcf.gz;
						if zcat $<.tmp.fixed0.vcf.gz | sed 's/[^\x00-\x7F]//gi' | grep -v '^##Prioritize list is' | sed s/Number=R/Number=./g | sed s/Number=G/Number=./g | $BCFTOOLS sort -T $<.sort2. > $<.tmp.fixed2.vcf; then \
							echo '#[INFO] VCF well-formed for $VCF (sedBCFToolsSort)' ; \
						else \
							echo '#[ERROR] VCF not well-formed for $VCF (sedBCFToolsSort)' ; \
							cp $VCF.empty.vcf $<.tmp.fixed2.vcf; \
						fi;
						# If not correctly fixed for merge
						$BGZIP -c $<.tmp.fixed2.vcf -l 0 > $<.tmp.fixed2.vcf.gz;
						$TABIX $<.tmp.fixed2.vcf.gz;
						if ! $BCFTOOLS merge $<.tmp.fixed2.vcf.gz $<.tmp.fixed2.vcf.gz --force-samples 2>/dev/null; then \
							echo '#[ERROR] VCF not well-formed for $VCF (merge test)' ; \
							cp $VCF.empty.vcf $<.tmp.fixed2.vcf; \
						fi;
						$BGZIP -c $<.tmp.fixed2.vcf -l 0 > $<.tmp.fixed.vcf.gz;
						$TABIX $<.tmp.fixed.vcf.gz
						if $BCFTOOLS annotate -x FILTER,QUAL,ID,INFO $<.tmp.fixed.vcf.gz 1>/dev/null 2>/dev/null; then \
							$BCFTOOLS annotate -x FILTER,QUAL,ID,INFO $<.tmp.fixed.vcf.gz | $BCFTOOLS norm -m -any -c s --fasta-ref $GENOMES/current/$ASSEMBLY.fa | $BCFTOOLS norm --rm-dup=exact | $BCFTOOLS +fixploidy  -- -f 2 | $BCFTOOLS +setGT  -- -t . -n 0 | $BCFTOOLS sort -T $<.sort. -o \$@ -O z 2>/dev/null; \
						else \
							cp $VCF.empty.vcf \$@.tmp; \
							$BGZIP -c \$@.tmp -l 0 > \$@; \
							rm -rf \$@.tmp; \
						fi;
						$TABIX \$@;
						rm -rf $<.sort.* $<.tmp*
					" >> $MK
					VCFGZ_LIST="$VCFGZ_LIST $VCF.simple.vcf.gz"
					((VCFGZ_NB++))
					#| sed s/ID=PL,Number=G/ID=PL,Number=./gi ,^FORMAT/GT
				fi;
			done

							#echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$SAMPLE_NAME' >> $<.tmp.fixed2.vcf; \
							#echo '##fileformat=VCFv4.1' > $<.tmp.fixed2.vcf; \
							#echo '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO' >> $<.tmp.fixed2.vcf; \

			#echo $VCFGZ_LIST > $TMP/$GROUP/$PROJECT/VCF_LIST
			
			# Minimum VCF
			echo "$TMP/$GROUP/$PROJECT/dejavu.simple.vcf: $VCFGZ_LIST $TMP/$GROUP/$PROJECT/dejavu.simple.empty.vcf" >> $MK
			if [ $VCFGZ_NB -gt 1 ]; then
				echo "	if ! $BCFTOOLS merge --force-samples $TMP/$GROUP/$PROJECT/*.simple.vcf.gz | $BCFTOOLS norm -m -any -c s --fasta-ref $GENOMES/current/$ASSEMBLY.fa | $BCFTOOLS norm --rm-dup=exact | $BCFTOOLS +setGT  -- -t . -n 0 | $BCFTOOLS +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS > \$@; then \
							cp $TMP/$GROUP/$PROJECT/dejavu.simple.empty.vcf \$@; \
							echo '#[ERROR] VCF not well-formed for \$@' ; \
						fi;
				" >> $MK
			else
				echo "	$BCFTOOLS norm -m -any $VCFGZ_LIST | $BCFTOOLS +setGT  -- -t . -n 0 | $BCFTOOLS +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS > \$@;" >> $MK
			fi;

			# BCFTOOLS stats
			echo "$TMP/$GROUP/$PROJECT/dejavu.stats.bcftools: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
				$BCFTOOLS stats $TMP/$GROUP/$PROJECT/dejavu.simple.vcf > $TMP/$GROUP/$PROJECT/dejavu.stats.bcftools
				$BCFTOOLS plugin counts $TMP/$GROUP/$PROJECT/dejavu.simple.vcf > $TMP/$GROUP/$PROJECT/dejavu.stats.bcftools.counts
			" >> $MK

			# VCFSTATS stats
			if (($DEJAVU_VCFSTATS )) && [ "$VCFSTATS" != "" ]; then
				echo "$TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats.tar.gz: $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi
					mkdir -p $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats
					java -jar $VCFSTATS --inputFile $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz --outputDir $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats --referenceFile $GENOMES/current/$ASSEMBLY.fa -t $THREADS \$\$($BGZIP -dc $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz |  perl -ne 'print \"\$\$1\n\" if /##INFO=<ID=(.*?),/' | awk '{print \"--infoTag \"\$\$1\":All\"}')
					tar -zvcf $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats.tar.gz $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats
					rm -rf $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats
				" >> $MK
			else
				echo "$TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats.tar.gz: $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi
					mkdir -p $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats
					tar -zvcf $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats.tar.gz $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats
					rm -rf $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats
				" >> $MK
			fi;

			# DejaVu
			echo "$TMP/$GROUP/$PROJECT/dejavu.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
				cp $TMP/$GROUP/$PROJECT/dejavu.simple.vcf $TMP/$GROUP/$PROJECT/dejavu.vcf
			" >> $MK

			# Annotation
			HOWARD_split=$(echo " ( 10000 / sqrt($NB_VCF+1) ) + 1 " | bc)	# Prevent high memory 
			if [ "$DEJAVU_ANNOTATION" != "" ]; then
				echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.howard.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.simple.vcf.gz.tbi
					$BCFTOOLS view $< | cut -f1-8 > $<.reduced.vcf
					+$HOWARD --input=$<.reduced.vcf --output=\$@.tmp.tmp.vcf --config=$HOWARD_CONFIG --config_annotation=$HOWARD_CONFIG_ANNOTATION --annotation=$DEJAVU_ANNOTATION --calculation=$DEJAVU_CALCULATION --nomen_fields=$DEJAVU_NOMEN_FIELDS --annovar_folder=$ANNOVAR --annovar_databases=$ANNOVAR_DATABASES --snpeff_jar=$SNPEFF --snpeff_databases=$SNPEFF_DATABASES --multithreading --threads=$THREADS --snpeff_threads=$THREADS --tmp=$TMP_FOLDER_TMP --env=$CONFIG_TOOLS
					$BCFTOOLS sort -T $<.sort.1. \$@.tmp.tmp.vcf > \$@.tmp.tmp.tmp.vcf
					$BGZIP -c \$@.tmp.tmp.tmp.vcf -l 0 > \$@.tmp.tmp.tmp.vcf.gz
					$TABIX \$@.tmp.tmp.tmp.vcf.gz
					$BCFTOOLS merge $< \$@.tmp.tmp.tmp.vcf.gz | $BCFTOOLS sort -T $<.sort.2. > \$@.tmp.ann.vcf
					+$HOWARD --input=\$@.tmp.ann.vcf --output=\$@.tmp.vcf --config=$HOWARD_CONFIG --config_annotation=$HOWARD_CONFIG_ANNOTATION --calculation=$DEJAVU_CALCULATION --nomen_fields=$DEJAVU_NOMEN_FIELDS --annovar_folder=$ANNOVAR --annovar_databases=$ANNOVAR_DATABASES --snpeff_jar=$SNPEFF --snpeff_databases=$SNPEFF_DATABASES --multithreading --threads=$THREADS --snpeff_threads=$THREADS --split=$HOWARD_split --tmp=$TMP_FOLDER_TMP --env=$CONFIG_TOOLS
					#
					#+$HOWARD --input=$< --output=\$@.tmp.vcf --config=$HOWARD_CONFIG --config_annotation=$HOWARD_CONFIG_ANNOTATION --annotation=$DEJAVU_ANNOTATION --calculation=$DEJAVU_CALCULATION --nomen_fields=$DEJAVU_NOMEN_FIELDS --annovar_folder=$ANNOVAR --annovar_databases=$ANNOVAR_DATABASES --snpeff_jar=$SNPEFF --snpeff_databases=$SNPEFF_DATABASES --multithreading --threads=$THREADS --snpeff_threads=$THREADS --split=$HOWARD_split --tmp=$TMP_FOLDER_TMP --env=$CONFIG_TOOLS
					mkdir $<.sort.
					$BCFTOOLS sort -T $<.sort. \$@.tmp.vcf | $BCFTOOLS annotate -x INFO/AN,INFO/AC,INFO/AF  | $BCFTOOLS +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Hom,AC_Het,ExcHet,HWE,MAF,NS > \$@;
					#rm \$@.tmp.vcf
					rm -rf \$@.tmp*
					rm -rf $<.sort.*
				" >> $MK
				
			else
				echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.howard.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
					cp $TMP/$GROUP/$PROJECT/dejavu.simple.vcf $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf
				" >> $MK
			fi;


			# snpEff classic
			if [ "$SNPEFF" != "" ] && [ -e $SNPEFF ]; then
				echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf.gz  $TMP/$GROUP/$PROJECT/dejavu.simple.vcf.gz.tbi
						$JAVA -Xmx4G -jar $SNPEFF -i vcf -classic -formatEff -o vcf $ASSEMBLY $TMP/$GROUP/$PROJECT/dejavu.simple.vcf.gz -stats $TMP/$GROUP/$PROJECT/dejavu.stats.snpeff.html -dataDir $SNPEFF_DATABASES > $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf
					" >> $MK
			else
				echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf: $TMP/$GROUP/$PROJECT/dejavu.simple.vcf
						$BCFTOOLS view $TMP/$GROUP/$PROJECT/dejavu.simple.vcf.gz > $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf
					" >> $MK
			fi;

			# Annotated
			# TODO: check if EFF exists
			echo "$TMP/$GROUP/$PROJECT/dejavu.annotated.vcf: $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.howard.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.howard.vcf.gz.tbi
				if ! $BCFTOOLS annotate -a $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz -c EFF $TMP/$GROUP/$PROJECT/dejavu.annotated.howard.vcf.gz --threads $THREADS --single-overlaps -k > $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf; then \
					$BGZIP -dc $TMP/$GROUP/$PROJECT/dejavu.annotated.howard.vcf.gz > $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf ; \
				fi;
			" >> $MK



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
				rm -rf \$@.tmp

			" >> $MK

			echo "$TMP/$GROUP/$PROJECT/dejavu.annovar.percent: $TMP/$GROUP/$PROJECT/dejavu.annovar $TMP/$GROUP/$PROJECT/dejavu.percent
				paste $^ | cut -f1-5,10-15 > \$@

			" >> $MK
			
			echo "$TMP/$GROUP/$PROJECT/dejavu.txt: $TMP/$GROUP/$PROJECT/dejavu.annovar.percent
				#perl $ANNOVAR/index_annovar.sh $< --outfile \$@
				perl $SCRIPT_DIR/index_annovar.pl $< --outfile \$@

			" >> $MK

			#cat $MK
			cp -p $MK $MK.$GROUP.$PROJECT.mk
			#(($VERBOSE)) && echo "#[INFO] DEJAVU database generation process..."
			#make -j $THREADS -f $MK $TMP/$GROUP/$PROJECT/dejavu.txt $TMP/$GROUP/$PROJECT/dejavu.vcf $TMP/$GROUP/$PROJECT/dejavu.tsv $TMP/$GROUP/$PROJECT/dejavu.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.tsv $TMP/$GROUP/$PROJECT/dejavu.stats.bcftools $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz.tbi 1>>$MK.$GROUP.$PROJECT.log 2>>$MK.$GROUP.$PROJECT.err
			make -j $THREADS -f $MK $TMP/$GROUP/$PROJECT/dejavu.txt $TMP/$GROUP/$PROJECT/dejavu.vcf $TMP/$GROUP/$PROJECT/dejavu.tsv.gz $TMP/$GROUP/$PROJECT/dejavu.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi $TMP/$GROUP/$PROJECT/dejavu.annotated.tsv.gz $TMP/$GROUP/$PROJECT/dejavu.stats.bcftools $TMP/$GROUP/$PROJECT/dejavu.stats.vcfstats.tar.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.eff.vcf.gz.tbi 1>>$MK.$GROUP.$PROJECT.log 2>>$MK.$GROUP.$PROJECT.err
			# grep "\*\*\*" $MK.err
			(($DEBUG)) && grep "\*\*\*" $MK.$GROUP.$PROJECT.log -B30
			(($DEBUG)) && grep "\*\*\*" $MK.$GROUP.$PROJECT.err -B30
			
			# $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz.tbi 

			# echo
			if (($(cat $MK.$GROUP.$PROJECT.log $MK.$GROUP.$PROJECT.err | grep "\*\*\*" -c))) || (($(cat $MK.$GROUP.$PROJECT.log $MK.$GROUP.$PROJECT.err | grep "^\[E::" -c))); then
				echo "#[ERROR] File '$DEJAVU/$RELEASE/dejavu.$GROUP.$PROJECT.txt' generation..."
				(($DEBUG)) && cat $MK.$GROUP.$PROJECT.log $MK.$GROUP.$PROJECT.err | grep "\*\*\*" -B 80
				(($DEBUG)) && cat $MK.$GROUP.$PROJECT.log $MK.$GROUP.$PROJECT.err | grep "^\[E::" -B 80
				#exit 1
			else
			
				# STATS
				> $TMP/$GROUP/$PROJECT/dejavu.stats.txt

				# SAMPLES
				NB_SAMPLES=$(echo $(grep "^#CHROM" $TMP/$GROUP/$PROJECT/dejavu.vcf | wc -w)" - 9" | bc)
				NB_VARIANTS=$(grep -cv "^#" $TMP/$GROUP/$PROJECT/dejavu.vcf)
				#NB_SAMPLES=$(echo $($BCFTOOLS view $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz | grep "^#CHROM" | wc -w)" - 9" | bc)
				#NB_VARIANTS=$($BCFTOOLS view $TMP/$GROUP/$PROJECT/dejavu.annotated.vcf.gz | grep -cv "^#")
				
				echo "### Number of samples: $NB_SAMPLES" >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
				echo "### Number of variants: $NB_VARIANTS" >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
				
				# VARIANTS
				echo "### variants frequency" >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
				# for p in $(seq 0 100); do
				#for p in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do awk -F"\t" -v p=$p '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print "# "SUM+0" variants out of "TOT" ("PERC"%) found more than "(p*100)"% in the set"}' $TMP/$GROUP/$PROJECT/dejavu.txt; done >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
				if (($NB_VARIANTS)); then
					for p in $(seq 1 100); do awk -F"\t" -v p=$( echo "scale=2; $p/100" | bc) '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print "# "SUM+0" variants out of "TOT" ("PERC"%) found more than "(p*100)"% in the set"}' $TMP/$GROUP/$PROJECT/dejavu.txt; done >> $TMP/$GROUP/$PROJECT/dejavu.stats.txt
				fi;
				#cat $TMP/$GROUP/$PROJECT/dejavu.stats.txt

				# TSV
				echo "$NB_SAMPLES" >> $TMP/$GROUP/$PROJECT/dejavu.stats.nb_samples
				echo "$NB_VARIANTS" >> $TMP/$GROUP/$PROJECT/dejavu.stats.nb_variants

				echo -e "#NB_variant\tTOT_variant\tPercent\tmore_than" > $TMP/$GROUP/$PROJECT/dejavu.stats.tsv
				#for p in 0.001 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do awk -F"\t" -v p=$p '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print SUM+0"\t"TOT"\t"PERC"\t"(p*100)}' $TMP/$GROUP/$PROJECT/dejavu.txt; done > $TMP/$GROUP/$PROJECT/dejavu.stats.frequency.tsv
				if (($NB_VARIANTS)); then
					for p in $(seq 1 100); do awk -F"\t" -v p=$( echo "scale=2; $p/100" | bc) '$6>p{SUM++} {TOT++} END {PERC=((SUM+0)/TOT)*100; print SUM+0"\t"TOT"\t"PERC"\t"(p*100)}' $TMP/$GROUP/$PROJECT/dejavu.txt; done >> $TMP/$GROUP/$PROJECT/dejavu.stats.tsv
				fi;
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

				cp $TMP/$GROUP/$PROJECT/dejavu.checksum $DEJAVU_FOLDER_VCF/$GROUP/$PROJECT/checksum
				
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


			fi;
		
		fi;

	else

		(($VERBOSE)) && echo "#[INFO] DEJAVU database '$GROUP/$PROJECT' without VCF files"

	fi;

	# CLeaning
	rm -rf $TMP/$GROUP/$PROJECT

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
if (($RELEASE_LATEST)); then
	rm -f $DEJAVU/latest
	ln -s $RELEASE/ $DEJAVU/latest
fi;

for SL in $RELEASE_SYMLINK; do
	rm -f $DEJAVU/$SL
	ln -s $RELEASE/ $DEJAVU/$SL
done;


(($VERBOSE)) && echo "#"
(($VERBOSE)) && echo "#[INFO] DEJAVU database release '$RELEASE' done."
(($VERBOSE)) && echo "$RELEASE" > $DEJAVU/$RELEASE/release

rm -rf $TMP

exit 0;
