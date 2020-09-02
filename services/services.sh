#!/bin/bash
#################################
##
## STARK Services Modules 
##
#################################

SCRIPT_NAME="STARKServicesModules"
SCRIPT_DESCRIPTION="STARK Services Modules"
SCRIPT_RELEASE="0.9b"
SCRIPT_DATE="18/05/2020"
SCRIPT_AUTHOR="Antony Le Bechec"
SCRIPT_COPYRIGHT="HUS/CPS"
SCRIPT_LICENCE="GNU GPLA V3"

# Realse note
RELEASE_NOTES=$RELEASE_NOTES"# 0.9b-8/05/2020: Script creation\n";

# Script folder
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


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
	echo "# USAGE: $(basename $0) [-h] [options...]";
	echo "#";
	echo "### This script setup STARK envronment.";
	echo "### Can be used with curl (e.g. curl https://gitlab.bioinfo-diag.fr/Strasbourg/vision/raw/master/setup.sh | bash)";
	echo "###    (e.g. curl https://gitlab.bioinfo-diag.fr/Strasbourg/STARK/raw/master/setup.sh | bash)";
	echo "#";
  	echo "# --env=<FILE>                 STARK Docker Compose environment file ";
	echo "#                              Default: '.env'";
  	echo "# --folder_services=<FOLDER>   STARK services modules folder ";
	echo "#                              Default: in env file, or 'services'";
  	echo "# --modules=<STRING>           STARK modules as a list";
	echo "#                              Format: 'module1,module2,...', use '*' as wildcard";
	echo "#                              Default: '*' all modules";
	echo "# --services=<STRING>          STARK modlule services as a list";
	echo "#                              Only for command: up, start, stop, restart'";
	echo "#                              Format: 'module1,module2,...'";
	echo "#                              Default: '*' all services";
	echo "# --command=<STRING>           Docker Compose Command ";
	echo "#                              - up: Create and start containers (in daemon mode)";
	echo "#                              - down: Stop and remove containers, networks, images, and volumes";
	echo "#                              - start: Start services";
	echo "#                              - stop: Stop services";
	echo "#                              - restart: Restart services";
	echo "#                              - config: Check config services";
	echo "#                          Default: 'up'";
	echo "# -v|--verbose             Verbose mode";
	echo "# -d|--debug               Debug mode";
	echo "# -n|--release             Script Release";
	echo "# -h|--help                Help message";
	echo "#";
	#echo -e "#\n# RUN Analysis\n################";
	#$STARK_FOLDER_BIN/launch.sh -h | grep "# [ |-]";
	#echo -e "#\n# SAMPLE Analysis\n###################";
	#$STARK_FOLDER_BIN/launch.sample.sh -h | grep "# [ |-]";
}




####################################################################################################################################
# Getting parameters from the input
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ":" tells that the option has a required argument, "::" tells that the option has an optional argument, no ":" tells no argument
ARGS=$(getopt -o "vdnh" --long "env:,folder_services:,modules:,services:,command:,verbose,debug,release,help" -- "$@" 2> /dev/null)
if [ $? -ne 0 ]; then
	:
	echo "#[ERROR] Error in the argument list:";
	echo "#[ERROR] $@"
	echo ""
	usage;
	exit;
fi;




PARAM=$@
DEBUG=0
VERBOSE=0

eval set -- "$ARGS"
while true
do
	#echo "$1=$2"
	#echo "Eval opts";
	case "$1" in
		--env)
			ENV="$2"
			shift 2
			;;
    --folder_services)
			FOLDER_SERVICES="$2"
			shift 2
			;;
    --modules)
			MODULES="$2"
			shift 2
			;;
    --services)
			SERVICES="$2"
			shift 2
			;;
    --command)
			COMMAND="$2"
			shift 2
			;;
  	-h|--help)
			usage
			exit 0
			;;
		-v|--verbose)
			VERBOSE=1
			shift 1
			;;
		-d|--debug)
			DEBUG=1
			shift 1
			;;
		--) shift
			break
			;;
		*) 	echo "Option $1 is not recognized. " "Use -h or --help to display the help." && \
			exit 1
			;;
	esac
done

# header
(($NO_HEADER)) || header;



####################################################################################################################################
# Checking the input parameter
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# if [ -z "$SOURCE_RUNS" ] && [ -z "$DEST_RUNS" ] && ((!$DEBUG)); then
# 	echo "#[ERROR] Required parameter: --sources and --dest. Use --help to display the help." && echo "" && usage && exit 1;
# fi
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


(($DEBUG)) && VERBOSE=1


## PARAM
##########

# Copy Command
#COMMAND_COPY="rsync -auczqAXhi --no-links --no-perms --no-owner --no-group"
COMMAND_COPY="rsync -az --update"


# TMP
TMP_FOLDER=/tmp/STARK/services.$RANDOM
mkdir -p $TMP_FOLDER


# DOCKER
DOCKER_VERSION=$(docker --version)
(($VERBOSE)) && echo "#[INFO] Docker version '$DOCKER_VERSION'"


# ENV
if [ -z "$ENV" ]; then
  if [ -e "$SCRIPT_DIR/../.env" ]; then
	   ENV="$SCRIPT_DIR/../.env";
  elif [ -e "$SCRIPT_DIR/.env" ]; then
    ENV="$SCRIPT_DIR/.env";
  elif [ -e ".env" ]; then
  	ENV=".env";
  else
    echo "#[ERROR] No STARK Docker Compose environment file"
    usage;
    exit 1;
  fi;
fi;

(($VERBOSE)) && echo "#[INFO] STARK Docker Compose environment file '$ENV'"

source $ENV;


# DOCKER_STARK_MAIN_FOLDER
if [ -z "$DOCKER_STARK_MAIN_FOLDER" ]; then
  if [ "$DOCKER_STARK_MAIN_FOLDER" != "" ] && [ -d $DOCKER_STARK_MAIN_FOLDER ]; then
    DOCKER_STARK_MAIN_FOLDER=$DOCKER_STARK_MAIN_FOLDER
  elif [ -d "~/STARK" ]; then
    DOCKER_STARK_MAIN_FOLDER="~/STARK"
  else
    echo "#[ERROR] STARK Main folder"
    exit 1;
  fi;
fi;

(($VERBOSE)) && echo "#[INFO] STARK main folder '$DOCKER_STARK_MAIN_FOLDER'"


# SERVICES
if [ -z "$FOLDER_SERVICES" ]; then
  if [ "$DOCKER_STARK_FOLDER_SERVICES" != "" ]; then
    FOLDER_SERVICES=$DOCKER_STARK_MAIN_FOLDER/$DOCKER_STARK_FOLDER_SERVICES
  else
  	FOLDER_SERVICES=$DOCKER_STARK_MAIN_FOLDER/services
  fi;
fi;

mkdir -p $FOLDER_SERVICES

(($VERBOSE)) && echo "#[INFO] STARK services modules folder '$FOLDER_SERVICES'"


# MODULE
if [ -z "$MODULES" ]; then
  MODULES='*'
fi;


# SERVICES
if [ -z "$SERVICES" ]; then
  SERVICES=''
fi;



# COMMAND
if [ -z "$COMMAND" ]; then
  COMMAND='up'
fi;

[ "$COMMAND" == "up" ] && COMMAND="up -d"
[ "$COMMAND" == "config" ] && SERVICES=""

(($VERBOSE)) && echo "#[INFO] STARK Docker Compose Command '$COMMAND'"


### Find services

for service_module in $(ls -d $SCRIPT_DIR/$MODULES 2>/dev/null); do

	if [ -d $service_module ]; then

		(($DEBUG)) && echo "#[INFO] Folder '$service_module'"

		# Module name
		module_name=$(basename $service_module)

		# Find YML file
		service_module_yml="";
		if [ -e "$service_module/STARK.docker-compose.yml" ]; then
			service_module_yml="$service_module/STARK.docker-compose.yml";
		elif [ -e "$service_module/STARK.yml" ]; then
			service_module_yml="$service_module/STARK.yml";
		elif [ -e "$service_module/docker-compose.yml" ]; then
			service_module_yml="$service_module/docker-compose.yml";
		else
			(($DEBUG)) && echo "#[WARNING] Folder '$module_name' does not have docker configuration file"
		fi;

		# Find ENV file
		service_module_env="";
		if [ -e "$service_module/STARK.env" ]; then
			service_module_env="$service_module/STARK.env";
		elif [ -e "$service_module/.env" ]; then
			service_module_env="$service_module/.env";
		else
			(($DEBUG)) && echo "#[WARNING] Folder '$module_name' does not have environment file"
		fi;


		# Find MODULE file
		service_module_module="";
		if [ -e "$service_module/STARK.module" ]; then
			service_module_module="$service_module/STARK.module";
		elif [ -e "$service_module/.module" ]; then
			service_module_module="$service_module/.module";
		else
			(($DEBUG)) && echo "#[WARNING] Folder '$module_name' does not have module configuration file"
		fi;


		# Find README file
		service_module_readme="";
		if [ -e "$service_module/STARK.README.md" ]; then
			service_module_readme="$service_module/STARK.README.md";
		elif [ -e "$service_module/README.md" ]; then
			service_module_readme="$service_module/README.md";
		else
			(($DEBUG)) && echo "#[WARNING] Folder '$module_name' does not have README file"
		fi;



		if [ "$service_module_yml" != "" ] && [ "$service_module_env" != "" ] && [ "$service_module_module" != "" ]; then

			(($VERBOSE)) && echo "#[INFO] Service module '$module_name'"
			(($VERBOSE)) && echo "#[INFO] Folder '$module_name' docker configuration file: '$service_module_yml'"
			(($VERBOSE)) && echo "#[INFO] Folder '$module_name' environment file: '$service_module_env'"
			(($VERBOSE)) && echo "#[INFO] Folder '$module_name' module configuration file: '$service_module_module'"
			(($VERBOSE)) && echo "#[INFO] Folder '$module_name' README file: '$service_module_readme'"

			# Alternative ENV
			(($DEBUG)) && echo "#[INFO] Service module '$module_name' - Alternate environment file creation"
			cat $ENV 1> $TMP_FOLDER/.env 2>$TMP_FOLDER/.env.err
			cat $service_module_env 1>> $TMP_FOLDER/.env 2>>$TMP_FOLDER/.env.err
			# Module configuration file
			(($DEBUG)) && echo "#[INFO] Service module '$module_name' - Module configuration file copy"
			#echo $COMMAND_COPY $service_module_module $SERVICES/$module_name/;
			$COMMAND_COPY $service_module_module $FOLDER_SERVICES/$module_name/;
			(($DEBUG)) && echo "#[INFO] Service module '$module_name' - Readme file copy"
			if [ "$service_module_readme" != "" ]; then
				$COMMAND_COPY $service_module_readme $FOLDER_SERVICES/$module_name/ 2>/dev/null;
			fi;
			
			# DEBUG
			(($DEBUG)) && cat $TMP_FOLDER/.env $TMP_FOLDER/.env.err
			(($DEBUG)) && echo "    docker-compose --file $service_module_yml --env-file $TMP_FOLDER/.env -p $module_name $COMMAND"
			(($DEBUG)) && docker-compose --file $service_module_yml --env-file $TMP_FOLDER/.env -p $module_name config


			# patch for --env-file unavailable
			if (( $(docker-compose --help | grep "\-\-env\-file" -c) )); then
				
				# Command
				if docker-compose --file $service_module_yml --env-file $TMP_FOLDER/.env -p $module_name config 1>$TMP_FOLDER/docker-compose.log 2>$TMP_FOLDER/docker-compose.err; then
				    docker-compose --file $service_module_yml --env-file $TMP_FOLDER/.env -p $module_name $COMMAND $SERVICES
				else
					echo "#[ERROR] docker-compose error";
					cat $TMP_FOLDER/docker-compose.log $TMP_FOLDER/docker-compose.err;
					exit 1;
				fi;

			else

				(($VERBOSE)) &&  echo "#[WARNING] Docker version does not provide '--env-file' parameter. STARK module is patched to work, but please update your docker version."

				# Command
				if env $(cat $TMP_FOLDER/.env | grep "#" -v) docker-compose --file $service_module_yml -p $module_name config 1>$TMP_FOLDER/docker-compose.log 2>$TMP_FOLDER/docker-compose.err; then
				    env $(cat $TMP_FOLDER/.env | grep "#" -v) docker-compose --file $service_module_yml -p $module_name $COMMAND $SERVICES
				else
					echo "#[ERROR] docker-compose error";
					cat $TMP_FOLDER/docker-compose.log $TMP_FOLDER/docker-compose.err;
					exit 1;
				fi;

			fi;
		else
			(($VERBOSE)) && echo "#[INFO] Folder '$module_name'  is not a Service module (no STARK.docker-compose.yml, STARK.env and STARK.module)"
		fi;

	fi;

done;

# TMP
rm -rf $TMP_FOLDER


exit 0
