#! /usr/bin/env python
# -*- coding: utf-8 -*-

# launcher for DECoN
#from __future__ import division
#from __future__ import print_function

import os
import re
import subprocess
import json

from datetime import datetime
from os.path import join as osj

def createContainerFile(containersFile, run, containerName):
	file = open(osj(containersFile,containerName+".log"), "w+")
	file.write("RUN: "+os.path.basename(run)+"\n")
	file.write("FOLDER: "+run+"\n")
	file.write("EXEC_DATE: "+datetime.now().strftime("%d%m%Y-%H%M%S")+"\n")
	file.write("ID: "+containerName+"\n")
	file.close()


# From listener.py == launch(run, serviceName, containersFile, os.getenv('MICROSERVICE_MONTAGE'), os.getenv('MICROSERVICE_IMAGE'), os.getenv('MICROSERVICE_LAUNCH'), configFile)
# montage, image and launchCommand came from the .env
# montage = MICROSERVICE_MONTAGE = ${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_PARAM_CONTAINER_MOUNT}
# images = MICROSERVICE_IMAGE = DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_IMAGE = =${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_IMAGE_REPOSITORY}${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_IMAGE_NAME}:${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_IMAGE_TAG}
# launchCommand = MICROSERVICE_LAUNCH = DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_CONTAINER_COMMAND = "bin/canoes run -g ${DOCKER_STARK_INNER_FOLDER_DATABASES}/genomes/current/hg19.fa"

# serviceName is the name of the module (ex DECON)
# run from listener.py (STARKComplexte.txt detection) = complete run directory ; /STARK/output/repository/GROUP/APP/{run} ie ${DOCKER_STARK_MAIN_FOLDER}/${DOCKER_STARK_FOLDER_OUTPUT_REPOSITORY}/GROUP/APP/{nameoftherun}
# containersFile is a log file created by the launcher.py
# configFile is the path to the config file module.conf containing also "image" and "launchCommand"


# /home1/data/STARK/config/structuralvariation/decon/listener/decon.conf
# configFile is the --config argument from listener.py, defined in the decon.env file
# we get image and launchCommand from the decon.conf file
# os.getenv('MICROSERVICE_MONTAGE') == DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_PARAM_JSON
configFile = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_PARAM_JSON')
json_config = json.load(configFile)
launchCommand = getDataFromJson(json_config)['services'][serviceName]['launch']
image = getDataFromJson(json_config)['services'][serviceName]['images']

# Add the default.yaml config to the cmd and specific yaml file depending on the APP and/or GROUP tag
# Find the GROUP or AP tag and then add an GROUP.yaml or a APP.yaml to the cmd after --config
# get group [4] and app [5] name from run
# /stark/output/repository/group/app/run
# run structure is  "/STARK/output/repository/group/app/run"
try:
	group_name =  run.split('/')[4]
	app_name = run.split('/')[5]
except IndexError: pass

# Path of the yaml config file
# (/home1/data)/STARK/config/structuralvariation/decon/cli/*.yaml
# option is only for GROUP NAME
yaml_path = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_INNER_FOLDER_CONFIG')

if group_name and yaml_path:
	yaml_config_file = yaml_path+"/"+group_name+".yaml"	
	if not os.path.exists(yaml_config_file) and os.path.getsize(yaml_config_file) == 0:
		yaml_config_file = yaml_path+"/default.yaml"
else:
	yaml_config_file = ""


# snakefile command options after --config should be separated by a space
# to exclude sample = EXCLUDE_SAMPLE_LIST=["SGT2101772","SGT2101762"]
# exemple : snakemake -s snakefile_decon -c1 --config run=/STARK/data/run EXCLUDESAMPLE_LIST=["SGT2101769","SGT2101772"] --configfile default.yaml

def launch(run, serviceName, containersFile, montage, image, launchCommand, configFile):
	containerName = serviceName+os.path.basename(run)
	cli_service_name = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_CONTAINER_NAME')
	# Should be stark-module-structuralvariation-submodule-decon-service-cli
	# docker-compose -f to specify an alternate compose file (default: docker-compose.yml) ; --env-file PATH to specify an alternate environment file
	# we use docker-compose run to launch an analysis, the -f argument can be use to specify a specific servicename.docker-compose.yml
	# we can pass arguements like a classic docker run
	cmd = "docker-compose run --rm -f "+serviceName+"docker-compose.yml " +cli_service_name+" "+launchCommand+" --config DATA_DIR="+os.path.basename(run)+" --configfile "+yaml_config_file
	subprocess.call(cmd, shell = True)
	# Create a log file
	createContainerFile(containersFile, run, containerName)