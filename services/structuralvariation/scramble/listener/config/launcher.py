#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import os
import re
import subprocess
import json

from datetime import datetime
from os.path import join as osj

### Start FUNCTIONS ###

# containersFile is a log file created by the launcher.py
def createContainerFile(containersFile, run, containerName):
	file = open(osj(containersFile,containerName+".log"), "w+")
	file.write("RUN: "+os.path.basename(run)+"\n")
	file.write("FOLDER: "+run+"\n")
	file.write("EXEC_DATE: "+datetime.now().strftime("%d%m%Y-%H%M%S")+"\n")
	file.write("ID: "+containerName+"\n")
	file.close()

# From listener.py == launch(run, serviceName, containersFile, os.getenv('MICROSERVICE_MONTAGE'), os.getenv('MICROSERVICE_IMAGE'), os.getenv('MICROSERVICE_LAUNCH'), configFile)
# run from listener.py (STARKComplexte.txt detection) = complete run directory ; /STARK/output/repository/GROUP/APP/{run} ie ${DOCKER_STARK_MAIN_FOLDER}/${DOCKER_STARK_FOLDER_OUTPUT_REPOSITORY}/GROUP/APP/{nameoftherun}
# we get montage, image and launchCommand from the .env or the the decon.conf file
# launchCommand = MICROSERVICE_LAUNCH = DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_CONTAINER_COMMAND
# images = MICROSERVICE_IMAGE = DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_IMAGE =${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_IMAGE_REPOSITORY}${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_IMAGE_NAME}:${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_IMAGE_TAG}
# montage = MICROSERVICE_MONTAGE = ${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_PARAM_CONTAINER_MOUNT}
# serviceName is the name of the module (ex DECON)
# here we get image and launchCommand from decon.conf file == configFile
# for ex /home1/data/STARK/config/structuralvariation/decon/listener/decon.conf

# Variables get from env
configFile = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_PARAM_CONF')
serviceName = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_PARAM_MICROSERVICE_NAME')

with open(configFile,'r') as f:
	json_config = json.load(f)
launchCommand = json_config['services'][serviceName]['launch']
image = json_config['services'][serviceName]['image']

def launch(run, serviceName, containersFile, montage, image, launchCommand, configFile):
	containerName = serviceName+os.path.basename(run)
	cli_service_name = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_CONTAINER_NAME')
	# run structure is  "/STARK/output/repository/group/app/run"
	try:
		group_name = run.split('/')[4]
		app_name = run.split('/')[5]
	except IndexError: pass
	# Add the default.yaml config to the cmd and specific yaml file depending on the APP and/or GROUP tag
	# Find the GROUP or AP tag and then add an GROUP.yaml or a APP.yaml to the cmd after --config
	# get group [4] and app [5] name from rundir == /stark/output/repository/group/app/run
	# Path of the yaml config file == /home1/data)/STARK/config/structuralvariation/decon/cli/*.yaml
	# option is only for GROUP NAME
	yaml_path = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_INNER_FOLDER_CONFIG')
	yaml_config_file = yaml_path+"/default.yaml"
	if group_name and yaml_path:
		yaml_config_file = yaml_path+"/"+group_name+".yaml"	
	#	if not os.path.exists(yaml_config_file) and os.path.getsize(yaml_config_file) == 0:
	#		yaml_config_file = yaml_path+"/default.yaml"
	# snakefile command options after --config should be separated by a space
	# to exclude sample = EXCLUDE_SAMPLE_LIST=["SGT2101772","SGT2101762"]
	# exemple : snakemake -s snakefile_decon -c1 --config run=/STARK/data/run EXCLUDE_SAMPLE_LIST=["SGT2101769","SGT2101772"] --configfile default.yaml
	# docker-compose -f to specify an alternate compose file (default: docker-compose.yml) ; --env-file PATH to specify an alternate environment file
	# Construct the docker-compose command
	# we use docker-compose run to launch an analysis, the -f argument can be use to specify a specific servicename.docker-compose.yml
	# we can pass arguements like a classic docker run
	# docker-compose -f test.docker-compose.yml  run  stark-module-structuralvariation-submodule-decon-service-listener
	cmd = "docker-compose -f "+serviceName+"docker-compose.yml run -rm " +cli_service_name+" "+launchCommand+" --config run="+os.path.basename(run)+" --configfile "+yaml_config_file
	# launch the cmd to the shell 
	subprocess.call(cmd, shell = True)
	# Create a log file
	createContainerFile(containersFile, run, containerName)