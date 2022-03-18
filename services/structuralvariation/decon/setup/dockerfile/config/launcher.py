##########################################################################
# Launcher Version:			0.1
# Description:				Launcher to run Snakemake module
##########################################################################

# DEV version 0.1 : 10/11/2021
# INT version 0.1 : 17/03/2022
# Authoring : Thomas LAVAUX

################## Context ##################
# launch command exemple (-p to display shell command)
# snakemake -p -s /app/bin/snakefile_decon -c1 --config run=/STARK/output/repository/GROUP/APP/run
####################################

# From listener.py
# launch(run, serviceName, containersFile, os.getenv('MICROSERVICE_MONTAGE'), os.getenv('MICROSERVICE_IMAGE'), os.getenv('MICROSERVICE_LAUNCH'), configFile, os.getenv('MICROSERVICE_REPOSITORY'))

# run = path of the run to analyse
# serviceName = name of the service (ex DECON)
# containersFile = part name of a log file
# configFile = config json with .conf ext containing launch command (= MICROSERVICE_LAUNCH) and image name of the docker ( = MICROSERVICE_IMAGE)

# MICROSERVICE_MONTAGE = additional volume to mount
# MICROSERVICE_IMAGE = cf. supra
# MICROSERVICE_LAUNCH = cf. supra
# MICROSERVICE_REPOSITORY = path to the repository folder

	# Add the default.yaml config to the cmd and specific yaml file depending on the APP and/or GROUP tag
	# Find the GROUP or AP tag and then add an GROUP.yaml or a APP.yaml to the cmd after --config
	# get group [4] and app [5] name from rundir == /stark/output/repository/group/app/run
	# Path of the yaml config file == /home1/data)/STARK/config/structuralvariation/decon/cli/*.yaml
	# option is only for GROUP NAME

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

################## Import libraries ##################

import os
import re
import subprocess
import json
import argparse
import doctest

from datetime import datetime
from os.path import join as osj

### FUNCTIONS ###

def createContainerFile(containersFile, run, containerName):
	""" Function to create a log file """
	file = open(osj(containersFile,containerName+".log"), "w+")
	file.write("RUN: "+os.path.basename(run)+"\n")
	file.write("FOLDER: "+run+"\n")
	file.write("EXEC_DATE: "+datetime.now().strftime("%d%m%Y-%H%M%S")+"\n")
	file.write("ID: "+containerName+"\n")
	file.close()

def readconfig(configfile, serviceName, configkey):
	""" Function to extract a specific key configuration variable from a json configuration file"""
	with open(configFile,'r') as f:
		json_config = json.load(f)
	outputconfig = json_config['services'][serviceName][configkey]
	return outputconfig

def getsomepathname(path, number, sep):
	""" Function to extract variable from a string with a separator """
	try:
		output = path.split(sep)[number]
	except IndexError: pass
	return output

# need to test the default argument with the listener.py
def launch(run, serviceName, configFile, image=None, launchCommand=None, montage=None, containersFile=None):
	""" Function to start a docker container with a specific command """
	# Variables get from env
	if not configFile:
		configFile = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_PARAM_CONF')
	if not serviceName:
		serviceName = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_PARAM_MICROSERVICE_NAME')
	if not yaml_path:
		yaml_path = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_INNER_FOLDER_CONFIG')
	if not image:
		image = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_CONTAINER_NAME')
	# Variables get from config file
	if configFile:
		launchCommand = readconfig(configfile, serviceName, 'launch')
		image = readconfig(configfile, serviceName, 'image')
	if serviceName and run:
		containerName = serviceName + os.path.basename(run)
		group_name = getsomepathname(run, 4, '/')
		app_name = getsomepathname(run, 5, '/')
	if yaml_path:
		yaml_config_file = yaml_path+"/default.yaml"
	if group_name and yaml_path:
		yaml_config_file = yaml_path+"/"+group_name+".yaml"
	# Construct the docker-compose command
	# --name assign a name to the container 
	# -v, --volume=[] Bind mount a volume (default []) can be empty
	# --configfile can't be empty ; default is yaml_path+"/default.yaml" ; if not set it will use the default yaml file in the docker container
	if yaml_path:
		cmd = "docker-compose -f " +serviceName+ ".docker-compose.yml run -rm -v "+ montage + " " +image+ " " +launchCommand+ " --config run=" +run+ " --configfile " +yaml_config_file
	else:
		cmd = "docker-compose -f " +serviceName+ ".docker-compose.yml run -rm -v "+ montage + " " +image+ " " +launchCommand+ " --config run=" +run
	# launch the cmd to the shell 
	subprocess.call(cmd, shell = True)
	# Create a log file
	createContainerFile(containersFile, run, containerName)


def myoptions():
	'''
	*arg parser*
	*return options*
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("-r", "--run", type = str, default = " ", help = "Run to launch", dest = 'run')
	parser.add_argument("-s", "--service", type = str, default = " ", help = "Service name", dest = 'serviceName')
	parser.add_argument("-log", "--log", type = str, default = " ", help = "Log file", dest = 'containersFile')
	parser.add_argument("-v", "--volume", type = str, default = " ", help = "Docker volume to add", dest = 'montage')
	parser.add_argument("-i", "--image", type = str, default = " ", help = "Docker image to use", dest = 'image')
	parser.add_argument("-l", "--launchcommand", type = str, default = " ", help = "Command to launch inside the container", dest = 'launchCommand')
	parser.add_argument("-c", "--config", type = str, default = " ", help = "Config file to read from", dest = 'configFile')
	return parser.parse_args()

if __name__ == "__main__":
	doctest.testmod()
	args = myoptions()
	launch(args.run, args.serviceName, args.containersFile, args.montage, args.image, args.launchCommand, args.configFile)
