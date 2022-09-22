##########################################################################
# Launcher Version:			0.1
# Description:				Launcher to run Snakemake module
##########################################################################

# DEV version 0.1 : 10/11/2021
# INT version 0.1 : 17/03/2022
# Authoring : Thomas LAVAUX

# yaml files are groupe_name and project_name dependant

################## Context ##################
# type python launcher.py -h for help
#
# ex of command :python launcher.py -r run
#
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

################## Import libraries ##################

import os
import re
import subprocess
import json
import argparse
import doctest

from datetime import datetime
from os.path import join as osj

# Variables
# set datetime to add to container name
date_time = datetime.now().strftime("%Y%m%d-%H%M%S")

### FUNCTIONS ###

def createlog(containersFile, run, containerName):
	""" Function to create a log file """
	file = open(osj(containersFile,containerName+".log"), "w+")
	file.write("RUN: "+os.path.basename(run)+"\n")
	file.write("FOLDER: "+run+"\n")
	file.write("EXEC_DATE: "+datetime.now().strftime("%d%m%Y-%H%M%S")+"\n")
	file.write("ID: "+containerName+"\n")
	file.close()

def readconfig(configFile, serviceName, configkey):
	""" Function to extract a specific key configuration variable from a json configuration file"""
	with open(configFile,'r') as f:
		json_config = json.load(f)
	outputconfig = json_config['services'][serviceName][configkey]
	return outputconfig

def launch(run, serviceName, containersFile, montage, image, launchCommand, configFile, microserviceRepo):
	""" Function to start a docker container with a specific command """
	""" See help (-h) for details """
	# Variables get from env
	if not configFile:
		configFile = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_PARAM_CONF')
	if not microserviceRepo:
		microserviceRepo = os.getenv('MICROSERVICE_REPOSITORY')
	if not serviceName:
		serviceName = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_PARAM_MICROSERVICE_NAME')
	if not image:
		image = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_CONTAINER_NAME')
	if not montage:
		montage = os.getenv('MICROSERVICE_MONTAGE')
	# Variables get from config file
	if configFile:
		launchCommand = readconfig(configFile, serviceName, 'launch')
		image = readconfig(configFile, serviceName, 'image')
	if run:
		containerName = serviceName + "_" +date_time+"_"+os.path.basename(run)
		group_name = run.split('/')[4]
		project_name = run.split('/')[5]
	# Config snakefile config file path
	yaml_path = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_INNER_FOLDER_CONFIG')
	if group_name and project_name:
		yaml_config_file = yaml_path+"/"+group_name+"_"+project_name+".yaml"
		if not yaml_config_file:
			yaml_config_file = yaml_path+"/"+group_name+".yaml"
	else:
		yaml_config_file = None
	# Construct the docker command
	# snakemake --configfile can't be empty ; if not, it will use the default yaml file in the docker container
	# /bin/bash -c 'source activate variantconvert is mandatory to activate the conda environment, you have to run bash
	if yaml_config_file and os.path.exists(yaml_config_file):
		cmd = "docker run --rm --name="+containerName+" "+montage+" "+image+" /bin/bash -c 'source activate variantconvert && "+launchCommand+" --config run="+run+" --configfile "+yaml_config_file+"'"
	else:
		cmd = "docker run --rm --name="+containerName+" "+montage+" "+image+" /bin/bash -c 'source activate variantconvert && "+launchCommand+" --config run="+run+"'"
	# launch the cmd to the shell 
	subprocess.call(cmd, shell = True)
	# Create a log file
	if not containersFile:
		containersFile = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_INNER_FOLDER_SERVICES')
	if os.path.exists(containersFile):
		createlog(containersFile, run, containerName)


def myoptions():
	'''
	*arg parser*
	*return options*
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("-r", "--run", type = str, default = "", help = "Run to launch", dest = 'run')
	parser.add_argument("-s", "--service", type = str, default = "", help = "Service name", dest = 'serviceName')
	parser.add_argument("-log", "--log", type = str, default = "", help = "Path for the log file", dest = 'containersFile')
	parser.add_argument("-v", "--volume", type = str, default = "", help = "Docker volume to add", dest = 'montage')
	parser.add_argument("-i", "--image", type = str, default = "", help = "Docker image to use", dest = 'image')
	parser.add_argument("-l", "--launchcommand", type = str, default = "", help = "Command to launch inside the container", dest = 'launchCommand')
	parser.add_argument("-c", "--config", type = str, default = "", help = "Config file to read from", dest = 'configFile')
	parser.add_argument("-repo", "--repo", type = str, default = "", help = "Microservice repository name", dest = 'microserviceRepo')
	return parser.parse_args()

if __name__ == "__main__":
	doctest.testmod()
	args = myoptions()
	launch(args.run, args.serviceName, args.containersFile, args.montage, args.image, args.launchCommand, args.configFile, args.microserviceRepo)
