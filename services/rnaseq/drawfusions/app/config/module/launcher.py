##########################################################################
# Launcher Version:			3.0
# Description:				Launcher to run Snakemake module
##########################################################################

# DEV version 0.1 : 10/11/2021
# INT version 0.1 : 17/03/2022
# PROD version 1 : 03/06/2022
# Authoring : Thomas LAVAUX

# PROD version 2 : 16/06/2022 changelog
	# yaml files can be defined by groupe_name + project_name or groupe_name only

# PROD version 3.0 : 28/11/2023 changelog
	# docker compose to run containers

################## Context ##################
# type python launcher.py -h for help
# ex of command: python launcher.py -r run
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

date_time = datetime.now().strftime("%Y%m%d-%H%M%S")

def readconfig(configFile, serviceName, configkey):
	""" Function to extract a specific key configuration variable from a json configuration file"""
	with open(configFile,'r') as f:
		json_config = json.load(f)
	outputconfig = json_config['services'][serviceName][configkey]
	return outputconfig

def launch(run, serviceName, containersFile=None, montage=None, image=None, launchCommand=None, configFile=None, microserviceRepo=None):
	""" Function to start a docker container with a specific command """
	if configFile:
		launchCommand = readconfig(configFile, serviceName, 'launch')
		image = readconfig(configFile, serviceName, 'image')
	if run:
		containerName = f"{serviceName}_{date_time}_{os.path.basename(run)}"
		group_name = run.split('/')[4]
		project_name = run.split('/')[5]
	
	if group_name and project_name:
		yaml_path = f"{os.getenv('DOCKER_STARK_MODULE_SUBMODULE_INNER_FOLDER_CONFIG')}/cli"
		yaml_config_file = f"{yaml_path}/{group_name}_{project_name}.yaml"
		if not os.path.exists(yaml_config_file):
			yaml_config_file = f"{yaml_path}/{group_name}.yaml"
	else:
		yaml_config_file = None

	COMPOSE_PATH = f"{os.getenv('DOCKER_STARK_MODULE_SUBMODULE_INNER_FOLDER_CONFIG')}/listener/"
	if yaml_config_file and os.path.exists(yaml_config_file):
		cmd = f"docker compose -f {COMPOSE_PATH}/STARK.docker-compose.yml run --rm --name={containerName} {image} '{launchCommand} --config run={run} --configfile {yaml_config_file}'"
		print(cmd)
	else:
		cmd = f"docker compose -f {COMPOSE_PATH}/STARK.docker-compose.yml run --rm --name={containerName} {image} '{launchCommand} --config run={run}'"
		print(cmd)
	subprocess.call(cmd, shell = True)
	
def myoptions():
	'''
	*arg parser*
	*return options*
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("-r", "--run", type = str, default = "", help = "Run to launch", dest = 'run')
	parser.add_argument("-s", "--service", type = str, default = "", help = "Service name", dest = 'serviceName')
	parser.add_argument("-log", "--log", type = str, default = "", help = "Path for the log file", dest = 'containersFile')
	parser.add_argument("-v", "--volume", type = str, default = "", help = "Docker volumes to add", dest = 'montage')
	parser.add_argument("-i", "--image", type = str, default = "", help = "Docker image to use", dest = 'image')
	parser.add_argument("-l", "--launchcommand", type = str, default = "", help = "Command to launch inside the container", dest = 'launchCommand')
	parser.add_argument("-c", "--config", type = str, default = "", help = "Config file to read from", dest = 'configFile')
	parser.add_argument("-repo", "--repo", type = str, default = "", help = "Microservice repository name", dest = 'microserviceRepo')
	return parser.parse_args()

if __name__ == "__main__":
	doctest.testmod()
	args = myoptions()
	launch(args.run, args.serviceName, args.containersFile, args.montage, args.image, args.launchCommand, args.configFile, args.microserviceRepo)
