##########################################################################
# SETUP Version:			2
# Description:				Setup to configure module
##########################################################################

# DEV version 0.1 : 10/11/2021
# INT version 0.1 : 17/03/2022
# Authoring : Thomas LAVAUX

# PROD version 1 : 17/03/2022
# PROD version 2 : 17/10/2023 : changelog
	# remove os.getenv
	# update AnnotSV to 3.3.6
	# switch to aria2

################## Context ##################
#
# This script will setup a cli ie checking and installing proper databases if needed ; copying launcher.py, .conf and .json files in the specific directory
#
####################################

import os
import subprocess
from datetime import datetime

def systemcall(command):
	""" Execute shell command and return stdout lines as a list """
	process = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return process.stdout.read().decode('utf8').strip().split('\n')

def log_file(logfile, text, sep, items_list=None, items=None):
	""" Function to log a variable value or a list of values into a log file """
	with open(logfile, 'a+') as f:
		f.write(f"{text}{sep}")
		if items_list:
			for item in items_list:
				f.write(f"{str(item) if item != '' else 'None'}{sep}")
		else:
			f.write(f"{str(items) if items != '' else 'None'}{sep}")


def installdatabase(destination, source, archive_name, logfile, errfile):
	""" Function to download and install an archive (zip or tar.gz) """
	os.makedirs(destination, exist_ok = True)
	systemcall(f"aria2c -c -s 16 -x 16 -k 1M -j 1 {source} 1>> {logfile} 2>> {errfile}")
	if archive_name.endswith('.zip'):
		systemcall(f"unzip -q {archive_name} -d {destination} 1>> {logfile} 2>> {errfile}")
	if archive_name.endswith('.tar.gz'):
		systemcall(f"tar xzf {archive_name} -C {destination} 1>> {logfile} 2>> {errfile}")

date_time = datetime.now().strftime("%Y%m%d-%H%M%S")


### INSTALL DATABASES ###

######################
# DATABASE VARIABLES #
######################

DATABASES = "/STARK/databases"
serviceName = "drawfusions"
moduleName = "rnaseq"

#serviceName = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_NAME')
#moduleName = os.getenv('DOCKER_STARK_MODULE_NAME')


#######################
# Copy config files (service.conf & service.json) and launcher.py from app/config to /config/module/servicename/listener/
#######################

if serviceName and moduleName:
	logfile = f"/STARK/config/{moduleName}/{serviceName}/listener/logs/{serviceName}.{date_time}.setup.log"
	errfile = f"/STARK/config/{moduleName}/{serviceName}/listener/logs/{serviceName}.{date_time}.setup.err"
	log_file(logfile, 'Setup copying configuration files:', "\n", items = date_time)
	systemcall(f"rsync -ar /app/config/ /STARK/config/{moduleName}/{serviceName}/listener 1>> {logfile} 2>> {errfile}")
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	log_file(logfile, 'Setup end:', "\n", items = date_time_end)

# SETUPComplete cli services (condition for healthy cli)
systemcall("touch ${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_INNER_FOLDER_SERVICES}/SETUPComplete.txt")