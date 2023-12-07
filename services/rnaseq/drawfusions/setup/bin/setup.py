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

###############
# Python Func #
###############

import os
import subprocess
from datetime import datetime


def systemcall(command):
	'''
	*passing command to the shell*
	*return list containing stdout lines*
	command - shell command (string)
	'''
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read().decode('utf8').strip().split('\n')

def logsomefile(logfile, text, sep, items_list=None, items=None):
	""" Function to log variable value or list values into a log file """
	pathlog = os.path.dirname(logfile)
	os.makedirs(pathlog, exist_ok = True)
	with open(logfile, 'a+') as f:
		f.write(text + sep)
		if items_list:
			for items in items_list:
				f.write(str(items) + sep)
		else:
			f.write(str(items) + sep)

def installdatabase(destination, source, archive_name, logfile, errfile):
	""" Function to download and install an archive (zip or tar.gz) """
	os.makedirs(destination, exist_ok = True)
	systemcall("aria2c -c -s 16 -x 16 -k 1M -j 1 "+ source+" 1>> "+logfile+" 2>> "+errfile+" ")
	if archive_name.endswith('.zip'):
		systemcall("unzip -q "+archive_name+" -d " +destination+" 1>> "+logfile+" 2>> "+errfile+" ")
	if archive_name.endswith('.tar.gz'):
		systemcall("tar xzf "+archive_name+" -C "+destination+" 1>> "+logfile+" 2>> "+errfile+" ")

# for wget
# systemcall("wget "+ source+" 1>> "+logfile+" 2>> "+errfile+" ")

# Variables initialisation
# set datetime to add to output file name
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
	logfile = "/STARK/config/"+moduleName+"/"+serviceName+"/listener/logs/" + serviceName + "." +date_time+".setup.log"
	errfile = "/STARK/config/"+moduleName+"/"+serviceName+"/listener/logs/" + serviceName + "." +date_time+".setup.err"
	logsomefile(logfile, 'Setup copying configuration files:', "\n", items = date_time)
	systemcall("rsync -ar /app/config/ /STARK/config/"+moduleName+"/"+serviceName+"/listener 1>> "+logfile+" 2>> "+errfile+" ")
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	logsomefile(logfile, 'Setup end:', "\n", items = date_time_end)

# SETUPComplete cli services (condition for healthy cli)
systemcall("touch ${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_INNER_FOLDER_SERVICES}/SETUPComplete.txt")