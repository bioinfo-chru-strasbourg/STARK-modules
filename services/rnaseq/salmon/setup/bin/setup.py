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
serviceName = "salmon"
moduleName = "rnaseq"

#serviceName = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_NAME')
#moduleName = os.getenv('DOCKER_STARK_MODULE_NAME')

###############
# GENCODE DB  #
###############

# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz

# Version must be the same for the gencode's db
GENCODE_VERSION = "41"

GENCODEREF_TARBALL = "GRCh38.primary_assembly.genome.fa.gz"
GENCODEREF_SOURCE_EXTERNAL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"+GENCODEREF_VERSION+"/"+GENCODEREF_TARBALL
GENCODEREF_PARAM_DATABASE_FOLDER_LINK = DATABASES+"/gencode/"+GENCODEREF_VERSION+"/"

logfile = GENCODEREF_PARAM_DATABASE_FOLDER_LINK + serviceName + "." +date_time+".database.setup.log"
errfile = GENCODEREF_PARAM_DATABASE_FOLDER_LINK + serviceName + "." +date_time+".database.setup.err"

# Check if a directory exist
# os.path.isdir('folder') will return true if exist
if os.path.isdir(GENCODEREF_PARAM_DATABASE_FOLDER_LINK) == False:
	os.makedirs(GENCODEREF_PARAM_DATABASE_FOLDER_LINK, exist_ok = True)
	logsomefile(logfile, 'Gencode version '+GENCODE_VERSION+' installation start:', "\n", items = date_time)
	installdatabase(GENCODEREF_PARAM_DATABASE_FOLDER_LINK, GENCODEREF_SOURCE_EXTERNAL, GENCODEREF_TARBALL, logfile, errfile)
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	logsomefile(logfile, 'Installation end:', "\n", items = date_time_end)


GENCODE_TARBALL = "gencode.v"+GENCODE_VERSION+".transcripts.fa.gz"
GENCODE_SOURCE_EXTERNAL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"+GENCODE_VERSION+"/"+GENCODE_TARBALL
GENCODE_PARAM_DATABASE_FOLDER_LINK = DATABASES+"/gencode/"+GENCODE_VERSION+"/"

logfile = GENCODE_PARAM_DATABASE_FOLDER_LINK + serviceName + "." +date_time+".database.setup.log"
errfile = GENCODE_PARAM_DATABASE_FOLDER_LINK + serviceName + "." +date_time+".database.setup.err"

# Check if a directory exist
# os.path.isdir('folder') will return true if exist
if os.path.isdir(GENCODE_PARAM_DATABASE_FOLDER_LINK) == False:
	os.makedirs(GENCODE_PARAM_DATABASE_FOLDER_LINK, exist_ok = True)
	logsomefile(logfile, 'Gencode version '+GENCODE_VERSION+' installation start:', "\n", items = date_time)
	installdatabase(GENCODE_PARAM_DATABASE_FOLDER_LINK, GENCODE_SOURCE_EXTERNAL, GENCODE_TARBALL, logfile, errfile)
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	logsomefile(logfile, 'Installation end:', "\n", items = date_time_end)


REFGENOME = DATABASES+"/gencode/"+GENCODE_VERSION+"/"+GENCODEREF_TARBALL
GENCODE = DATABASES+"/gencode/"+GENCODE_VERSION+"/"+GENCODE_TARBALL

# recipe for Salmon index to perform count extraction
# GRch38 assembly and gencode.vxx.transcripts.fa.gz can be found on https://www.gencodegenes.org/human/
# release can be change (actually v41 is used, relase is 07.2022)
# for GRch38 primary assembly v42 : https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
# for gencode transcripts : https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz
# grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
# sed -i.bak -e 's/>//g' decoys.txt
# cat gencode.v41.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
# salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode
# eventually use --keepDuplicates for isoform analysis

systemcall("grep "^>" <(gunzip -c "+REFGENOME+") | cut -d " " -f 1 > decoys.txt")
systemcall("sed -i.bak -e 's/>//g' decoys.txt")
systemcall("cat "+GENCODE+" "+REFGENOME+" > gentrome.fa.gz")
systemcall("salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode")


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