##########################################################################
# SETUP Version:			0.1
# Description:				Setup to configure module
##########################################################################

# DEV version 0.1 : 10/11/2021
# INT version 0.1 : 17/03/2022
# Authoring : Thomas LAVAUX

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
	systemcall("wget "+ source+" 1>> "+logfile+" 2>> "+errfile+" ")
	if archive_name.endswith('.zip'):
		systemcall("unzip -q "+archive_name+" -d " +destination+" 1>> "+logfile+" 2>> "+errfile+" ")
	if archive_name.endswith('.tar.gz'):
		systemcall("tar xzf "+archive_name+" -C "+destination+" 1>> "+logfile+" 2>> "+errfile+" ")


# Variables initialisation
# set datetime to add to output file name
date_time = datetime.now().strftime("%Y%m%d-%H%M%S")

### INSTALL DATABASES ###

######################
# DATABASE VARIABLES #
######################

STARK_FOLDER = "/STARK"
DATABASES = STARK_FOLDER+ "/databases"

serviceName = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_NAME')
moduleName = os.getenv('DOCKER_STARK_MODULE_NAME')

####################
# DATABASE ANNOTSV #
####################

# https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.1.tar.gz

ANNOTSV_NAME = "AnnotSV"
ANNOTSV_VERSION = "3.1"
ANNOTSV_TARBALL = "Annotations_Human_"+ANNOTSV_VERSION+".tar.gz"
ANNOTSV_SOURCE_EXTERNAL = "https://www.lbgi.fr/~geoffroy/Annotations/"+ANNOTSV_TARBALL
ANNOTSV_PARAM_DATABASE_FOLDER_LINK = DATABASES+"/AnnotSV/"+ANNOTSV_VERSION+"/"

logfile = ANNOTSV_PARAM_DATABASE_FOLDER_LINK + serviceName + "." +date_time+".database.setup.log"
errfile = ANNOTSV_PARAM_DATABASE_FOLDER_LINK + serviceName + "." +date_time+".database.setup.err"

# export http_proxy=
# export https_proxy=
# export ftp_proxy=

systemcall("export https_proxy=http://hux144_proxy:proxysurf@cyclope:8080 ")

# Check if a directory exist
# os.path.isdir('folder') will return true if exist
if os.path.isdir(ANNOTSV_PARAM_DATABASE_FOLDER_LINK) == False:
	os.makedirs(ANNOTSV_PARAM_DATABASE_FOLDER_LINK, exist_ok = True)
	logsomefile(logfile, 'AnnotSV version '+ANNOTSV_VERSION+' installation start:', "\n", items = date_time)
	installdatabase(ANNOTSV_PARAM_DATABASE_FOLDER_LINK, ANNOTSV_SOURCE_EXTERNAL, ANNOTSV_TARBALL, logfile, errfile)

#####################
# DATABASE EXOMISER #
#####################

# Mount point for exomiser jar file
# /share/AnotSV/jar/exomiser-rest-prioritiser-12.1.0.jar in /STARK/databases/AnnotSV/$TOOL_VERSION/jar/
# with $TOOL_VERSION = AnnotSV $TOOL_VERSION (ie 3.1)

###############################################

# https://www.lbgi.fr/~geoffroy/Annotations/2109_hg19.tar.gz

TOOL_NAME="hg19"
TOOL_VERSION="2109"
TOOL_TARBALL= TOOL_VERSION+"_"+TOOL_NAME+".tar.gz"
TOOL_SOURCE_EXTERNAL="https://www.lbgi.fr/~geoffroy/Annotations/"+TOOL_TARBALL
TOOL_PARAM_DATABASE_FOLDER_LINK= DATABASES+"/AnnotSV/"+ANNOTSV_VERSION+"/Annotations_Exomiser/"+ TOOL_VERSION+"/"

logfile = TOOL_PARAM_DATABASE_FOLDER_LINK + serviceName + "." +date_time+".database.setup.log"
errfile = TOOL_PARAM_DATABASE_FOLDER_LINK + serviceName + "." +date_time+".database.setup.err"

if os.path.isdir(TOOL_PARAM_DATABASE_FOLDER_LINK) == False:
	installdatabase(TOOL_PARAM_DATABASE_FOLDER_LINK, TOOL_SOURCE_EXTERNAL, TOOL_TARBALL, logfile, errfile)
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	logsomefile(logfile, 'Setup end:', "\n", items = date_time_end)

# https://data.monarchinitiative.org/exomiser/data/2109_phenotype.zip

	TOOL_NAME = "phenotype"
	TOOL_TARBALL = TOOL_VERSION+"_"+TOOL_NAME+".zip"
	TOOL_SOURCE_EXTERNAL = "https://data.monarchinitiative.org/exomiser/data/"+TOOL_TARBALL
	installdatabase(TOOL_PARAM_DATABASE_FOLDER_LINK, TOOL_SOURCE_EXTERNAL, TOOL_TARBALL, logfile, errfile)
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	logsomefile(logfile, 'Setup end:', "\n", items = date_time_end)


##########
# COSMIC #
##########

# AnnotSV needs the “CosmicCompleteCNA.tsv.gz” file from https://cancer.sanger.ac.uk/cosmic/download (Copy Number Variants part, Download Whole file)
# Put the “CosmicCompleteCNA.tsv.gz” file in the corresponding directory: “$ANNOTSV/share/AnnotSV/Annotations_Human/FtIncludedInSV/COSMIC/GRCh37/
# create the directory is necessary (./Annotations_Human/FtIncludedInSV/COSMIC/GRCh37/)
# These files will be reprocessed and then removed the first time AnnotSV is executed.

## Auto install COSMIC Annotation for SV
# 1 - Generate an authentication string
# Your first request needs to supply your registered email address and COSMIC password. 
# We use HTTP Basic Auth to check your credentials, which requires you to combine your email address and password and then Base64 encode them. 
# For example, using standard Unix command line tools: 
	# echo "email@example.com:mycosmicpassword" | base64
	# yougetthekey64
# You can use the same authentication string for all of your downloads. You only need to re-generate the string if you change your COSMIC password. 

# 2 - Get a download link
#Make a request to https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v95/CosmicCompleteCNA.tsv.gz. 
# You need to pass the authentication string to the server when you make your request. For example: 
	# curl -H "Authorization: Basic yougetthekey64" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v95/CosmicCompleteCNA.tsv.gz
# That request will return a snippet of JSON containing the link that you need to use to download your file. For example: 
	#{"url" : "https://cog.sanger.ac.uk/cosmic/GRCh37/cosmic/v95/CosmicCompleteCNA.tsv.gz?AWSAccessKeyId=yourID&Expires=yourData" }

# 3 - Download the data file
# Extract the URL from the JSON response and make another request to that URL to download the file. For example:
	# curl "https://cog.sanger.ac.uk/cosmic/GRCh37/cosmic/v95/CosmicCompleteCNA.tsv.gz?AWSAccessKeyId=yourID&Expires=yourData"
# You do not need to supply the authentication header on this request. Download links are valid for one hour. 


# Install COSMIC SV database (file is CosmicCompleteCNA.tsv.gz placed in /STARK/config/structuralvariation/serviceName/setup/COSMIC/ directory)
COSMIC_source = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_SERVICE_SETUP_INNER_FOLDER_CONFIG')+"/setup/COSMIC/CosmicCompleteCNA.tsv.gz" # /STARK/config/structuralvariation/serviceName/setup
COSMIC_install_path = "/STARK/databases/AnnotSV/3.1/Annotations_Human/FtIncludedInSV/COSMIC/GRCh37"
if not os.path.exists(COSMIC_install_path) and os.path.exists(COSMIC_source):
	os.makedirs(COSMIC_install_path, exist_ok = True)
	systemcall("cp "+COSMIC_source+" "+COSMIC_install_path+" ")
# AnnotSV dummy vcf in /app/src/dummy.vcf to process COSMIC database (need rw database for that)
	setup_config_path = "/STARK/databases/AnnotSV/3.1/Annotations_Human/"
	bedtools_path = "/tools/bedtools/current/bin/bedtools"
	genomeBuild_version = "GRCh37"
	systemcall("AnnotSV -SVinputFile /app/src/dummy.vcf -outputFile "+setup_config_path+"/dummy.tsv -bedtools "+bedtools_path+" -genomeBuild "+genomeBuild_version+" ")


##############
# GeneHancer #
##############

# GeneHancer data is under a specific licence that prevent the systematic availability in AnnotSV sources. Users need to request the up-to-date GeneHancer data dedicated to AnnotSV " #(“GeneHancer_<version>_for_annotsv.zip“) by contacting directly the GeneCards team:
# Academic users: genecards@weizmann.ac.il

# unzip -q GeneHancer_<version>_for_annotsv.zip -d ($ANNOTSV/share/AnnotSV)/Annotations_Human/FtIncludedInSV/RegulatoryElements/


#######################
# Copy config files (service.conf & service.json) and launcher.py from app/config to /config/module/servicename/listener/
#######################

if serviceName and moduleName:
	logfile = "/STARK/config/"+moduleName+"/"+serviceName+"/listener/" + serviceName + "." +date_time+".setup.log"
	errfile = "/STARK/config/"+moduleName+"/"+serviceName+"/listener/" + serviceName + "." +date_time+".setup.err"
	logsomefile(logfile, 'Setup start:', "\n", items = date_time)
	systemcall("\cp -r /app/config/* /STARK/config/"+moduleName+"/"+serviceName+"/listener/ 1>> "+logfile+" 2>> "+errfile+" ")
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	logsomefile(logfile, 'Setup end:', "\n", items = date_time_end)

# SETUPComplete cli services (condition for healthy cli)
systemcall("touch /STARK/config/${DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_INNER_FOLDER_SERVICES}/SETUPComplete.txt")