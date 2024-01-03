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
serviceName = "decon"
moduleName = "structuralvariation"

#serviceName = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_NAME')
#moduleName = os.getenv('DOCKER_STARK_MODULE_NAME')

####################
# DATABASE ANNOTSV #
####################

# https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.1.tar.gz
# include GrCH37 & 38


ANNOTSV_VERSION = "3.3.6"
ANNOTSV_TARBALL = f"Annotations_Human_{ANNOTSV_VERSION}.tar.gz"
ANNOTSV_SOURCE_EXTERNAL = f"https://www.lbgi.fr/~geoffroy/Annotations/{ANNOTSV_TARBALL}"
ANNOTSV_PARAM_DATABASE_FOLDER_LINK = f"{DATABASES}/AnnotSV/{ANNOTSV_VERSION}/"

logfile = f"{ANNOTSV_PARAM_DATABASE_FOLDER_LINK}{serviceName}.{date_time}.database.setup.log"
errfile = f"{ANNOTSV_PARAM_DATABASE_FOLDER_LINK}{serviceName}.{date_time}.database.setup.err"

if not os.path.isdir(ANNOTSV_PARAM_DATABASE_FOLDER_LINK):
	os.makedirs(ANNOTSV_PARAM_DATABASE_FOLDER_LINK, exist_ok = True)
	log_file(logfile, 'AnnotSV version '+ANNOTSV_VERSION+' installation start:', "\n", items = date_time)
	installdatabase(ANNOTSV_PARAM_DATABASE_FOLDER_LINK, ANNOTSV_SOURCE_EXTERNAL, ANNOTSV_TARBALL, logfile, errfile)
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	log_file(logfile, 'Installation end:', "\n", items = date_time_end)

#####################
# DATABASE EXOMISER #
#####################

# Mount point for exomiser jar file
# /share/AnotSV/jar/exomiser-rest-prioritiser-12.1.0.jar in /STARK/databases/AnnotSV/$TOOL_VERSION/jar/
# with $TOOL_VERSION = AnnotSV $TOOL_VERSION (ie 3.1)

# GrCH37 only

###############################################

# https://www.lbgi.fr/~geoffroy/Annotations/2109_hg19.tar.gz

TOOL_NAME="hg19"
TOOL_VERSION="2202"
TOOL_TARBALL = f"{TOOL_VERSION}_{TOOL_NAME}.tar.gz"
TOOL_SOURCE_EXTERNAL = f"https://www.lbgi.fr/~geoffroy/Annotations/{TOOL_TARBALL}"
TOOL_PARAM_DATABASE_FOLDER_LINK = f"{DATABASES}/AnnotSV/{ANNOTSV_VERSION}/Annotations_Exomiser/{TOOL_VERSION}/"

logfile = f"{TOOL_PARAM_DATABASE_FOLDER_LINK}{serviceName}.{date_time}.database.setup.log"
errfile = f"{TOOL_PARAM_DATABASE_FOLDER_LINK}{serviceName}.{date_time}.database.setup.err"

if os.path.isdir(TOOL_PARAM_DATABASE_FOLDER_LINK) == False:
	log_file(logfile, 'Exomiser version '+TOOL_VERSION+' installation start:', "\n", items = date_time)
	installdatabase(TOOL_PARAM_DATABASE_FOLDER_LINK, TOOL_SOURCE_EXTERNAL, TOOL_TARBALL, logfile, errfile)

# https://data.monarchinitiative.org/exomiser/data/2109_phenotype.zip

	TOOL_NAME = "phenotype"
	TOOL_TARBALL = TOOL_VERSION+"_"+TOOL_NAME+".zip"
	TOOL_SOURCE_EXTERNAL = "https://data.monarchinitiative.org/exomiser/data/"+TOOL_TARBALL
	installdatabase(TOOL_PARAM_DATABASE_FOLDER_LINK, TOOL_SOURCE_EXTERNAL, TOOL_TARBALL, logfile, errfile)
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	log_file(logfile, 'Installation end:', "\n", items = date_time_end)

##########
# COSMIC #
##########

# AnnotSV needs the "CosmicCompleteCNA.tsv.gz" file from https://cancer.sanger.ac.uk/cosmic/download (Copy Number Variants part, Download Whole file)
# Put the "CosmicCompleteCNA.tsv.gz" file in the corresponding directory: "$ANNOTSV/share/AnnotSV/Annotations_Human/FtIncludedInSV/COSMIC/GRCh37/
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

##################
# COSMIC install #
##################

# TODO Check if install is also GrCH38 ready

# Install COSMIC SV database (file is CosmicCompleteCNA.tsv.gz placed in /STARK/config/structuralvariation/serviceName/setup/COSMIC/ directory)
# /STARK/config/structuralvariation/serviceName/setup
genomeBuild_version = "GRCh37"
COSMIC_source = f"/STARK/config/{moduleName}/{serviceName}/setup/COSMIC/CosmicCompleteCNA.tsv.gz"
COSMIC_install_path = f"{DATABASES}/AnnotSV/{ANNOTSV_VERSION}/Annotations_Human/FtIncludedInSV/COSMIC/{genomeBuild_version}"
if not os.path.exists(COSMIC_install_path) and os.path.exists(COSMIC_source):
	os.makedirs(COSMIC_install_path, exist_ok = True)
	systemcall(f"mv {COSMIC_source} {COSMIC_install_path}")
	# AnnotSV dummy vcf in /app/src/dummy.vcf to process COSMIC database (need rw databases access for that)
	setup_config_path = f"{DATABASES}/AnnotSV/{ANNOTSV_VERSION}"
	systemcall(f"AnnotSV -SVinputFile /app/dummy/dummy.vcf -outputFile {setup_config_path}/AnnotSV.dummyannotation.tsv -genomeBuild {genomeBuild_version} 1> {setup_config_path}/AnnotSV.dummyannotation.log")


##############
# GeneHancer #
##############

# GeneHancer data is under a specific licence that prevent the systematic availability in AnnotSV sources. Users need to request the up-to-date GeneHancer data dedicated to AnnotSV " #(“GeneHancer_<version>_for_annotsv.zip“) by subscribing to the Genecard web site https://www.genecards.org/
# unzip -q GeneHancer_<version>_for_annotsv.zip -d ($ANNOTSV/share/AnnotSV)/Annotations_Human/FtIncludedInSV/RegulatoryElements/

GENEHANCER_version = "v5.9.zip"
GENEHANCER_source = f"/STARK/config/{moduleName}/{serviceName}/setup/GENEHANCER/GeneHancer_hg19_{GENEHANCER_version}"
GENEHANCER_install_path = f"{DATABASES}/AnnotSV/{ANNOTSV_VERSION}/Annotations_Human/FtIncludedInSV/RegulatoryElements/"
if not os.path.exists(GENEHANCER_install_path) and os.path.exists(GENEHANCER_source):
	os.makedirs(GENEHANCER_install_path, exist_ok = True)
	systemcall(f"unzip -q {GENEHANCER_source} -d {GENEHANCER_install_path}")

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