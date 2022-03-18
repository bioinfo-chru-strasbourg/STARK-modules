##########################################################################
# SETUP Version:			0.1
# Description:				Setup to configure module
##########################################################################

# DEV version 0.1 : 10/11/2021
# INT version 0.1 : 17/03/2022
# Authoring : Thomas LAVAUX

################## Context ##################


####################################


########################
###   Python Func   ####
########################

import os
import subprocess
#import os.path
#import glob


def systemcall(command):
	'''
	*passing command to the shell*
	*return list containing stdout lines*
	command - shell command (string)
	'''
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read().decode('utf8').strip().split('\n')


### INSTALL DATABASES ###

######################
# DATABASE VARIABLES #
######################

STARK_FOLDER = "/STARK"
DATABASES = STARK_FOLDER+ "/databases"

####################
# DATABASE ANNOTSV #
####################

# https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.0.9.tar.gz

ANNOTSV_NAME = "AnnotSV"
ANNOTSV_VERSION = "3.1"
ANNOTSV_TARBALL = ANNOTSV_VERSION+".tar.gz"
ANNOTSV_SOURCE_EXTERNAL = "https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_"+ANNOTSV_VERSION
ANNOTSV_PARAM_DATABASE_FOLDER_LINK = DATABASES+"/AnnotSV/"+ANNOTSV_VERSION+"/"

# Check if a directory exist
# os.path.isdir('folder') will return true if exist
if os.path.isdir(ANNOTSV_PARAM_DATABASE_FOLDER_LINK) == False:
	try:
		systemcall("wget "+ ANNOTSV_SOURCE_EXTERNAL+" ")
		systemcall("tar xf Annotations_Human_"+ANNOTSV_TARBALL+" -C "+ANNOTSV_PARAM_DATABASE_FOLDER_LINK+" ")
	except: pass

## OMIM database ln will be set in the dockercompose.yml
# Create symlink for OMIM 
# OMIM is /STARK/databases/OMIMannotations/current/
# ln -s /STARK/databases/OMIMannotations/current/ /STARK/databases/AnnotSV/3.0.9/Annotations_Human/Genes-based/OMIM/
# ln -s $DATABASES/OMIMannotations/current/$TOOL_PARAM_DATABASE_FOLDER_LINK/Annotations_Human/Genes-based/OMIM/

##########
# COSMIC #
##########

# AnnotSV needs the “CosmicCompleteCNA.tsv.gz” file from https://cancer.sanger.ac.uk/cosmic/download (Copy Number Variants part, Download Whole file)
# Put the “CosmicCompleteCNA.tsv.gz” file in the corresponding directory: “$ANNOTSV/share/AnnotSV/Annotations_Human/FtIncludedInSV/COSMIC/GRCh37/
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

##############
# GeneHancer #
##############

# GeneHancer data is under a specific licence that prevent the systematic availability in AnnotSV sources. Users need to request the up-to-date GeneHancer data dedicated to AnnotSV " #(“GeneHancer_<version>_for_annotsv.zip“) by contacting directly the GeneCards team:
# Academic users: genecards@weizmann.ac.il

# unzip -q GeneHancer_<version>_for_annotsv.zip -d $ANNOTSV/share/AnnotSV/Annotations_Human/FtIncludedInSV/RegulatoryElements/  

#####################
# DATABASE EXOMISER #
#####################

# Mount point for exomiser jar file
# /share/AnotSV/jar/exomiser-rest-prioritiser-12.1.0.jar in /STARK/databases/AnnotSV/$TOOL_VERSION/jar/
# with $TOOL_VERSION = AnnotSV $TOOL_VERSION (ie 3.0.9)

###############################################

# https://www.lbgi.fr/~geoffroy/Annotations/2109_hg19.tar.gz

TOOL_NAME="hg19"
TOOL_VERSION="2109"
TOOL_TARBALL= TOOL_VERSION+"_"+TOOL_NAME+".tar.gz"
TOOL_SOURCE_EXTERNAL="https://www.lbgi.fr/~geoffroy/Annotations/"+TOOL_TARBALL
TOOL_PARAM_DATABASE_FOLDER_LINK= DATABASES+"/AnnotSV/"+ANNOTSV_VERSION+"/Annotations_Exomiser/"+ TOOL_VERSION+"/"

if os.path.isdir(TOOL_PARAM_DATABASE_FOLDER_LINK) == False:
	try:
		systemcall("wget "+ TOOL_SOURCE_EXTERNAL+" ")
		systemcall("tar -xf "+TOOL_TARBALL+"-C " +TOOL_PARAM_DATABASE_FOLDER_LINK+" ")
	except: pass

# https://data.monarchinitiative.org/exomiser/data/2109_phenotype.zip

	TOOL_NAME = "phenotype"
	TOOL_TARBALL = TOOL_VERSION+"_"+TOOL_NAME+".zip"
	TOOL_SOURCE_EXTERNAL = "https://data.monarchinitiative.org/exomiser/data/"+TOOL_TARBALL

	try:
		systemcall("wget "+ TOOL_SOURCE_EXTERNAL+" ")
		systemcall("unzip -q "+TOOL_TARBALL+"-d " +TOOL_PARAM_DATABASE_FOLDER_LINK+" ")
	except : pass

#######################
# Copy config files (service.conf & service.json) and launcher.py from app/condig to /config/module/servicename/listener/
#######################

if not serviceName:
	serviceName = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_NAME')
if not moduleName:
	moduleName = os.getenv('DOCKER_STARK_MODULE_NAME')

if serviceName and moduleName:
	systemcall("cp /app/config/* /config/"+moduleName+"/"+serviceName+"/listener/ ")
