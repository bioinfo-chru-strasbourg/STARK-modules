##########################################################################
# SETUP Version:			3.0
# Description:				Setup to configure module
##########################################################################

# DEV version 0.1 : 10/11/2021
# INT version 0.1 : 17/03/2022
# Authoring : Thomas LAVAUX

# PROD version 1 : 17/03/2022
# PROD version 2 : 17/10/2023
# PROD version 3 : 06/09/2024
	# switch to aria2c & snakemake

################## Context ##################
#
# This snakemake will setup a cli ie checking and installing proper databases if needed ; copying launcher.py, .conf and .json files into specific directory
####################################

################## Import libraries ##################
import os
import yaml
import logging
from datetime import datetime
import json

################## Configuration ##################
configfile: "/app/config/snakefile/canoes_default.yaml"
###################################################

ANNOTSV_VERSION = config['ANNOTSV_VERSION']
GENEHANCER_VERSION = config['GENEHANCER_VERSION']
EXOMISER_VERSION = config['EXOMISER_VERSION']
assembly = config['ASSEMBLY']
genomeBuild=config['genomeBuild']
serviceName = config['serviceName']
moduleName = config['module']
db = config['databases']
services = config['services']
config = config['config']
services_folder = {services}/{moduleName}/{serviceName}
config_folder = {config}/{moduleName}/{serviceName}

date_time = config['DATE_TIME'] if config['DATE_TIME'] else datetime.now().strftime("%Y%m%d-%H%M%S")

# Set up logging
logfile = f"{services_folder}/cli/{serviceName}.{date_time}.parameters.log"
logging.basicConfig(
	filename=logfile,
	level=config['LOG_LEVEL'],
	format='%(asctime)s - %(levelname)s - %(message)s'
)
log_items = [
	('Start of the analysis:', date_time),
	('Database:', db),
	('serviceName:', serviceName),
	('moduleName:', moduleName)
]
for item in log_items:
	if isinstance(item[1], list):
		logging.info(f"{item[0]}\n{'\n'.join(map(str, item[1]))}")
	else:
		logging.info(f"{item[0]} {item[1]}")

################## Snakemake Rules ##################

if os.path.exists(f"{services_folder}/setup/COSMIC/CosmicCompleteCNA.tsv.gz"):
	run_cosmic = True
else:
	run_cosmic = False

if os.path.exists(f"{services_folder}/setup/GENEHANCER/{GENEHANCER_VERSION}.zip"):
	run_genehancer = True
else:
	run_genehancer = False

rule all:
	input: f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV.dummyannotation.tsv"

rule download_and_extract_annotSV:
	output: f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV_DB_install.success"
	params:
		folder=f"{db}/AnnotSV/{ANNOTSV_VERSION}",
		command=config['COMMAND'],
		download_link=f"https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_{ANNOTSV_VERSION}.tar.gz"
	shell: """
		mkdir -p {params.folder}
		{params.command} {params.download_link}
		tar xzf Annotations_Human_{ANNOTSV_VERSION}.tar.gz -C {params.folder}
		touch {output}
	"""
# https://data.monarchinitiative.org/exomiser/data/2309_hg19.zip
# https://data.monarchinitiative.org/exomiser/data/2309_hg38.zip
# https://github.com/lgmgeo/AnnotSV/blob/d20fd35999902f41b9d60d32ea26d0e4599cde49/share/AnnotSV/jar/exomiser-rest-prioritiser-12.1.0.jar

rule download_exomiser_data:
	input: rules.download_and_extract_annotSV.output
	output: f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV_Exomiser_install.success"
	params:
		folder=f"{db}/AnnotSV/{ANNOTSV_VERSION}/Annotations_Exomiser/{EXOMISER_VERSION}",
		command=config['COMMAND'],
		download_link1=f"https://data.monarchinitiative.org/exomiser/data/{EXOMISER_VERSION}_{assembly}.zip",
		download_link2=f"https://data.monarchinitiative.org/exomiser/data/{EXOMISER_VERSION}_phenotype.zip",
	shell: """
		mkdir -p {params.folder_exomiser}
		{params.command} {params.download_link1}
		{params.command} {params.download_link2}
		unzip -q {EXOMISER_VERSION}_{assembly}.zip -d {params.folder_exomiser}
		unzip -q {EXOMISER_VERSION}_phenotype.zip -d {params.folder_exomiser}
		mkdir -p {params.folder}/jar
		cp -r /app/config/annotsv/exomiser-rest-prioritiser-12.1.0.jar {params.folder}/jar
		touch {output}
	"""

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
if run_cosmic:
	rule cosmic:
		input: f"{services_folder}/setup/COSMIC/CosmicCompleteCNA.tsv.gz"
		output: f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV_COSMIC_install.success"
		params: f"{db}/AnnotSV/{ANNOTSV_VERSION}/Annotations_Human/FtIncludedInSV/COSMIC/{genomeBuild}/"
		shell: " mkdir -p {params} && unzip -q {input} -d {params} && touch {output} "

# GeneHancer data is under a specific licence that prevent the systematic availability in AnnotSV sources. Users need to request the up-to-date GeneHancer data dedicated to AnnotSV " #(“GeneHancer_<version>_for_annotsv.zip“) by subscribing to the Genecard web site https://www.genecards.org/
# unzip -q GeneHancer_<version>_for_annotsv.zip -d ($ANNOTSV/share/AnnotSV)/Annotations_Human/FtIncludedInSV/RegulatoryElements/
if run_genehancer:
	rule genehancer:
		input: f"{services_folder}/setup/GENEHANCER/{GENEHANCER_VERSION}.zip"
		output: f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV_GENEHANCER_install.success"
		params: f"{db}/AnnotSV/{ANNOTSV_VERSION}/Annotations_Human/FtIncludedInSV/RegulatoryElements/"
		shell: " mkdir -p {params} && unzip -q {input} -d {params} && touch {output} "

rule annotSV_dummy_vcf:
	input:
		success=rules.download_exomiser_data.output,
		file="/app/scripts/dummy/dummy.vcf"
	output: f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV.dummyannotation.tsv"
	params:
		annotsvdir=f"{db}/AnnotSV/{ANNOTSV_VERSION}/",
		genomeBuild=config['genomeBuild']

	log: f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV.dummyannotation.log"
	shell: """
		AnnotSV -SVinputFile {input.file} -annotationsDir {params.annotsvdir} -outputFile {output} -genomeBuild {params.genomeBuild} 1>{log}
	"""

onstart:
	shell(f"touch {services_folder}/cli/SETUPRunning.txt")
	with open(logfile, "a+") as f:
		f.write("\nGlobal parameters of the setup for debug only\n")
		json.dump(config, f, ensure_ascii=False, indent=2)
		f.write("\n")

onsuccess:
	shell(f"rm -f {services_folder}/cli/SETUPRunning.txt")
	shell(f"touch {services_folder}/cli/SETUPComplete.txt")
	shell(f"mkdir -p {config_folder}/listener && cp -r /app/config/module/* {config_folder}/listener && cp -r /app/config/snakefile/* {config_folder}/cli")
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	with open(logfile, "a+") as f:
		f.write(f"End of the setup: {date_time_end}\n")

onerror:
	shell(f"rm -f {services_folder}/cli/SETUPRunning.txt")
	shell(f"touch {services_folder}/cli/SETUPFailed.txt")