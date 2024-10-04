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
configfile: "/app/config/snakefile/scramble_default.yaml"
###################################################

# Format versions for download
ANNOTSV_VERSION = config['ANNOTSV_VERSION']
GENEHANCER_VERSION = config['GENEHANCER_VERSION']
EXOMISER_VERSION = config['EXOMISER_VERSION']
# Format paths
db = config['databases']
services_folder = f"{config['services']}/{config['moduleName']}/{config['serviceName'].lower()}"
config_folder = f"{config['config']}/{config['moduleName']}/{config['serviceName'].lower()}"
date_time = config['DATE_TIME'] if config['DATE_TIME'] else datetime.now().strftime("%Y%m%d-%H%M%S")

# Set up logging
logfile = f"{config['serviceName']}_SETUP.{date_time}.parameters.log"
logging.basicConfig(
	filename=logfile,
	level=config['LOG_LEVEL'],
	format='%(asctime)s - %(levelname)s - %(message)s'
)
log_items = [
	('Start of the analysis:', date_time),
	('Database:', db),
	('serviceName:', config['serviceName']),
	('moduleName:', config['moduleName'])
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


rule install_db:
	output:
		annotSV_success=f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV_DB_install.success",
		exomiser_success=f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV_Exomiser_install.success"
	params:
		folder_annotSV= f"{db}/AnnotSV/{ANNOTSV_VERSION}",
		folder_exomiser= f"{db}/AnnotSV/{ANNOTSV_VERSION}/Annotations_Exomiser/{EXOMISER_VERSION}",
		command=config['COMMAND'],
		assembly=config['ASSEMBLY'],
		annotSV_link= config['ANNOTSV_DOWNLOAD_LINK'].format(ANNOTSV_VERSION=ANNOTSV_VERSION),
		exomiser_link1= config['EXOMISER_DOWNLOAD_LINK1'].format(EXOMISER_VERSION=EXOMISER_VERSION, assembly=config['ASSEMBLY']),
		exomiser_link2= config['EXOMISER_DOWNLOAD_LINK2'].format(EXOMISER_VERSION=EXOMISER_VERSION)
	shell:
		"""
		echo 'AnnotSV download and extraction'
		mkdir -p {params.folder_annotSV}
		{params.command} {params.annotSV_link}
		tar xzf Annotations_Human_{ANNOTSV_VERSION}.tar.gz -C {params.folder_annotSV}
		touch {output.annotSV_success}

		echo 'Exomiser data download and extraction'
		mkdir -p {params.folder_exomiser}
		{params.command} {params.exomiser_link1}
		{params.command} {params.exomiser_link2}
		unzip -q {EXOMISER_VERSION}_{params.assembly}.zip -d {params.folder_exomiser}
		unzip -q {EXOMISER_VERSION}_phenotype.zip -d {params.folder_exomiser}
		mkdir -p {params.folder_exomiser}/jar
		cp -r /app/config/annotsv/exomiser-rest-prioritiser-12.1.0.jar {params.folder_exomiser}/jar
		touch {output.exomiser_success}
		"""

if run_cosmic:
	rule cosmic:
		input: f"{services_folder}/setup/COSMIC/CosmicCompleteCNA.tsv.gz"
		output: f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV_COSMIC_install.success"
		params: f"{db}/AnnotSV/{ANNOTSV_VERSION}/Annotations_Human/FtIncludedInSV/COSMIC/{genomeBuild}/"
		shell: " mkdir -p {params} && unzip -q {input} -d {params} && touch {output} "

if run_genehancer:
	rule genehancer:
		input: f"{services_folder}/setup/GENEHANCER/{GENEHANCER_VERSION}.zip"
		output: f"{db}/AnnotSV/{ANNOTSV_VERSION}/AnnotSV_GENEHANCER_install.success"
		params: f"{db}/AnnotSV/{ANNOTSV_VERSION}/Annotations_Human/FtIncludedInSV/RegulatoryElements/"
		shell: " mkdir -p {params} && unzip -q {input} -d {params} && touch {output} "

rule annotSV_dummy_vcf:
	input:
		success=rules.install_db.output.exomiser_success,
		file="/app/scripts/dummy/dummy.vcf"
	output: f"{db}/AnnotSV/{{ANNOTSV_VERSION}}/AnnotSV.dummyannotation.tsv"
	params:
		annotsvdir=f"{db}/AnnotSV/{{ANNOTSV_VERSION}}/",
		genomeBuild=config['genomeBuild']
	log: f"{db}/AnnotSV/{{ANNOTSV_VERSION}}/AnnotSV.dummyannotation.log"
	shell: """
		AnnotSV -SVinputFile {input.file} -annotationsDir {params.annotsvdir} -outputFile {output} -genomeBuild {params.genomeBuild} 1>{log}
	"""

onstart:
	shell(f"rm -f {services_folder}/cli/SETUPComplete.txt")
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
	shell(f"cp {logfile} {services_folder}/cli/{logfile}")
	shell(f"cp {logfile} {db}/AnnotSV/{ANNOTSV_VERSION}/{logfile}")

onerror:
	shell(f"rm -f {services_folder}/cli/SETUPRunning.txt")
	shell(f"touch {services_folder}/cli/SETUPFailed.txt")