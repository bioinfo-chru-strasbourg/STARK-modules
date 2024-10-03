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
configfile: "/app/config/snakefile/fusions_default.yaml"
###################################################

services_folder = f"{config['services']}/{config['moduleName']}/{config['serviceName'].lower()}"
config_folder = f"{config['config']}/{config['moduleName']}/{config['serviceName'].lower()}"
db = config['databases']
date_time = config['DATE_TIME'] if config['DATE_TIME'] else datetime.now().strftime("%Y%m%d-%H%M%S")

# Set up logging
logfile = {config['serviceName']}.{date_time}.parameters.log"
logging.basicConfig(
	filename=logfile,
	level=config['LOG_LEVEL'],
	format='%(asctime)s - %(levelname)s - %(message)s'
)
log_items = [
	('Start of the analysis:', date_time),
	('Database:', config['databases']),
	('serviceName:', config['serviceName']),
	('moduleName:', config['moduleName'])
]
for item in log_items:
	if isinstance(item[1], list):
		logging.info(f"{item[0]}\n{'\n'.join(map(str, item[1]))}")
	else:
		logging.info(f"{item[0]} {item[1]}")

################## Snakemake Rules ##################
rule all:
		input: f"{db}/CTAT_LIB/CTAT_DB_install.success"

rule install_gencode_db:
    output: f"{db}/gencode/{assembly}.v{gencode_version}/gencode_DB_download.success"
    params:
        genome_link=config['GENCODE_GENOME_LINK'].format(GENCODE_VERSION=config['GENCODE_VERSION'], ASSEMBLY=config['ASSEMBLY']),
        transcripts_link=config['GENCODE_TRANSCRIPTS_LINK'].format(GENCODE_VERSION=config['GENCODE_VERSION']),
        readme_link=config['GENCODE_README_LINK'].format(GENCODE_VERSION=config['GENCODE_VERSION']),
        command=config['COMMAND'],
        db_folder=f"{db}/gencode/{config['ASSEMBLY']}.v{gencode_version}"
    shell:
        """
        mkdir -p {params.db_folder}
        {params.command} {params.genome_link} -O {params.db_folder}/genome.fa.gz
        {params.command} {params.transcripts_link} -O {params.db_folder}/transcripts.fa.gz
        {params.command} {params.readme_link} -O {params.db_folder}/README.txt
        touch {output}
        """

rule install_CTAT_DB:
	input: rules.install_gencode_db.output,
	output:	f"{db}/CTAT_LIB/CTAT_DB_install.success"
	params:
		command=config['COMMAND'],
		ctlib=config['CTAT_LIB'],
		ctat_download=config['CTAT_LINK'],
		ctat_filter_pm_download=config['CTAT_FILTER_PM_LINK']

	shell:
		"""
		echo 'Download CTAT library'
		{params.command} {params.ctat_download} -o {params.ctlib}/CTAT_lib.tar.gz
		{params.command} {params.ctat_filter_pm_download} -o {params.ctlib}/AnnotFilterRule.pm

		echo 'Install CTAT database'
		mkdir -p {params.ctlib}
		tar -xzf {params.ctlib}/CTAT_lib.tar.gz -C {params.ctlib} --strip-components=1
		touch {output}
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
	shell(f"cp {logfile} {services_folder}/cli/{logfile}")

onerror:
	shell(f"rm -f {services_folder}/cli/SETUPRunning.txt")
	shell(f"touch {services_folder}/cli/SETUPFailed.txt")