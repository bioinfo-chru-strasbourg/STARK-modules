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
configfile: "/app/config/snakefile/filt3r_default.yaml"


db = os.getenv('DOCKER_STARK_INNER_FOLDER_DATABASES')  # DATABASES = "/STARK/databases"
serviceName = os.getenv('DOCKER_STARK_MODULE_SUBMODULE_NAME')  # serviceName
moduleName = os.getenv('DOCKER_STARK_MODULE_NAME')  # moduleName = "structuralvariation"
services_folder = f"{os.getenv('DOCKER_STARK_INNER_FOLDER_SERVICES')}/{moduleName}/{serviceName}"
config_folder = f"{os.getenv('DOCKER_STARK_INNER_FOLDER_CONFIG')}/{moduleName}/{serviceName}"

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
rule cp:
	output: f"{services_folder}/cli/SETUPComplete.txt"
	shell: " mkdir -p {config_folder}/listener && cp -r /app/config/module/* {config_folder}/listener && cp -r /app/config/snakefile/* {config_folder}/cli && touch {output} " 

onstart:
	shell(f"touch {services_folder}/cli/SETUPRunning.txt")
	with open(logfile, "a+") as f:
		f.write("\nGlobal parameters of the setup for debug only\n")
		json.dump(config, f, ensure_ascii=False, indent=2)
		f.write("\n")

onsuccess:
	shell(f"rm -f {services_folder}/cli/SETUPRunning.txt")
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S") 
	with open(logfile, "a+") as f:
		f.write(f"End of the setup: {date_time_end}\n")

onerror:
	shell(f"rm -f {services_folder}/cli/SETUPRunning.txt")
	shell(f"touch {services_folder}/cli/SETUPFailed.txt")