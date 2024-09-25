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
configfile: "/app/config/snakefile/salmon_default.yaml"
######################################################

gencode_version = config['GENCODE_VERSION']
assembly = config['ASSEMBLY']
species = config['SPECIES']

db = config['databases']
services_folder = f"{config['services']}/{config['moduleName']}/{config['serviceName']}
config_folder = f"{config['config']}/{config['moduleName']}/{config['serviceName']}

date_time = config['DATE_TIME'] if config['DATE_TIME'] else datetime.now().strftime("%Y%m%d-%H%M%S")

# Set up logging
logfile = f"{services_folder}/cli/{config['serviceName']}.{date_time}.parameters.log"
logging.basicConfig(
	filename=logfile,
	level=config['LOG_LEVEL'],
	format='%(asctime)s - %(levelname)s - %(message)s'
)
log_items = [
	('Start of the analysis:', date_time),
	('Database:', db),
	('serviceName:', config['serviceName']),
	('moduleName:', config['module'])
]
for item in log_items:
	if isinstance(item[1], list):
		logging.info(f"{item[0]}\n{'\n'.join(map(str, item[1]))}")
	else:
		logging.info(f"{item[0]} {item[1]}")

################## Snakemake Rules ##################
rule all:
		input:
			success=f"{db}/salmon/{assembly}.v{gencode_version}/salmon_index/salmon_install.success"

# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/_README.TXT

rule download_gencode:
	output: f"/app/gencode_DB_download.success"
	params:
		genome_link=f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{species}/release_{gencode_version}/{assembly}.primary_assembly.genome.fa.gz",
		transcripts_link=f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{species}/release_{gencode_version}/gencode.v{gencode_version}.transcripts.fa.gz",
		readme_link=f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{species}/release_{gencode_version}/_README.TXT",
		command=config['COMMAND'],
		db_folder=f"{db}/salmon/{assembly}.v{gencode_version}"
	shell:
		"""
		mkdir -p {params.db_folder}
		{params.command} {params.genome_link} 
		{params.command} {params.transcripts_link}
		{params.command} {params.readme_link}
		touch {output}
		"""

# recipe for Salmon index to perform count extraction
# GRch38 assembly and gencode.vxx.transcripts.fa.gz can be found on https://www.gencodegenes.org/human/
# release can be change (actually v41 is used, relase is 07.2022)
# for GRch38 primary assembly v42 : https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
# for gencode transcripts : https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.transcripts.fa.gz

rule create_salmon_index:
	input: rules.download_gencode.output
	output:
		txp2gene=f"{db}/salmon/{assembly}.v{gencode_version}/salmon_index/txp2gene.tsv",
		decoy=f"{db}/salmon/{assembly}.v{gencode_version}/salmon_index/decoys.txt",
		gentrome=f"{db}/salmon/{assembly}.v{gencode_version}/salmon_index/gentrome.fa.gz",
		success=f"{db}/salmon/{assembly}.v{gencode_version}/salmon_index/salmon_install.success"
	params:
		index=f"{db}/salmon/{assembly}.v{gencode_version}/salmon_index",
		genome=f"/app/{assembly}.primary_assembly.genome.fa.gz",
		transcripts=f"/app/gencode.v{gencode_version}.transcripts.fa.gz"
	shell:
		"""
		gunzip -c {params.genome} | grep '^>' | cut -d ' ' -f 1 > {output.decoy}
		sed -i.bak -e 's/>//g' {output.decoy}
		cat {params.transcripts} {params.genome} > {output.gentrome}
		salmon index -t {output.gentrome} -d {output.decoy} -p 12 -i {params.index} --gencode
		zgrep '^>' {output.gentrome}| cut -d '|' -f 1,6 --output-delimiter=$'\\t' - | sed 's/>//g; s/gene_symbol://g; s/"//g' > {output.txp2gene}
		touch {output.success}
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