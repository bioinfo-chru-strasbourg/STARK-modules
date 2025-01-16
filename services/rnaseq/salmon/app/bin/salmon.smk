##########################################################################
# Snakemakefile Version:   1.0
# Description:             Snakemake file to run Salmon module
##########################################################################

################## Context ##################
# launch snakemake -s  snakefile_count -c(numberofthreads) --config run=absolutepathoftherundirectory without / at the end of the path
# to launch the snakemake file, use --config to replace variables that must be properly set for the pipeline to work ie run path directory
# every variable defined in the yaml file can be change
# separate multiple variable with a space (ex  --config run=runname var1=0.05 var2=12)
# also use option --configfile another.yaml to replace and merge existing config.yaml file variables

# use -p to display shell commands
# use --lt to display docstrings of rules

# input file = fastq files
# output file = salmon count files (quant.sf)
################## Import libraries ##################

########## Note ########################################################################################
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
########################################################################################################

import os
import glob
import pandas as pd
import json
import csv
import logging
from shutil import copy2
from datetime import datetime
from itertools import product
from collections import defaultdict
from jinja2 import Environment, FileSystemLoader

################## Configuration file ##################
configfile: "/app/config/snakefile/salmon_default.yaml"

####################### FUNCTIONS #####################
def parse_samplesheet(samplesheet_path):
	"""
	samplesheet_path: absolute path of a samplesheet file, Illumina format
	return: a dataframe containing 9 columns :
	Sample_ID, Sample_Plate, Sample_Well, I7_Index_ID, index, Manifest, GenomeFolder, Sample_Project, Description
	The description field contains tags separated by ! ; the name of the tag and the value is separated by # (ex: SEX#F!APP#DIAG.BBS_RP)
	"""
	header_line = next(line.strip().split(',') for line in open(samplesheet_path) if 'Sample_ID' in line)
	df = pd.read_csv(samplesheet_path, skiprows=1, names=header_line)
	df['Description'] = df['Description'].apply(lambda x: dict(item.split('#') for item in x.split('!')))
	
	return df

def getSampleInfos(samplesheet_path, exclude_samples):
	"""
	samplesheet_path: absolute path of a samplesheet file, Illumina format
	return a dictionary with Sample_ID from samplesheet as key and 'gender': 'F or M or NULL'
	"""
	result_dict = {}
	
	if samplesheet_path:
		samplesheet = parse_samplesheet(samplesheet_path)
		
		for _, rows in samplesheet.iterrows():
			sampleID = rows["Sample_ID"]
			
			if any(exclude in sampleID for exclude in exclude_samples):
				continue
			
			result_dict[sampleID] = {'gender': next((tag.split('_')[-1] for tag in rows['Description'].split('!') if 'SEX' in tag), '')}
	
	return result_dict

def populate_dictionary(samples_list, extensions_list, files_list, pattern_include=None, pattern_exclude=None, split_index=0):
	
	dictionary = defaultdict(dict)

	for sample in samples_list:
		for ext in extensions_list:
			found = False
			for file in files_list:
				file_parts = os.path.basename(file).split(".")
				
				if split_index >= len(file_parts):
					continue  # Skip if split_index is out of range
				
				file_base = file_parts[split_index]
				
				if file_base != sample or not os.path.basename(file).endswith(ext):
					continue

				if pattern_exclude and any(exclude in file for exclude in pattern_exclude):
					continue

				if pattern_include and pattern_include not in file:
					continue

				dictionary.setdefault(sample, {})[ext] = file
				found = True
				break  # Stop after the first match

			if not found:
				dictionary.setdefault(sample, {})[ext] = None
	
	return dictionary


def filter_files(files_list, filter_in=None, filter_out=None, extensions=None):
	filtered_files = []

	for file_name in files_list:
		if extensions and any(file_name.endswith(ext) for ext in extensions):
			if filter_out and filter_out in file_name:
				continue  # Skip files with both the specified extension and substring
			if filter_in and filter_in not in file_name:
				continue  # Skip files that do not contain the specified substring
		filtered_files.append(file_name)

	return filtered_files

			
def find_item_in_dict(sample_list, ext_list, dictionary, include_ext, exclude_ext=None):
	""" Function to search in a dictionary for a non-empty file path by iterating through sample_list and ext_list with inclusion and exclusion filters """
	search_result = ""
	
	for sample in sample_list:
		for ext in ext_list:
			try:
				items = dictionary.get(sample, {}).get(ext, [])
				if items is not None:
					if include_ext in items and (exclude_ext is None or exclude_ext not in items):
						if os.path.exists(items) and os.path.getsize(items) != 0:
							search_result = items
			except KeyError as e:
				print(f"KeyError encountered: {e}")
	
	return search_result

def searchfiles(directory, search_args, recursive_arg):
	""" Function to search all files in a directory with multiple search patterns and an option for recursive search """
	results = []
	for search_arg in search_args:
		results.extend(filter(os.path.isfile, glob.glob(directory + search_arg, recursive=recursive_arg)))
	return sorted(results)

def extractlistfromfiles(file_list, ext_list, sep, position):
	""" Function for creating list from a file list, with a specific extension, a separator and the position of the string we want to extract """
	return list(set(os.path.basename(files).split(sep)[position] for files in file_list if any(files.endswith(ext) for ext in ext_list)))

def replace_path(file_paths, old_substring, new_substring):
	return [path.replace(old_substring, new_substring).lstrip("/") for path in file_paths]

def generate_html_report(result_dict, run_name, service_name, sample_list, template_name, output_file='report.html'):
	env = Environment(loader=FileSystemLoader(config['TEMPLATE_DIR']))
	template = env.get_template(template_name)

	rendered_html = template.render(
		runDict=result_dict,
		runName=run_name,
		serviceName=service_name,
		sample_list=sample_list
	)

	with open(output_file, 'w') as f:
		f.write(rendered_html)

	print(f"HTML report generated successfully: {output_file}")

def update_results(dictionary, update_dictionary, keys, exclude_samples=None, exclude_keys=None):
	"""
	Update a dictionary and remove specified keys for certain samples.

	Parameters:
	- dictionary: Dictionary containing sample results.
	- update_dictionary: Dictionary to be updated with results.
	- keys: A list of keys to update.
	- exclude_samples: List of specific samples for which to remove keys.
	- exclude_keys: List of specific keys to remove for certain samples.
	"""
	if exclude_samples is None:
		exclude_samples = []
	if exclude_keys is None:
		exclude_keys = []

	# First, update the dictionary with results
	for sample, results in dictionary.items():
		if sample in update_dictionary:
			for key, value in results.items():
				if key in keys:
					if key not in update_dictionary[sample]:
						update_dictionary[sample][key] = value
					elif update_dictionary[sample][key] is None:
						update_dictionary[sample][key] = value
	
	# Then, remove the specified keys for samples in exclude_samples
	for sample in exclude_samples:
		if sample in update_dictionary:
			for key in exclude_keys:
				update_dictionary[sample].pop(key, None)  # Remove the key if it exists

### END OF FUNCTIONS ###
serviceName = config['serviceName']
runName = os.path.basename(os.path.normpath(config['run']))
date_time = config['DATE_TIME'] if config['DATE_TIME'] else datetime.now().strftime("%Y%m%d-%H%M%S")
resultDir = f"/app/res/{runName}/{date_time}"
outputDir = config['OUTPUT_DIR'] if config['OUTPUT_DIR'] else config['run']

depotDir = config['DEPOT_DIR']
if config['DEPOT_DIR'] == "depository":
	depotDir = outputDir.replace("repository", "depository")

directories = [resultDir, outputDir]
if depotDir:  # Ensure depotDir is not an empty string
	directories.append(depotDir)
for directory in directories:
	os.makedirs(directory, exist_ok=True)

# Search files in repository 
files_list = searchfiles(os.path.normpath(config['run']), config['SEARCH_ARGUMENT'],  config['RECURSIVE_SEARCH'])
# Create sample  list
sample_list = extractlistfromfiles(files_list, config['PROCESS_FILE'], '.', 0)
# Exclude samples from the exclude_list , case insensitive
sample_list = [sample for sample in sample_list if not any(sample.upper().startswith(exclude.upper()) for exclude in config['EXCLUDE_SAMPLE'])]
# If filter_sample_list variable is not empty, it will force the sample list
if config['FILTER_SAMPLE']:
	sample_list = list(config['FILTER_SAMPLE'])
runDict = populate_dictionary(sample_list, config['EXT_INDEX_LIST'], filtered_files, None, None)
print(dict(runDict))

# Log
logfile = f"{resultDir}/{serviceName}.{date_time}.parameters.log"
logging.basicConfig(filename=logfile, level=config['LOG_LEVEL'], format='%(asctime)s %(message)s')
log_items = [
	('Start of the analysis:', date_time),
	('Analysing run:', runName),
	('List of all samples:', sample_list)
]

for item in log_items:
	if isinstance(item[1], list):
		logging.info(f"{item[0]}\n{'\n'.join(map(str, item[1]))}")
	else:
		logging.info(f"{item[0]} {item[1]}")

################################################## RULES ##################################################
ruleorder: copy_fastq > bcl_convert

rule all:
	input:
		expand(f"{resultDir}/{{sample}}/{{sample}}.quant.sf",sample=sample_list)


rule help:
	"""
	General help for salmon module
	Launch snakemake -s salmon.smk -c(numberofthreads) --config DATA_DIR=absolutepathoftherundirectory
	Separate multiple variable with a space (ex --config DATA_DIR=runname transProb=0.05 var1=0.05 var2=12)
	Use -p to display shell commands, use --lt to display docstrings of rules, 
	Input file = fastq PE (SE don't work for now)
	Output file = quant.sf files
	"""

rule bcl_convert:
	""" Convert BCL files to FASTQ using bcl-convert """
	input:
		bcl_dir = config['BCL_DIRECTORY']
	output:
		fastqR1 = f"{resultDir}/{{sample}}.R1.fastq.gz",
		fastqR2 = f"{resultDir}/{{sample}}.R2.fastq.gz"
	params:
		sample_sheet = config['SAMPLE_SHEET'], 
		output_dir = resultDir,  
		sample_name = lambda wildcards: wildcards.sample
		#config_file = config['BCL_CONVERT_CONFIG']  # Path to bcl-convert config JSON file
	threads: workflow.cores
	shell:
		"""
		bcl-convert --input-dir {input.bcl_dir} \
					--output-dir {params.output_dir} \
					--sample-sheet {params.sample_sheet} \
					--threads {threads}

		R1_fastq=$(find {params.output_dir} -type f -name "{params.sample_name}*_R1_*.fastq.gz" | head -n 1)
		R2_fastq=$(find {params.output_dir} -type f -name "{params.sample_name}*_R2_*.fastq.gz" | head -n 1)
		
		mv $R1_fastq {output.fastqR1}
		mv $R2_fastq {output.fastqR2}
		"""


rule copy_fastq:
	output:
		fastqR1=temp(f"{resultDir}/{{sample}}.R1.fastq.gz"),
		fastqR2=temp(f"{resultDir}/{{sample}}.R2.fastq.gz")
	params:
		process = config['PROCESS_CMD'],
		download_link1 = lambda wildcards: runDict[wildcards.sample]['.R1.fastq.gz'],
		download_link2 = lambda wildcards: runDict[wildcards.sample].get('.R2.fastq.gz', None)  # Use .get() to avoid KeyError
	shell:
		"""
		[ \"{params.process}\" = \"ln\" ] && \
		ln -sfn {params.download_link1} {output.fastqR1} && \
		{output.fastqR2} && \
		if [ -f {params.download_link2} ]; then
			ln -sfn {params.download_link2} {output.fastqR2} || rsync -azvh {params.download_link2} {output.fastqR2};
		else
			echo "Warning: {params.download_link2} does not exist, proceeding with single-end data.";
		fi
		|| rsync -azvh {params.download_link1} {output.fastqR1}
		"""


# https://github.com/OpenGene/fastp
rule fastp:
	input:
		fastqR1=f"{resultDir}/{{sample}}.R1.fastq.gz",
		fastqR2=lambda wildcards: f"{resultDir}/{wildcards.sample}.R2.fastq.gz" if runDict[wildcards.sample].get('.R2.fastq.gz') else None  # Make fastqR2 optional
	output:
		fastqR1=temp(f"{resultDir}/{{sample}}.fastp.R1.fastq.gz"),
		fastqR2=temp(f"{resultDir}/{{sample}}.fastp.R2.fastq.gz") if runDict[wildcards.sample].get('.R2.fastq.gz') else None
	params:
		miscs=config['FASTP_GLOBALS_PARAMS'],
		compression=config['FASTP_COMPRESSION'],
		trim=config['FASTP_TRIM']
	threads: workflow.cores
	shell:	
		"""
		if [ -f "{input.fastqR2}" ]; then
			fastp --thread={threads} {params.miscs} {params.trim} --compression={params.compression} --html={wildcards.sample}.QC.html --report_title={wildcards.sample} --in1={input.fastqR1} --in2={input.fastqR2} --out1={output.fastqR1} --out2={output.fastqR2}
		else
			fastp --thread={threads} {params.miscs} {params.trim} --compression={params.compression} --html={wildcards.sample}.QC.html --report_title={wildcards.sample} --in1={input.fastqR1} --out1={output.fastqR1}
		fi
		"""

rule salmon:
	input:
		fastqR1=rules.fastp.output.fastqR1,
		fastqR2=rules.fastp.output.fastqR2  # Ensure fastqR2 is still defined, but may not exist
	output:
		f"{resultDir}/tmp/{{sample}}.salmon.quant/quant.sf"
	params:
		fastquantref = config['REFSALMONINDEX_PATH'],
		misc_options = config['MISC_SALMON_OPTIONS'],
		salmondir= f"{resultDir}/tmp/{{sample}}.salmon.quant/"
	threads: workflow.cores
	shell:
		"""
		if [ -f {input.fastqR2} ]; then
			salmon quant -p {threads} -i {params.fastquantref} -l A -1 {input.fastqR1} -2 {input.fastqR2} {params.misc_options} -o {params.salmondir};
		else
			echo "Warning: {input.fastqR2} does not exist, proceeding with single-end quantification.";
			salmon quant -p {threads} -i {params.fastquantref} -l A -r {input.fastqR1} {params.misc_options} -o {params.salmondir};
		fi
		"""

rule rename:
	input: rules.salmon.output
	output: f"{resultDir}/{{sample}}/{{sample}}.quant.sf"
	shell:
		"""
		mv {input} {output}
		"""


onstart:
	shell(f"touch {os.path.join(outputDir, f'{serviceName}Running.txt')}")
	with open(logfile, "a+") as f:
		f.write("\n")
		f.write("Global parameters of the analysis for debug only")
		json.dump(config, f, ensure_ascii=False, indent=2)
		f.write("\n")

onsuccess:
	include = config['INCLUDE_RSYNC']
	shell(f"rm -f {outputDir}/{serviceName}Running.txt")
	shell(f"touch {outputDir}/{serviceName}Complete.txt")
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	with open(logfile, "a+") as f:
		f.write(f"End of the analysis : {date_time_end}\n")
	
	# Generate dictionary for results
	search_args = [arg.format(serviceName=serviceName) if '{serviceName}' in arg else arg for arg in ["/*/{serviceName}/*/*", "/*"]]
	result_files_list = searchfiles(resultDir, search_args, False)	
	replaced_paths = replace_path(result_files_list, resultDir, "")
	resultDict = populate_dictionary(sample_list, config['RESULT_EXT_LIST'], replaced_paths, pattern_include=serviceName, split_index=2)	
	update_results(runDict, resultDict, config['RESULT_EXT_LIST'])
	generate_html_report(resultDict, runName, serviceName, sample_list, f"{serviceName}.template.html" , f"{resultDir}/{serviceName}.{date_time}.report.html")
	copy2(config['TEMPLATE_DIR'] + '/' + serviceName + '.style.css', resultDir)
	shell("rm -rf {resultDir}/tmp/")

	# Clear existing output directories & copy results
	for sample in sample_list:
		shell(f"rm -f {outputDir}/{sample}/{serviceName}/* || true")
	shell("rsync -azvh --include={include} --exclude='*' {resultDir}/ {outputDir}")
	for sample in sample_list:
		shell(f"cp {outputDir}/{sample}/{serviceName}/{sample}_{date_time}_{serviceName}/* {outputDir}/{sample}/{serviceName}/ || true")

	# Optionally, perform DEPOT_DIR copy
	if config['DEPOT_DIR'] and outputDir != depotDir:
		for sample in sample_list:
			shell(f"rm -f {depotDir}/{sample}/{serviceName}/* || true")
		shell("rsync -azvh --include={include} --exclude='*' {resultDir}/ {depotDir}")
		for sample in sample_list:
			shell(f"cp {outputDir}/{sample}/{serviceName}/{sample}_{date_time}_{serviceName}/* {depotDir}/{sample}/{serviceName}/ || true")


onerror:
	include_log = config['INCLUDE_LOG_RSYNC']
	shell(f"touch {outputDir}/{serviceName}Failed.txt")
	shell(f"rm -f {outputDir}/{serviceName}Running.txt")
	shell("rsync -azvh --include={include_log} --exclude='*' {resultDir}/ {outputDir}")