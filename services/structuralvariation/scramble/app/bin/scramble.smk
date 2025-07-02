##########################################################################
# Snakemakefile Version:   3.0
# Description:             Snakemake file to run SCRAMBLE module
##########################################################################

# DEV version 0.1 : 10/11/2021
# INT version 0.1 : 17/03/2022
# PROD version 1 : 03/06/2022
# Authoring : Thomas LAVAUX

# PROD version 2 : 14/10/2022 changelog
	# remove panel vcf filtering (output was essentially empty) ; rename full to unfiltered
	# exclude samples list case insensitive
	# copying new analysis per sample in the root sample dir & removing the old ones
	# add vcf2tsv converter (https://github.com/sigven/vcf2tsvpy) & corresponding rules : each vcf will be convert to tsv
	# correct run path by removing ending '/' if exist
	# keep "Running.txt" file if failed, avoiding multiple analysis launch by the listener

# PROD version 3 : 19/09/2023 changelog
	# AnnotSV version 3.4.2 (include the vcf converter) ; add -snvIndelFiles option
	# Add the deletion mode to SCRAMBLE
	# refactor snakemake code a lot ; searchfiles function use a list of folder 
	
################## Import libraries ##################
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
configfile: "/app/config/snakefile/scramble_default.yaml"

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
	""" Function to search all files in a directory list with multiple search patterns and an option for recursive search """
	results = []
	for search_arg in search_args:
		results.extend(filter(os.path.isfile, glob.glob(directory + search_arg, recursive=recursive_arg)))
	return sorted(results)

def extractlistfromfiles(file_list, ext_list, sep, position):
	""" Function for creating list from a file list, with a specific extension, a separator and the position of the string we want to extract """
	return list(set(os.path.basename(files).split(sep)[position] for files in file_list if any(files.endswith(ext) for ext in ext_list)))

def replace_path(file_paths, old_substring, new_substring):
	return [path.replace(old_substring, new_substring).lstrip("/") for path in file_paths]

def extract_tag(tagfile, tag, tagsep, sep):
	""" Function to extract a tag """
	row = open(tagfile, 'r').readline().strip().split(tagsep)
	output_tag = next((items.split(sep)[-1] for items in row if tag in items and sep in items), "")
	return output_tag

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

def pedigree_into_dict(ped_file, sample_list, dictionary, individual_id_column='Individual ID', hpo_list_column='HPOList', sex_column='Sex'):
	"""
	Add dictionary based on data from a pedigree file.

	Parameters:
		ped_file (str): Path to the ped file.
		sample_list (list): List of sample IDs to filter.
		dictionary (dict): Dictionary to update with HPO and sex information.
		individual_id_column (str): Column name for individual IDs.
		hpo_list_column (str): Column name for HPO lists.
		sex_column (str): Column name for sex information.

	Returns:
		None
	"""
	df = pd.read_csv(ped_file, sep='\t', dtype=str)
	for index, row in df.iterrows():
		individual_id = row[individual_id_column]
		hpo_list = row[hpo_list_column]
		gender = row[sex_column]

		if isinstance(hpo_list, float):
			hpo_list = str(hpo_list)
		hpo_list = hpo_list.replace(" ", "")

		# Convert Sex to M/F/None
		if gender == '1':
			gender = 'M'
		elif gender == '2':
			gender = 'F'
		else:
			gender = None

		if individual_id in sample_list:
			dictionary.setdefault(individual_id, {})['hpo'] = hpo_list
			dictionary.setdefault(individual_id, {})['gender'] = gender

def process_gene_list(file, dest_dir, res_list):
	# The truncated panel name is the name of the file without the first part (the sample name) and the .genes extension
	panel_name_trunc = file.split('.', 1)[-1].replace('.genes', '')
	
	# Copy the file to the destination directory with the truncated name
	copy2(file, os.path.join(dest_dir, panel_name_trunc)) 
	
	# Append the truncated panel name to the provided result list
	res_list.append(panel_name_trunc)
	
	# Return the updated list
	return res_list

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
date_time = config['DATE_TIME'] if config['DATE_TIME'] else datetime.now().strftime("%Y%m%d-%H%M%S")
runName = os.path.basename(os.path.normpath(config['run']))
resultDir = f"/app/res/{runName}/{date_time}"
outputDir = config['OUTPUT_DIR'] if config['OUTPUT_DIR'] else config['run']
vcf_extension = config['vcf_extension']

depotDir = config['DEPOT_DIR']
if config['DEPOT_DIR'] == "depository":
	depotDir = outputDir.replace("repository", "depository")

directories = [resultDir, outputDir]
if depotDir:  # Ensure depotDir is not an empty string
	directories.append(depotDir)
for directory in directories:
	os.makedirs(directory, exist_ok=True)

print('[INFO] Starting searching files with the parameters provided')
files_list = searchfiles(os.path.normpath(config['run']), config['SEARCH_ARGUMENT'],  config['RECURSIVE_SEARCH'])
print('[INFO] Searching files done')

# Create sample and aligner list
sample_list = extractlistfromfiles(files_list, config['PROCESS_FILE'], '.', 0)
aligner_list = extractlistfromfiles(files_list, config['PROCESS_FILE'], '.', 1)

# Exclude samples from the exclude_list , case insensitive
sample_list = [sample for sample in sample_list if not any(sample.upper().startswith(exclude.upper()) for exclude in config['EXCLUDE_SAMPLE'])]

# If filter_sample_list variable is not empty, it will force the sample list
if config['FILTER_SAMPLE']:
	print('[INFO] Filtering samples list')
	sample_list = list(config['FILTER_SAMPLE'])
	print('[INFO] Filtering of samples list done')

# For validation analyse bam will be sample.aligner.validation.bam, so we append .validation to all the aligner strings
if config['VALIDATION_ONLY']:
	filtered_files = filter_files(files_list, filter_in='validation', extensions=config['PROCESS_FILE'])
	append_aligner = '.validation'
	aligner_list = [sub + append_aligner for sub in aligner_list]
else:
	filtered_files = filter_files(files_list, None ,filter_out='validation', extensions=config['PROCESS_FILE'])

print('[INFO] Construct the dictionary for the run')
runDict = populate_dictionary(sample_list, config['EXT_INDEX_LIST'], filtered_files, None, ['analysis'])
print('[INFO] Dictionary done')

# Set a filelist with all the files tag ; file format is sample.tag
# tag ex SEX#M!POOL#POOL_HFV72AFX3_M_10#POOL_HFV72AFX3_F_11!
for individual in runDict.values():
	individual.update({'hpo': None, 'gender': None})

print('[INFO] Seaching for pedigree files')
ped_file = find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.ped')

if os.path.exists(ped_file) and os.path.getsize(ped_file) > 0:
	print('[INFO] Reading pedigree files')
	pedigree_into_dict(ped_file, sample_list, runDict, individual_id_column='Individual ID', hpo_list_column='HPOList', sex_column='Sex')
	print('[INFO] Gender found in the pedigree file')
	df = pd.read_csv(ped_file, sep='\t', dtype=str)
	for _, row in df.iterrows():
		individual_id = row['Individual ID']
		hpo_list = str(row['HPOList']).replace(" ", "")
		if individual_id in sample_list:
			runDict[individual_id]['hpo'] = hpo_list
else:
	print('[INFO] No pedigree files found, reading gender from tag files')
	tagfile_list = [runDict[sample].get('.tag') for sample in sample_list if '.tag' in runDict.get(sample, {})]
	tagfile_list = [tag for tag in tagfile_list if tag is not None]
	if not tagfile_list:
		print('[INFO] No tag files found to find the genders of the samples')
		print('[WARN] The calling of chrXX & chrXY will not be done')
	else:
		for tagfile in tagfile_list:
			sample_name = os.path.basename(tagfile).split(".")[0]
			runDict[sample_name]['gender'] = extract_tag(tagfile, 'SEX', '!', '#')
		print('[INFO] Gender found in the tag files')

# Find bed file (Design)
config['BED_FILE'] = config['BED_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.design.bed', '.genes.bed')
# Find genes file (Panel); we can't use .genes files because .list.genes and .genes are not distinctable from the indexing we made
config['GENES_FILE'] = config['GENES_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.genes.bed', '.list.genes')
# Find list.genes files (a list of panel files)
config['LIST_GENES'] = config['LIST_GENES'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.list.genes', '.list.transcripts')
# Find transcripts files (file containing NM)
config['TRANSCRIPTS_FILE'] = config['TRANSCRIPTS_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.transcripts', '.list.transcripts')

# Find transcripts files (NM)
config['TRANSCRIPTS_FILE'] = config['TRANSCRIPTS_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.transcripts', '.list.transcripts')
# If transcript file exist, create the annotation file for AnnotSV
annotation_file = f"{resultDir}/{serviceName}.{date_time}.AnnotSV.txt"
if os.path.exists(config['TRANSCRIPTS_FILE']) and os.path.getsize(config['TRANSCRIPTS_FILE']):
	print('[INFO] Processing transcript file')
	df = pd.read_csv(config['TRANSCRIPTS_FILE'], sep='\t', names=["NM", "Gene"])
	with open(annotation_file, 'w+') as f:
		f.write('\t'.join(df['NM']))
	print('[INFO] Processing transcript file done')
else:
	print('[WARN] No transcript file found')
	with open(annotation_file, 'w') as f:
		f.write("No NM found")

# Transform list_genes into a list if list_genes exist, else use genes_file if exist
panels_list = []
file_source = config['LIST_GENES'] or config['GENES_FILE']
if file_source:
	print('[INFO] Processing panel bed files')
	with open(file_source) as f:
		# files is either a list of files in the LIST_GENES file or a single file == GENES_FILE ; ext is .genes
		files = f.read().splitlines() if config['LIST_GENES'] else [os.path.basename(file_source)]
	for file in files:
		if not os.path.isabs(file): # if the file name not contains the path we append it
			file = os.path.join(os.path.dirname(file_source), file)
		
		# Process the gene list and accumulate results
		res_list = []
		res_list = process_gene_list(file, resultDir, res_list)
		panels_list.extend(res_list)  # Append results to panels_list

	print('[INFO] Processing panel bed files done')

# Log
logfile = f"{resultDir}/{serviceName}.{date_time}.parameters.log"
logging.basicConfig(filename=logfile, level=config['LOG_LEVEL'], format='%(asctime)s %(message)s')

log_items = [
	('Start of the analysis:', date_time),
	('Analysing run:', runName),
	('List of all samples:', sample_list),
	('Aligner list from files:', aligner_list),
	('Design bed file:', config['BED_FILE']),
	('Panel bed file:', config['GENES_FILE']),
	('Transcripts file:', config['TRANSCRIPTS_FILE']),
	('Genes list file', config['LIST_GENES'])
]
for item in log_items:
	if isinstance(item[1], list):
		logging.info(f"{item[0]}\n{'\n'.join(map(str, item[1]))}")
	else:
		logging.info(f"{item[0]} {item[1]}")

print(dict(runDict))

################################################## RULES ##################################################
# Priority order
ruleorder: copy_bam > copy_cram > cramtobam > indexing

# To avoid expansion of aligner into bwamem.AnnotSV for some weird reason we constraint the wildcard
wildcard_constraints:
	aligner="|".join(aligner_list)

rule all:
	"""Rule will create an unfiltered vcf.gz and corresponding tsv, with an optional design/panel vcf.gz & tsv if a bed/gene file is provided"""
	input:
		expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full.vcf.gz", sample=sample_list, aligner=aligner_list) +
		expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Full.vcf.gz", aligner=aligner_list) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design.vcf.gz", aligner=aligner_list)  if config['BED_FILE'] else []) +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Design.vcf.gz", sample=sample_list, aligner=aligner_list) if config['BED_FILE'] else []) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel.{{panel}}.vcf.gz", panel=panels_list, aligner=aligner_list) if config['GENES_FILE'] else []) +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Panel.{{panel}}.vcf.gz", aligner=aligner_list, sample=sample_list, panel=panels_list) if config['GENES_FILE'] else []) +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.vcf.gz", sample=sample_list, aligner=aligner_list) if config['ANNOTSV_ANNOTATION'] else []) +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.tsv", sample=sample_list, aligner=aligner_list) if config['ANNOTSV_ANNOTATION'] else [])  +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.AnnotSV.Design.vcf.gz", aligner=aligner_list) if config['ANNOTSV_ANNOTATION'] else []) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.AnnotSV.Design.tsv", aligner=aligner_list) if config['ANNOTSV_ANNOTATION'] else [])  +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.vcf.gz", aligner=aligner_list, sample=sample_list, panel=panels_list) if config['GENES_FILE'] and config['ANNOTSV_ANNOTATION'] else []) +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv", aligner=aligner_list, sample=sample_list, panel=panels_list) if config['GENES_FILE'] and config['ANNOTSV_ANNOTATION']
		else []) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.AnnotSV.Panel.{{panel}}.vcf.gz", aligner=aligner_list, panel=panels_list) if config['GENES_FILE'] and config['ANNOTSV_ANNOTATION'] else []) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv", aligner=aligner_list, panel=panels_list) if config['GENES_FILE'] and config['ANNOTSV_ANNOTATION'] else [])

rule help:
	"""
	General help for SCRAMBLE module
	Launch snakemake -s scramble.smk -c(numberofthreads) --config run=absolutepathoftherundirectory
	To launch the snakemake file, use --config to replace variables that must be properly set for the pipeline to work ie run path directory
	Every variable defined in the yaml file can be change
	Separate multiple variable with a space (ex  --config DATA_DIR=runname transProb=0.05 var1=0.05 var2=12)
	Use option --configfile another.yaml to replace and merge existing config.yaml file variables
	Use -p to display shell commands
	Use --lt to display docstrings of rules
	Use -n for a dry run
	Input file = cram or bam files (if bai is needed, it will be generate)
	Output file = tsv & vcf annoted with AnnotSV 3.x for each sample/bam, and a global tsv & vcf file will all the samples ; a set of tsv & vcf files by design/panel
	"""

rule copy_bam:
	output: temp(f"{resultDir}/{{sample}}.{{aligner}}.bam")
	params:
		process=config['PROCESS_CMD'],
		download_link=lambda wildcards: runDict[wildcards.sample]['.bam']
	message: """ Copying {params.download_link} file """
	shell: "[ \"{params.process}\" = \"ln\" ] && ln -sfn {params.download_link} {output} || rsync -azvh {params.download_link} {output}"

rule copy_cram:
	output: temp(f"{resultDir}/{{sample}}.{{aligner}}.cram")
	params:
		process=config['PROCESS_CMD'],
		download_link=lambda wildcards: runDict[wildcards.sample]['.cram']
	message: """ Copying {params.download_link} file """
	shell: "[ \"{params.process}\" = \"ln\" ] && ln -sfn {params.download_link} {output} || rsync -azvh {params.download_link} {output}"

rule cramtobam:
	""" Extract bam from a cram file with samtools, need a reference genome """
	input: f"{resultDir}/{{sample}}.{{aligner}}.cram"
	output: temp(f"{resultDir}/{{sample}}.{{aligner}}.bam")
	params: config['REFGENEFA_PATH']
	shell: "samtools view -b -T {params} -o {output} {input}"

rule indexing:
	""" Indexing bam files with samtools if .bam.bai is not in PROCESS_FILE """
	input: f"{resultDir}/{{sample}}.{{aligner}}.bam"
	output: temp(f"{resultDir}/{{sample}}.{{aligner}}.bam.bai")
	params:
		process=config['PROCESS_CMD'],
		download_link=lambda wildcards: runDict[wildcards.sample]['.bam.bai'],
		process_file_str=" ".join(config['PROCESS_FILE']) 
	threads: workflow.cores
	shell: "[[ \"{params.process_file_str}\" =~ \"bam.bai\" ]] && ( [ \"{params.process}\" = \"ln\" ] && ln -sfn {params.download_link} {output} || rsync -azvh {params.download_link} {output} ) || samtools index -b -@ {threads} {input} {output}"

rule cluster_identifier:
	"""
	Cluster identifier will identify soft clipped clusters.
	The output is a tab delimited text file with clipped cluster consensus sequences. The columns are as follows:
	1. Coordinate
	2. Side of read where soft-clipped occurred
	3. Clipped read consensus
	4. Anchored read consensus
	Requirement : .bam file must have a .bai file associated in the same folder
	"""
	input:
		bam=f"{resultDir}/{{sample}}.{{aligner}}.bam",
		bai=f"{resultDir}/{{sample}}.{{aligner}}.bam.bai",
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.clusters.txt")
	params:
		mini = config['m'],
		soft = config['s'],
		region = config['r']
	message: """ Starting extraction of clusters for {input.bam} file """
	shell: "cluster_identifier -m {params.mini} -s {params.soft} -r {params.region} {input.bam} > {output}"

# Indels option : you will need a reference genome indexed with the following command for blast (indels mode) to work : makeblastdb -in refgene.fa -dbtype nucl
rule scramble:
	"""
	Calling SCRAMble.R with --eval-mei produces a tab delimited file. If a genomereference.fa file is provided, then a VCF 1 based is produced as well.
	The <out-name>_MEIs.txt output is a tab delimited text file with MEI calls. If no MEIs are present an output txt file will still be produced with only the header, and a dummy vcf will be output as well.
	Calling SCRAMble.R with --eval-dels produced a VCF and a tab delimted file. The <out-name>_PredictedDeletions.txt output is a tab delimited text file with deletion calls. If no deletions are present an output file will still be produced with only the header.
	"""
	input:
		rules.cluster_identifier.output
	params:
		scrambledir = config['SCRAMBLE_PATH'],
		refmei = config['REFMEI_PATH'],
		nCluster = config['nCluster'],
		meiscore = config['mei-score'],
		polyafrac = config['poly-a-frac'],
		polyadist = config['poly-a-dist'],
		mode = config['scramble_mode'],
		refgene = config['REFGENEFA_PATH'],
		dummypath = config['DUMMY_FILES'],
		output = f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full_raw"
	output:
		temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full_raw.vcf")
	log: 
		log = f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.ScrambleR.log", 
		err = f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.ScrambleR.err"
	message: """ Starting Scramble on {input} file """
	shell:
		"""
		Rscript --vanilla {params.scrambledir}/SCRAMble.R \
		--cluster-file {input} \
		--out-name {params.output} \
		--install-dir {params.scrambledir} \
		--mei-refs {params.refmei} \
		--ref {params.refgene} \
		--mei-score {params.meiscore} --nCluster {params.nCluster} --poly-a-dist {params.polyadist} --poly-a-frac {params.polyafrac} {params.mode} > {log.log} 2> {log.err} && [[ -s {output} ]] || cat {params.dummypath}/empty.vcf > {output}
		"""

rule correct_vcf_scramble:
	"""	Correction of vcf output, add sample name and genotype to be consistent with the vcf format specification """
	input: rules.scramble.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full_corr.vcf")
	params:
		mode = config['scramble_mode'],
		dummypath = config['DUMMY_FILES']
	shell:
		"""
		(grep "^##" {input} && echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' && grep "^#CHROM" {input} | awk -v SAMPLE={wildcards.sample} '{{print $0"\tFORMAT\t"SAMPLE}}' && grep "^#" -v {input} | awk '{{print $0"\tGT\t0/1"}}') > {output} && [[ -s {output} ]] || cat {params.dummypath}/empty.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output}
		"""

rule vcf2gz:
	input: rules.correct_vcf_scramble.output
	output: 
		vcfgz=temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full_unfiltered.vcf.gz"),
		vcfgztbi=temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full_unfiltered.vcf.gz.tbi")
	params: config['DUMMY_FILES']
	shell: "bgzip -c {input} > {output.vcfgz} ; tabix {output.vcfgz}"


rule bcftools_filter:
	"""	Filter with bcftools """
	input: rules.vcf2gz.output.vcfgz	
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full.vcf.gz"
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full.bcftoolsfilter.log"
	params: bed=config['BCFTOOLS_FILTER'],
			dummy=config['DUMMY_FILES']
	shell: "bcftools view {params.bed} {input} -o {output}  2> {log} && [[ -s {output} ]] || cat {params.dummy}/empty.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' | bgzip > {output} ; tabix {output}" 

rule merge_vcf:
	""" Merge multiple vcfs using bcftools """
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full.vcf.gz", sample=sample_list, aligner=aligner_list)
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Full.vcf.gz"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Full.bcftoolsmerge.log"
	shell: "bcftools merge --force-single -W=tbi {input} -O z -o {output} 2> {log} && [[ -s {output} ]] || echo -e '## Dummy file created because there's no variant found in any samples\n## You can check individual samples for confirmation' | gzip > {output}; tabix {output} || true"

# we filter Full to Design
# Design individual samples vcf.gz no annotation
rule filter_vcf:
	"""	Filter vcf with a bed file """
	input: rules.bcftools_filter.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Design.vcf.gz"
	params: config['BED_FILE']
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Design.bedtoolsfilter.log"
	shell: "bedtools intersect -header -a {input} -b {params} 2> {log} | bgzip > {output} ; tabix {output}"

use rule merge_vcf as merge_vcf_design with:
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Design.vcf.gz", sample=sample_list, aligner=aligner_list)
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design.vcf.gz"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design.bcftoolsmerge.log"

# We filter non annoted design to get panels
rule filter_vcf_panel:
	"""	Filter vcf with a bed file """
	input: rules.filter_vcf.output
	output: 
		vcfgz=temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Panel_unnorm.{{panel}}.vcf.gz"),
		vcfgztbi=temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Panel_unnorm.{{panel}}.vcf.gz.tbi")
	params: lambda wildcards: f"{resultDir}/{wildcards.panel}"
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Panel.bedtoolsfilter.{{panel}}.log"
	shell: "bedtools intersect -header -a {input} -b {params} 2> {log} | bgzip > {output.vcfgz} ; tabix {output.vcfgz}"

# Panel vcf.gz individual samples no annotation
rule vcf_normalization:
	input: 
		vcfgz=rules.filter_vcf_panel.output.vcfgz,
		vcfgztbi=rules.filter_vcf_panel.output.vcfgztbi
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Panel.{{panel}}.vcf.gz"
	shell: "bcftools norm  -W=tbi -d all -o {output} -Oz {input.vcfgz}"

# Panel vcf.gz all samples no annotation
use rule merge_vcf as merge_vcf_panel_noannotation with:
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Panel.{{panel}}.vcf.gz", sample=sample_list, aligner=aligner_list, panel=panels_list)
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel.{{panel}}.vcf.gz"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel.{{panel}}.bcftoolsmerge.log"


rule AnnotSV:
	"""
	Annotate and rank Structural Variations from a vcf file
	Output a vcf file, if no annotation is done output an empty vcf
	- annotationMode can be split by exons/introns or full by genes
	- txtFile: path to a file containing a list of preferred genes transcripts for annotation
	- genomeBuild must be specified if not hg38
	"""
		input: rules.filter_vcf.output
		output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.tsv"
		log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.log"
		params:
			dummypath=config['DUMMY_FILES'],
			genome=config['genomeBuild'],
			overlap=config['overlap'],
			mode=config['annotationMode'],
			annotation=config['annotationdir'],
			hpo=lambda wildcards: runDict[wildcards.sample]['hpo'],
			samplevcf=lambda wildcards: runDict[wildcards.sample][vcf_extension]
		resources:
			AnnotSVjobs=1
		shell:
			"""
			AnnotSV -SVinputFile {input} -outputFile {output} -snvIndelFiles {params.samplevcf} -annotationMode {params.mode} \
					-annotationsDir {params.annotation} -hpo {params.hpo} -txFile {annotation_file} \
					-genomeBuild {params.genome} -overlap {params.overlap} > {log} && \
					([[ -s {output} ]] || (cat {params.dummypath}/emptyAnnotSV.tsv | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output}))
			"""

rule wait_for_AnnotSV:
	"""
	AnnotSV is not compatible with multiple jobs running at the same time, so we need to wait for it to finish before moving on to the next sample
	"""
	input:
		output_from_AnnotSV=rules.AnnotSV.output,
		log_file= f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.log"
	output:
		ready=f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.ready"
	params:
		limit=config['annotSV_limit']
	shell:
		"""
		count=0
		limit={params.limit}
		while true; do
			if grep -q "Exit without error" {input.log_file} || grep -q "AnnotSV is done with the analysis" {input.log_file}; then
				touch {output.ready}
				break
			fi

			count=$((count + 1))

			if [ $count -ge $limit ]; then
				echo "[ERROR] AnnotSV failed to annotate {wildcards.sample} after $limit attempts. Check the log file for more information."
				exit 1
			fi
			
			sleep 10  # in sec
		done
		"""

# Design tsv all samples AnnotSV
rule merge_tsv:
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.tsv", sample=sample_list, aligner=aligner_list)
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.AnnotSV.Design.tsv"
	params:	config['MERGE_SCRIPT']
	shell: " python {params} -i {input} -o {output} "


rule correct_tsv:
	input: 
		AnnotSV=rules.AnnotSV.output,
		ready=rules.wait_for_AnnotSV.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_fix.tsv")
	shell: 
		"""
		sed '1s/[()/\\]/_/g' {input.AnnotSV} > {output}
		"""

# Design tsv individual samples AnnotSV
rule correct_chr:
	"""
	Add chr to the SV_chrom column
	"""
	input: rules.correct_tsv.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_chr.tsv")
	shell:
		"""
		awk '{{{{FS=OFS="\t"}};if(NR==1){{print; next}}; $2="chr"$2; print}}' {input} > {output}
		"""

rule variantconvert:
	input: rules.correct_chr.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_unfix.vcf")
	params:
			variantconvertannotsv=config['VARIANT_CONVERT_ANNOTSV'],
			dummypath=config['DUMMY_FILES']
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.variantconvert.log"
	shell:
		"""
		variantconvert convert -i {input} -o {output} -c {params.variantconvertannotsv} 2> {log} && [[ -s {output} ]] || cat {params.dummypath}/emptyAnnotSV.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output}
		"""

rule correct_vcf:
	input: rules.variantconvert.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_fix.vcf")
	shell:
		""" 
		sed 's/SAMPLENAME/{wildcards.sample}/g' {input} > {output}
		"""

rule sort_vcf:
	input: rules.correct_vcf.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_sort.vcf")
	params: config['DUMMY_FILES']
	shell: " {{ grep \'^#\' {input} && grep -v \'^#\' {input} | sort -k1,1V -k2,2g; }} > {output} && [[ -s {output} ]] || cat {params}/empty.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output} "

rule fix_vcf:
	input: rules.sort_vcf.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.vcf.gz"
	params:	config['MULTIFIX_SCRIPT']
	shell: "python {params} -i {input} -o {output} -z ; tabix {output}"

# Design vcf.gz all samples AnnotSV
use rule merge_vcf as merge_vcf_annotation with:
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.vcf.gz", sample=sample_list, aligner=aligner_list)
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.AnnotSV.Design.vcf.gz"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.AnnotSV.Design.bcftoolsmerge.log"


# Panel tsv individual samples AnnotSV
use rule AnnotSV as AnnotSV_panel with:
		input: rules.vcf_normalization.output
		output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv",
		log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.log"

use rule wait_for_AnnotSV as wait_for_AnnotSV_panel with:
	input:
		output_from_AnnotSV=rules.AnnotSV.output,
		log_file=f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.log"
	output:
		ready=f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.ready"

# Panel tsv all samples AnnotSV
use rule merge_tsv as merge_tsv_panel with:
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv", sample=sample_list, aligner=aligner_list, panel=panels_list)
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv"

use rule correct_tsv as correct_tsv_panel with:
	input:
		AnnotSV=rules.AnnotSV_panel.output,
		ready=rules.wait_for_AnnotSV_panel.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel_corr.{{panel}}.tsv")

use rule correct_chr as correct_chr_panel with:
	input: rules.correct_tsv_panel.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel_chr.{{panel}}.tsv")

use rule variantconvert as variantconvert_panel with:
	input: rules.correct_chr_panel.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel_uncorr.{{panel}}.vcf")
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.varianconvert.log"

use rule correct_vcf as correct_vcf_panel with:
	input: rules.variantconvert_panel.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel_unfix.{{panel}}.vcf")

use rule sort_vcf as sort_vcf_panel with:
	input: rules.correct_vcf_panel.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel_sort.{{panel}}.vcf")

# Panel vcf.gz individual samples AnnotSV
use rule fix_vcf as fix_vcf_panel with:
	input: rules.sort_vcf_panel.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.vcf.gz"

# Panel vcf.gz all samples AnnotSV
use rule merge_vcf as merge_vcf_panel with:
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.vcf.gz", sample=sample_list, aligner=aligner_list, panel=panels_list)
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.AnnotSV.Panel.{{panel}}.vcf.gz"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.AnnotSV.Panel.{{panel}}.bcftoolsmerge.log"


onstart:
	print('[INFO] Starting SCRAMBLE pipeline, brace yourselves!')
	shell(f"touch {os.path.join(outputDir, f'{serviceName}Running.txt')}")
	with open(logfile, "a+") as f:
		f.write("\n")
		f.write("Global parameters of the analysis for debug only")
		json.dump(config, f, ensure_ascii=False, indent=2)
		f.write("\n")

onsuccess:
	print('[INFO] SCRAMBLE pipeline finished successfully!')
	include = config['INCLUDE_RSYNC']
	shell(f"touch {outputDir}/{serviceName}Complete.txt")
	shell(f"rm -f {outputDir}/{serviceName}Running.txt")
	date_time_end = datetime.now().strftime("%Y%m%d-%H%M%S")
	with open(logfile, "a+") as f:
		f.write(f"End of the analysis : {date_time_end}\n")
	
	# Generate dictionary for results
	search_args = [arg.format(serviceName=serviceName) if '{serviceName}' in arg else arg for arg in ["/*/{serviceName}/*/*", "/*"]]
	result_files_list = searchfiles(resultDir, search_args, False)
	replaced_paths = replace_path(result_files_list, resultDir, "")
	sample_list.insert(0,"allsamples")
	resultDict = populate_dictionary(sample_list, config['RESULT_EXT_LIST'], replaced_paths, pattern_include=serviceName, split_index=2)	
	update_results(runDict, resultDict, config['RESULT_EXT_LIST'], exclude_samples=["allsamples"], exclude_keys=["hpo", "gender"])
	print('[INFO] Generating html report')
	generate_html_report(resultDict, runName, serviceName, sample_list, f"{serviceName}.template.html" , f"{resultDir}/{serviceName}.{date_time}.report.html")
	copy2(config['TEMPLATE_DIR'] + '/' + serviceName + '.style.css', resultDir)
	print('[INFO] Generating html report done')

	print('[INFO] Removing old results')
	for sample in sample_list:
		shell(f"rm -f {outputDir}/{sample}/{serviceName}/* || true")
	print('[INFO] Copying files')
	shell("rsync -azvh --include={include} --exclude='*' {resultDir}/ {outputDir}")
	for sample in sample_list:
		shell(f"cp {outputDir}/{sample}/{serviceName}/{sample}_{date_time}_{serviceName}/* {outputDir}/{sample}/{serviceName}/ || true")
	print('[INFO] Copying files done')

	# Optionally, perform DEPOT_DIR copy
	if config['DEPOT_DIR'] and outputDir != depotDir:
		print('[INFO] Removing old results from archives')
		for sample in sample_list:
			shell(f"rm -f {depotDir}/{sample}/{serviceName}/* || true")
		print('[INFO] Copying files into archives')
		shell("rsync -azvh --include={include} --exclude='*' {resultDir}/ {depotDir}")
		for sample in sample_list:
			shell(f"cp {outputDir}/{sample}/{serviceName}/{sample}_{date_time}_{serviceName}/* {depotDir}/{sample}/{serviceName}/ || true")
		print('[INFO] Copying files into archives done')

onerror:
	print('[ERROR] SCRAMBLE pipeline did not end well, grab a cup of coffee and check for log and err files')
	include_log = config['INCLUDE_LOG_RSYNC']
	shell(f"touch {outputDir}/{serviceName}Failed.txt")
	shell(f"rm -f {outputDir}/{serviceName}Running.txt")
	shell("rsync -azvh --include={include_log} --exclude='*' {resultDir}/ {outputDir}")