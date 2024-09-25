##########################################################################
# Snakemakefile Version:   2.0
# Description:             Snakemake file to run CANOES module https://github.com/ShenLab/CANOES
##########################################################################

# PROD version 1 : 03/06/2022
# Authoring : Jean Baptiste LAMOUCHE

# PROD version 2 : 14/03/2024 changelog
# Authoring : Thomas LAVAUX
	# AnnotSV version 3.4 (include the vcf converter) HPO ready with exomiser
	# options to process bed, k-merisation
	# convert R to R scripts and rewrite some code 
	# re-arrange/simplify/compact code a lot, input files copying is within the rules, f-string, etc.
	# add an html report using jinja2
	# remove header from tsv converted
	# addruns option

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
from os.path import join as osj

################## Configuration file ##################
configfile: "/app/config/snakefile/canoes_default.yaml"

####################### FUNCTIONS ######################

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

def create_list(txtlistoutput, file_list, ext_list, pattern):
	""" Function to create a txt file from a list of files filtered by extension and pattern """
	with open(txtlistoutput, 'a+') as f:
		f.writelines(f"{files}\n" for files in file_list if any(files.endswith(ext) and pattern in files for ext in ext_list))

def extract_tag(tagfile, tag, tagsep, sep):
	""" Function to extract a tag """
	row = open(tagfile, 'r').readline().strip().split(tagsep)
	output_tag = next((items.split(sep)[-1] for items in row if tag in items and sep in items), "")
	return output_tag

def kmerisation(kmerSize, bedFile, kmerBedFile):
	""" Bed kmerisation """
	print(f"{kmerSize} kmerisation of your bed {bedFile} in progress")
	
	with open(bedFile, 'r') as readBed, open(kmerBedFile, 'w+') as writeKbed:
		for line in readBed:
			if line.startswith('#'):
				continue
			
			chr, start, end, gene = line.split()[:4]
			diff = int(end) - int(start)
			
			while diff >= kmerSize:
				newEnd = int(start) + kmerSize - 1
				writeKbed.write(f"{chr}\t{start}\t{newEnd}\t{gene}\n")
				start, diff = newEnd + 1, diff - kmerSize
			
			if diff > 0:
				writeKbed.write(f"{chr}\t{start}\t{end}\t{gene}\n")
	
	print('Kmerisation done')

def replace_path(file_paths, old_substring, new_substring):
	return [path.replace(old_substring, new_substring).lstrip("/") for path in file_paths]

def generate_html_report(result_dict, run_name, service_name, sample_list, template_name, output_file='report.html', sample_list_added=None, gender_list=None, list_XX=None, list_XY=None, none_gender_samples=None):
	env = Environment(loader=FileSystemLoader(config['TEMPLATE_DIR']))
	template = env.get_template(template_name)

	rendered_html = template.render(
		runDict=result_dict,
		runName=run_name,
		serviceName=service_name,
		sample_list=sample_list,
		sample_list_added=sample_list_added,
		gender_list=gender_list,
		list_XX=list_XX,
		list_XY=list_XY,
		none_gender_samples=none_gender_samples
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

def process_gene_list(file, source_dir, dest_dir, panels_list):
	inputfile = os.path.join(source_dir, f"{file}.bed") 
	panel_trunc = file.split(".", 1)[1].split(".genes", 1)[0] if "." in file else file 
	outputfile = os.path.join(dest_dir, panel_trunc) 
	copy2(inputfile, outputfile) 
	panels_list.append(panel_trunc)
	return panels_list

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
# Definition of paths and name files
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

# Search files in repository 
files_list = searchfiles(os.path.normpath(config['run']), config['SEARCH_ARGUMENT'],  config['RECURSIVE_SEARCH'])

# Create sample and aligner list
sample_list = extractlistfromfiles(files_list, config['PROCESS_FILE'], '.', 0)
aligner_list = extractlistfromfiles(files_list, config['PROCESS_FILE'], '.', 1)

# Exclude samples from the exclude_list, case insensitive
sample_list = [sample for sample in sample_list if not any(sample.upper().startswith(exclude.upper()) for exclude in config['EXCLUDE_SAMPLE'])]

# If filter_sample_list variable is not empty, it will force the sample list
if config['FILTER_SAMPLE']:
	sample_list = list(config['FILTER_SAMPLE'])

# For validation analyse bam will be sample.aligner.validation.bam, so we append .validation to all the aligner strings
if config['VALIDATION_ONLY']:
	filtered_files = filter_files(files_list, filter_in='validation', extensions=config['PROCESS_FILE'])
	append_aligner = '.validation'
	aligner_list = [sub + append_aligner for sub in aligner_list]
else:
	filtered_files = filter_files(files_list, None ,filter_out='validation', extensions=config['PROCESS_FILE'])

# Populate dictionary
runDict = populate_dictionary(sample_list, config['EXT_INDEX_LIST'], filtered_files, None, ['analysis'])

# Add other runs/samples to improve sensibility
if config['addruns']:
	files_list_addrun = []  
	for run in config['addruns']:
		result = searchfiles(os.path.normpath(run), config['SEARCH_ARGUMENT'], config['RECURSIVE_SEARCH'])
		files_list_addrun.extend(result)
	# We generate the sample list with the ext of the file in the PROCESS_FILE varialbe (bam/cram for ex)
	sample_list_addruns = extractlistfromfiles(files_list_addrun, config['PROCESS_FILE'], '.', 0)
	if config['VALIDATION_ONLY']: 	# Filter the list if validation only
		filter_files(files_list_addrun, filter_in='validation', extensions=config['PROCESS_FILE'])
	# We populate the dictionary with the same EXT_INDEX_LIST variable as the original list, to get HPO and gender is avalaible
	addrunDict = populate_dictionary(sample_list_addruns, config['EXT_INDEX_LIST'], files_list_addrun, None, ['analysis'])

	for key, value in addrunDict.items(): 	# Merge runDict dictionary with the addrunDict
		runDict[key].update(value)
		sample_list = sample_list + sample_list_addruns # Concatenate all samples for the analysis
		sample_list_to_copy = [x for x in sample_list if x not in sample_list_addruns]
else:
	sample_list_addruns = None
	sample_list_to_copy = sample_list.copy()

# Set a filelist with all the files tag ; file format is sample.tag
# tag ex SEX#M!POOL#POOL_HFV72AFX3_M_10#POOL_HFV72AFX3_F_11!
for individual in runDict.values():
	individual.update({'hpo': None, 'gender': None})
ped_file = find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.ped')

if os.path.exists(ped_file) and os.path.getsize(ped_file) > 0:
	pedigree_into_dict(ped_file, sample_list, runDict, individual_id_column='Individual ID', hpo_list_column='HPOList', sex_column='Sex')
	print('[INFO] Gender found in the pedigree file')

	df = pd.read_csv(ped_file, sep='\t', dtype=str)
	for _, row in df.iterrows():
		individual_id = row['Individual ID']
		hpo_list = str(row['HPOList']).replace(" ", "")
		if individual_id in sample_list:
			runDict[individual_id]['hpo'] = hpo_list
else:
	tagfile_list = [runDict[sample].get('.tag') for sample in sample_list if '.tag' in runDict.get(sample, {})]

	if not tagfile_list:
		print('[INFO] No gender found for the samples')
	else:
		for tagfile in tagfile_list:
			sample_name = os.path.basename(tagfile).split(".")[0]
			runDict[sample_name]['gender'] = extract_tag(tagfile, 'SEX', '!', '#')
		print('[INFO] Gender found in the tag files')

# Extract the gender_list from dictionary with the key gender
gender_list = ['A'] + ['XX' if runDict[sample]['gender'] == 'F' else 'XY' for sample in sample_list if sample in runDict and 'gender' in runDict[sample]]
gender_list = list(set(gender_list)) # Removing duplicate

# Find bed file (Design)
config['BED_FILE'] = config['BED_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.design.bed', '.genes.bed')
if not config['BED_FILE']:
	print('No bed found, CANOES cannot continue, exiting')
	exit()
# Find genes file (Panel); we can't use .genes files because .list.genes and .genes are not distinctable from the indexing we made
config['GENES_FILE'] = config['GENES_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.genes.bed', '.list.genes')
# Find transcripts files (NM)
config['TRANSCRIPTS_FILE'] = config['TRANSCRIPTS_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.transcripts', '.list.transcripts')

# If transcript file exist, create the annotation file for AnnotSV
annotation_file = f"{resultDir}/{serviceName}.{date_time}.AnnotSV.txt"
if os.path.exists(config['TRANSCRIPTS_FILE']) and os.path.getsize(config['TRANSCRIPTS_FILE']):
	df = pd.read_csv(config['TRANSCRIPTS_FILE'], sep='\t', names=["NM", "Gene"])
	with open(annotation_file, 'w+') as f:
		f.write('\t'.join(df['NM']))
else:
	with open(annotation_file, 'w') as f:
		f.write("No NM found")

# Find list.genes files 
config['LIST_GENES'] = config['LIST_GENES'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.list.genes', '.list.transcripts')

# Transform list_genes into a list if list_genes exist, else use genes_file if exist
panels_list = []
file_source = config.get('LIST_GENES') or config.get('GENES_FILE')

if file_source:
	with open(file_source) as f:
		files = f.read().splitlines() if config['LIST_GENES'] else [os.path.splitext(os.path.basename(file_source))[0]]
	for file in files:
		panels_list = process_gene_list(file, os.path.dirname(file_source), resultDir, panels_list)

# Kmerisation
canoesbed_file = f"{resultDir}/{serviceName}.{date_time}.bed"
if config['KMER']:
	kmerisation(config['KMER'], config['BED_FILE'], canoesbed_file)
else:
	copy2(config['BED_FILE'], canoesbed_file)

# Create file_list of bam by gender with the runDict, depending on the gender_list
# A = all ; XX = Female only ; XY = Male only
# files_list_A contains the full path of all files (bam files for ex)
# Warning the key dictionary for sexe is A/M/F but the gender_list is A/XY/XX
files_list_A = list(set([runDict[sample]['.bam'] for sample in sample_list if sample in runDict and '.bam' in runDict[sample]]))
files_list_XX = list(set([files for files in files_list_A if runDict.get(os.path.basename(files).split(".")[0], {}).get('gender') == 'F']))
files_list_XY = list(set([files for files in files_list_A if runDict.get(os.path.basename(files).split(".")[0], {}).get('gender') == 'M']))

# Creating a txt list for the bam files per aligner per gender (A, M & F)
for aligner in aligner_list:
	for gender in gender_list:
		bamlist = f"{resultDir}/{serviceName}.{date_time}.{gender}.list.txt"
		create_list(bamlist, globals()[f"files_list_{gender}"], config['PROCESS_FILE'], aligner)

# Option to remove gender in the gender_list
gender_list = sorted([gender for gender in gender_list if gender not in config['REMOVE_GENDER']])

# check the ref bam list and presence of gender, and update the gender_list based on the presence of 'F' and 'M'
if config['REF_BAM_LIST']:
	df = pd.read_csv(config['REF_BAM_LIST'], sep='\t')
	if not df['gender'].str.contains('F').any():
		gender_list.remove('XX')
	if not df['gender'].str.contains('M').any():
		gender_list.remove('XY')

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
# check the number of sample for copy or merge vcf rule
sample_count = len(sample_list) 

# Priority order
ruleorder: copy_bam > copy_cram > cramtobam > indexing

rule all:
	"""Output a design vcf.gz with the bams and the bed provided, optionally using a reference set and generating plots if specified"""
	input:
		expand(f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.CNVCall.tsv", aligner=aligner_list, gender=gender_list) +
		expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.vcf.gz", aligner=aligner_list, sample=sample_list) +
		[f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Design.tsv"] +
		(expand(f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.refsamples.reads.tsv", aligner=aligner_list) if config['REF_BAM_LIST'] else []) +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.barplot.pdf", sample=sample_list, aligner=aligner_list) if config['plots'] else []) +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.boxplot.pdf", sample=sample_list, aligner=aligner_list) if config['plots'] else []) +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv", aligner=aligner_list, sample=sample_list, panel=panels_list) if config['GENES_FILE'] else []) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel.{{panel}}.tsv", panel=panels_list) if config['GENES_FILE'] else []) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel.{{panel}}.vcf.gz", aligner=aligner_list, panel=panels_list) if config['GENES_FILE'] else [])


rule help:
	"""
	General help for CANOES
	Launch snakemake -s snakefile_canoes -c(numberofthreads) --config run=absolutepathoftherundirectory
	To launch the snakemake file, use --config to replace variables that must be properly set for the pipeline to work ie run path directory
	Every variable defined in the yaml file can be change
	Separate multiple variable with a space (ex  --config run=runname transProb=0.05 var1=0.05 var2=12)
	Use option --configfile another.yaml to replace and merge existing default config.yaml file variables
	Use -p to display shell commands
	Use --lt to display docstrings of rules
	Use -n for a dry run
	Input file = bam or cram files (if bai is needed, it will be generate)
	Output file = vcf annoted with AnnotSV 3.x for each sample/bam, and a global vcf file will all the samples ; a set of vcf files by design/panel
	"""

rule copy_bam:
	output: temp(f"{resultDir}/{{sample}}.{{aligner}}.bam")
	params:
		process=config['PROCESS_CMD'],
		download_link=lambda wildcards: runDict[wildcards.sample]['.bam']
	shell: "[ \"{params.process}\" = \"ln\" ] && ln -sfn {params.download_link} {output} || rsync -azvh {params.download_link} {output}"

rule copy_cram:
	output: temp(f"{resultDir}/{{sample}}.{{aligner}}.cram")
	params:
		process=config['PROCESS_CMD'],
		download_link=lambda wildcards: runDict[wildcards.sample]['.cram']
	shell: "[ \"{params.process}\" = \"ln\" ] && ln -sfn {params.download_link} {output} || rsync -azvh {params.download_link} {output}"

rule cramtobam:
	""" Extract bam from a cram file with samtools, need a reference genome """
	input: rules.copy_cram.output
	output: temp(f"{resultDir}/{{sample}}.{{aligner}}.bam")
	params: config['REFGENEFA_PATH']
	shell: "samtools view -b -T {params} -o {output} {input}"

rule indexing:
	""" Indexing bam files with samtools or ln """
	input: f"{resultDir}/{{sample}}.{{aligner}}.bam"
	output: temp(f"{resultDir}/{{sample}}.{{aligner}}.bam.bai")
	params:
		process=config['PROCESS_CMD'],
		download_link=lambda wildcards: runDict[wildcards.sample]['.bam.bai']
	threads: workflow.cores
	shell: "[ \"{params.process}\" = \"ln\" ] && ln -sfn {params.download_link} {output} || samtools index -b -@ {threads} {input} {output}"

rule gc_percent:
	""" Compute GC percent from the reference genome and a bed file, removing header 
		Output a tsv with 4 columns : chromosome start stop gc	
	"""
	input: config['REFGENEFA_PATH']
	output: 
		gc_percent_header=temp(f"{resultDir}/GCpercent_header.tsv"),
		gc_percent = f"{resultDir}/GCpercent.tsv"
	params: 
		java_option = config['JAVA_OPTION']
	log: f"{resultDir}/GCpercent.log"
	shell: 
		"""
			gatk --java-options {params.java_option} AnnotateIntervals -L {canoesbed_file} -R {input} -imr OVERLAPPING_ONLY -O {output.gc_percent_header} 2> {log}
			grep -v '^@' {output.gc_percent_header} > {output.gc_percent}
		 """

rule bam_to_multicov:
	""" From each sample bam file generate a coverage file filtered with the bed provided """
	input: 
		bam = f"{resultDir}/{{sample}}.{{aligner}}.bam",
		bai = f"{resultDir}/{{sample}}.{{aligner}}.bam.bai"
	output:	f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.multicov.tsv"
	params: pyscripts = config['PY_SCRIPTS']
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.multicov.log"
	shell:" python {params.pyscripts}/pybedtools.py --bam {input.bam} --bed {canoesbed_file} -q 20 -o {output} 1> {log} "


rule bam_to_multicov_for_refbam:
	""" From a tsv list of bams generate a coverage file filtered with the bed provided, including the header """
	input: config['REF_BAM_LIST']
	output:	f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.refsamples.reads.tsv"
	params: pyscripts = config['PY_SCRIPTS']
	log: f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.refsamples.reads.log"
	shell:" python {params.pyscripts}/pybedtools.py --tsv {input} --bed {canoesbed_file} -q 20 -o {output} 1> {log} "


rule merge_multicov:
	""" Merge all coverage files and add header """
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.multicov.tsv", aligner=aligner_list, sample=sample_list)
	output: f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.allsamples.reads.tsv"
	log: f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.allsamples.reads.log"
	params: config['PY_SCRIPTS']
	shell: " python {params}/merge_multicov.py {input} -o {output} "


rule plot_coverage_stats:
	""" Generate coverage calculation """
	input:  rules.merge_multicov.output
	output:	f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.allsamples.multicoverage.stats.tsv"
	params: config['PY_SCRIPTS']
	shell: " python {params}/plot_coverage_stats.py {input} {output} "


rule plot_coverage:
	""" Generate coverage plots per sample in pdf """
	input: rules.plot_coverage_stats.output
	output: 
		barplot = f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.barplot.pdf",
		boxplot = f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.boxplot.pdf"
	params: config['R_SCRIPTS']
	log:
		log=f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.plotCoverage.log",
		err=f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.plotCoverage.err"
	shell: " Rscript {params}/plotCoverage.R --input {input} --output1 {output.barplot} --output2 {output.boxplot} 1> {log.log} 2> {log.err} "

rule canoes_calling:
	""" CANOES analysis """
	input: 
		read = rules.merge_multicov.output,
		gc = rules.gc_percent.output.gc_percent
	params:
		Rscripts = config['R_SCRIPTS'],
		chromosome= "{gender}",
		bamlist= f"{resultDir}/{serviceName}.{date_time}.{{gender}}.list.txt",
		pvalue = config['pval'],
		distance = config['dist'],
		tnumeric = config['tnum'],
		numreference = config['numref'],
		hom = config['homdel'],
		refbamlist= "--refbams {}".format(config['REF_BAM_LIST']) if config['REF_BAM_LIST'] else "",
		refmulticovtsv = f" --readsrefs {resultDir}/{serviceName}.{date_time}.{{aligner}}.refsamples.reads.tsv" if config['REF_BAM_LIST'] else ""
	output:
		cnvcall = f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.CNVCall.tsv",
		success = f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.CANOEScalling.success",
		rdata = f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.CANOES.Rdata"
	log: 
		log = f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.CANOEScalling.log",
		err = f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.CANOEScalling.err"
	shell:
		"""
		Rscript {params.Rscripts}/CANOES.v2.2.R --gcfile {input.gc} --readsfile {input.read} --chromosome {params.chromosome} \
		--samples {params.bamlist} --homdel {params.hom} --numref {params.numreference} --tnum {params.tnumeric} --distance {params.distance} \
		--pvalue {params.pvalue} {params.refbamlist} {params.refmulticovtsv} --output {output.cnvcall} --rdata {output.rdata} 1> {log.log} 2> {log.err} && touch {output.success}
		"""

rule merge_makeCNVcalls:
	""" Merge all CNV call results """
	input: expand(f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.CNVCall.tsv", gender=gender_list, aligner=aligner_list)
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.CNVCall_1based.tsv")
	shell: " (head -n 1 {input[0]} && tail -n +2 -q {input} | sort -u | sort -r) > {output} "

rule correct_1based:
	input: rules.merge_makeCNVcalls.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.CNVCall.tsv"
	shell:
		"""
		awk 'BEGIN {{FS=OFS="\\t"}} NR==1 {{print $0}} NR>1 {{split($6,a,":|-"); $2=$2+1; $6=a[1]":"a[2]+1"-"a[3]; print $0}}' {input} > {output}
		"""

rule bedtovcf:
	""" Convert CANOES bed to vcf """
	input: rules.correct_1based.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.vcf")
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.bedtovcf.log"
	params:	config['VARIANT_CONVERT_CONFIG']
	shell: " variantconvert convert -i {input} -o {output} -c {params}/config_canoes_bed.json 2> {log} "


rule sortvcf:
	""" Sort vcf """
	input: rules.bedtovcf.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design.sort.vcf")
	shell: " {{ grep '^#' {input} && grep -v '^#' {input} | sort -k1,1V -k2,2g; }} > {output} "


rule add_chr_to_vcf:
	input: rules.sortvcf.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design.chr.vcf")
	shell: "awk -F '\\t' 'BEGIN {{OFS=\"\\t\"}} /^##/ {{print}} /^#CHROM/ {{print}} /^[^#]/ {{$1 = \"chr\" $1; print}}' {input} > {output}"


rule fix_svlen:
	input: rules.add_chr_to_vcf.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design.vcf.gz"
	params:	config['PY_SCRIPTS']
	shell: " python {params}/multifix_vcf.py -i {input} -o {output} -z ; tabix {output} "


rule split_vcf:
	"""
	Split vcf into each individual sample with bcftools
	-s {sample}: comma-separated list of samples to include
	-Oz: output vcf compressed
	-c1: minimum allele count (INFO/AC) of sites to be printed
	"""
	input: rules.fix_svlen.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Design.Noannotation.vcf.gz")
	params: config['DUMMY_FILES']
	shell: "mkdir -p {resultDir}/{wildcards.sample}/{serviceName} && bcftools view -c1 -Oz -s {wildcards.sample} -o {output} {input} && [[ -s {output} ]] || cat {params}/empty.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' | bgzip > {output} ; tabix {output}"


rule AnnotSV:
	"""
	Annotate and rank Structural Variations from a vcf file
	Output a vcf file, if no annotation is done output an empty vcf
	- annotationMode can be split by exons/introns or full by genes
	- txtFile: path to a file containing a list of preferred genes transcripts for annotation
	- genomeBuild must be specified if not hg38
	"""
	input: rules.split_vcf.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.tsv"
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.log"
	params:
		dummypath=config['DUMMY_FILES'],
		outputfile=f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design",
		genome=config['genomeBuild'],
		overlap=config['overlap'],
		mode=config['annotationMode'],
		annotation=config['annotationdir'],
		hpo=lambda wildcards: runDict[wildcards.sample]['hpo'],
		samplevcf=lambda wildcards: runDict[wildcards.sample][vcf_extension] # -snvIndelSamples to specify samples in the vcf
	shell: 
		"""
		AnnotSV -SVinputFile {input} -outputFile {params.outputfile} -snvIndelFiles {params.samplevcf} -annotationMode {params.mode} -annotationsDir {params.annotation} -hpo {params.hpo} -txFile {annotation_file} -genomeBuild {params.genome} -overlap {params.overlap} > {log} && [[ -s {output} ]]  || cat {params.dummypath}/emptyAnnotSV.tsv | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output}
		"""

rule correct_tsv:
	input: rules.AnnotSV.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_fix.tsv")
	shell: "sed '1s/[()/\\]/_/g' {input} > {output}"

rule correct_chr:
	"""
	Add chr to the SV_chrom column
	"""
	input: rules.correct_tsv.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_chrfix.tsv")
	shell:
		"""
		awk '{{{{FS=OFS="\t"}};if(NR==1){{print; next}}; $2="chr"$2; print}}' {input} > {output}
		"""

rule variantconvert:
	input: rules.correct_chr.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_unfix.vcf")
	params: 
			variantconvertconfig=config['VARIANT_CONVERT_CONFIG'],
			dummypath=config['DUMMY_FILES']
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.varianconvert.log"
	shell:
		"""
		variantconvert convert -i {input} -o {output} -c {params.variantconvertconfig}/annotsv3_from_vcf.json 2> {log} && [[ -s {output} ]] || cat {params.dummypath}/emptyAnnotSV.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output}
		"""

rule correctvcf:
	input: rules.variantconvert.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_fix.vcf")
	shell: " sed 's/SAMPLENAME/{wildcards.sample}/g' {input} > {output} "

rule sortvcf_sample:
	input: rules.correctvcf.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_sort.vcf")
	params: config['DUMMY_FILES']
	shell: " {{ grep \'^#\' {input} && grep -v \'^#\' {input} | sort -k1,1V -k2,2g; }} > {output} && [[ -s {output} ]] || cat {params}/empty.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output} "

rule fix_vcf:
	input: rules.sortvcf_sample.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.vcf.gz"
	params:	config['PY_SCRIPTS']
	shell: "python {params}/multifix_vcf.py -i {input} -o {output} -z ; tabix {output}"

rule merge_vcf:
	""" Copy or merge multiple vcfs using bcftools """
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.vcf.gz", sample=sample_list, aligner=aligner_list)
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Design.vcf.gz"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Design.bcftoolsmerge.log"
	shell: "if [ {sample_count} -eq 1 ]; then cp {input} {output}; else bcftools merge {input} -O z -o {output} 2> {log}; fi; tabix {output}"

rule vcf2tsv:
	""" vcf to tsv conversion """
	input: rules.merge_vcf.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Design.tsv"
	shell:
		"""
		if [ $(gunzip -c {input} | wc -c) -eq 0 ]; then
			echo -e "No data available in the input file" > {output}
		else
			vcf2tsvpy --keep_rejected_calls --input_vcf {input} --out_tsv {output}.tmp
			grep -v '^#' {output}.tmp > {output}
		fi
		"""

use rule vcf2tsv as vcf2tsv_sample with:
	input: rules.fix_vcf.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.tsv"

rule filter_vcf:
	"""	Filter vcf with a bed file """
	input: rules.merge_vcf.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel_unnorm.{{panel}}.vcf.gz")
	params: lambda wildcards: f"{resultDir}/{wildcards.panel}"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel.{{panel}}.bedtoolsfilter.log"
	shell: "bedtools intersect -header -a {input} -b {params} 2> {log} | bgzip > {output}"

rule vcf_normalisation:
	input: rules.filter_vcf.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel.{{panel}}.vcf.gz"
	shell: "bcftools norm -d all -o {output} -Oz {input} ; tabix {output}"

use rule vcf2tsv as vcf2tsv_panel with:
	input: rules.vcf_normalisation.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel.{{panel}}.tsv"

use rule split_vcf as split_vcf_panel with:
	input: rules.fix_vcf.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.vcf.gz"

use rule vcf2tsv as vcf2tsv_panel_sample with:
	input: rules.split_vcf_panel.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv"

use rule filter_vcf as filter_vcf_panel with:
	input: rules.fix_svlen.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel_unnorm.{{panel}}.vcf.gz")
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel.{{panel}}.bedtoolsfilter.log"

use rule vcf_normalisation as vcf_normalisation_all with:
	input: rules.filter_vcf_panel.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel.{{panel}}.vcf.gz"


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
	
	# Removing additional samples from the results
	if sample_list_addruns:
		for sample in sample_list_addruns:
			shell(f"rm -rf {resultDir}/{sample} || true")
			#shell(f"find {outputDir}/{sample}/{serviceName} -type f -exec rm -f {{}} +")
		# We extract the samples that are present in the sample_list_addruns & the vcf
		shell(f"bcftools query -l {resultDir}/{serviceName}.{date_time}.allsamples.Design.unfiltered.vcf.gz > vcf_samples.txt")
		with open('final_sample_list.txt', 'w') as file:
			file.write('\n'.join(sample_list_addruns))
		# Then we filter out the vcf
		shell(f"bcftools view -S final_sample_list.txt -Oz -o {resultDir}/{serviceName}.{date_time}.allsamples.Design.filtered.vcf.gz {resultDir}/{serviceName}.{date_time}.allsamples.Design.unfiltered.vcf.gz")
		shell(f"bcftools view -c 1 -Oz -o {resultDir}/{serviceName}.{date_time}.allsamples.Design.vcf.gz {resultDir}/{serviceName}.{date_time}.allsamples.Design.filtered.vcf.gz")
		# And convert vcf to tsv
		shell(f"vcf2tsvpy --keep_rejected_calls --input_vcf {resultDir}/{serviceName}.{date_time}.allsamples.Design.vcf.gz --out_tsv {resultDir}/{serviceName}.{date_time}.allsamples.Design.tsv.tmp && cat {resultDir}/{serviceName}.{date_time}.allsamples.Design.tsv.tmp | grep -v '^#' > {resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Design.tsv ")


	# Generate dictionary from outputdir
	search_args = [arg.format(serviceName=serviceName) if '{serviceName}' in arg else arg for arg in ["/*/{serviceName}/*/*", "/*"]]
	result_files_list = searchfiles(resultDir, search_args, False)
	replaced_paths = replace_path(result_files_list, resultDir, "")
	sample_list_to_copy.insert(0,"allsamples")
	resultDict = populate_dictionary(sample_list_to_copy, config['RESULT_EXT_LIST'], replaced_paths, pattern_include=serviceName, split_index=2)	

	# Update resultDict with values from runDict
	update_results(runDict, resultDict, config['RESULT_EXT_LIST'], exclude_samples=["allsamples"], exclude_keys=["hpo", "gender"])
	list_XX = [os.path.basename(file).split('.')[0] for file in files_list_XX]
	list_XY = [os.path.basename(file).split('.')[0] for file in files_list_XY]
	none_gender_samples = [sample for sample, data in runDict.items() if data.get('gender') == 'None']
	generate_html_report(resultDict, runName, serviceName, sample_list, f"{serviceName}.template.html" , f"{resultDir}/{serviceName}.{date_time}.report.html", sample_list_added=sample_list_addruns, gender_list=gender_list, list_XX=list_XX, list_XY=list_XY, none_gender_samples=none_gender_samples)
	copy2(config['TEMPLATE_DIR'] + '/' + serviceName + '.style.css', resultDir)

	# Clear existing output directories & copy
	for sample in sample_list_to_copy:
		shell(f"rm -f {outputDir}/{sample}/{serviceName}/* || true")

	shell("rsync -azvh --include={include} --exclude='*' {resultDir}/ {outputDir}")
	for sample in sample_list_to_copy:
		shell(f"cp {outputDir}/{sample}/{serviceName}/{sample}_{date_time}_{serviceName}/* {outputDir}/{sample}/{serviceName}/ || true")

	# Optionally, perform DEPOT_DIR copy
	if config['DEPOT_DIR']:
		if outputDir != depotDir:
			shell("rsync -azvh --include={include} --exclude='*' {resultDir}/ {depotDir}")
			for sample in sample_list:
				shell(f"cp {outputDir}/{sample}/{serviceName}/{sample}_{date_time}_{serviceName}/* {depotDir}/{sample}/{serviceName}/ || true")

onerror:
	include_log = config['INCLUDE_LOG_RSYNC']
	shell(f"touch {outputDir}/{serviceName}Failed.txt")
	shell(f"rm -f {outputDir}/{serviceName}Running.txt")
	shell("rsync -azvh --include={include_log} --exclude='*' {resultDir}/ {outputDir}")