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
	# remove header from tsv converted
	
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

def process_gene_list(file, source_dir, dest_dir, panels_list):
	"""
    Processes a gene list file by copying it from the source directory to the destination directory 
    and appending a truncated version of the filename to a list of panels.

    Args:
        file (str): The name of the gene list file (without the '.bed' extension).
        source_dir (str): The directory where the input file is located.
        dest_dir (str): The directory where the processed file should be saved.
        panels_list (list): A list to which the truncated panel name will be appended.

    Returns:
        list: The updated list of panels with the newly appended panel name.
    """
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

# Exclude samples from the exclude_list , case insensitive
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

# Find bed file (Design)
config['BED_FILE'] = config['BED_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.design.bed', '.genes.bed')
# Find genes file (Panel); we can't use .genes files because .list.genes and .genes are not distinctable from the indexing we made
config['GENES_FILE'] = config['GENES_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.genes.bed', '.list.genes')
# Find list.genes files (a list of panel files)
config['LIST_GENES'] = config['LIST_GENES'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.list.genes', '.list.transcripts')
# Find transcripts files (file containing NM)
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

# Transform list_genes into a list if list_genes exist, else use genes_file if exist
panels_list = []
file_source = config.get('LIST_GENES') or config.get('GENES_FILE')
if file_source:
	with open(file_source) as f:
		files = f.read().splitlines() if config['LIST_GENES'] else [os.path.splitext(os.path.basename(file_source))[0]]
	for file in files:
		panels_list = process_gene_list(file, os.path.dirname(file_source), resultDir, panels_list)

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
ruleorder: copy_bam > copy_cram > cramtobam > indexing > fix_vcf_full > bcftools_filter

rule all:
	"""Rule will create an unfiltered vcf.gz and corresponding tsv, with an optional design/panel vcf.gz & tsv if a bed/gene file is provided"""
	input:
		expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full.vcf.gz", sample=sample_list, aligner=aligner_list) +
		[f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Full.tsv"] +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.vcf.gz", sample=sample_list, aligner=aligner_list) if config['BED_FILE'] else []) +
		([f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Design.tsv"] if config['BED_FILE'] else []) +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv", aligner=aligner_list, sample=sample_list, panel=panels_list) if config['GENES_FILE'] else []) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel.{{panel}}.tsv", panel=panels_list) if config['GENES_FILE'] else []) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel.{{panel}}.vcf.gz", aligner=aligner_list, panel=panels_list) if config['GENES_FILE'] else [])

rule help:
	"""
	General help for SCRAMBLE module
	Launch snakemake -s  snakefile_scramble -c(numberofthreads) --config run=absolutepathoftherundirectory
	To launch the snakemake file, use --config to replace variables that must be properly set for the pipeline to work ie run path directory
	Every variable defined in the yaml file can be change
	Separate multiple variable with a space (ex  --config DATA_DIR=runname transProb=0.05 var1=0.05 var2=12)
	Use option --configfile another.yaml to replace and merge existing config.yaml file variables
	Use -p to display shell commands
	Use --lt to display docstrings of rules
	Use -n for a dry run
	Input file = cram or bam files (if bai is needed, it will be generate)
	Output file = vcf annoted with AnnotSV 3.x for each sample/bam, and a global vcf file will all the samples ; a set of vcf files by design/panel
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
	""" Indexing bam files with samtools or ln """
	input: f"{resultDir}/{{sample}}.{{aligner}}.bam"
	output: temp(f"{resultDir}/{{sample}}.{{aligner}}.bam.bai")
	params:
		process=config['PROCESS_CMD'],
		download_link=lambda wildcards: runDict[wildcards.sample]['.bam.bai']
	threads: workflow.cores
	message: """ Indexing {params.download_link} file """
	shell: "[ \"{params.process}\" = \"ln\" ] && ln -sfn {params.download_link} {output} || samtools index -b -@ {threads} {input} {output}"

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
		bam = f"{resultDir}/{{sample}}.{{aligner}}.bam",
		bai = f"{resultDir}/{{sample}}.{{aligner}}.bam.bai",
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
	Calling SCRAMble.R with --eval-meis produces a tab delimited file. If a genomereference.fa file is provided, then a VCF is produced as well.
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
		output = f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.raw"
	output:
		f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.raw.vcf"
	log: 
		log = f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.ScrambleR.log", 
		err = f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.ScrambleR.err"
	message: """ Starting calling of {input} file """
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

rule correct_vcf:
	"""	Correction of vcf output, add sample name and genotype to be consistent with the vcf format specification """
	input: rules.scramble.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.corr.vcf")
	params:
		mode = config['scramble_mode'],
		dummypath = config['DUMMY_FILES']
	shell:
		"""
		(grep "^##" {input} && echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' && grep "^#CHROM" {input} | awk -v SAMPLE={wildcards.sample} '{{print $0"\tFORMAT\t"SAMPLE}}' && grep "^#" -v {input} | awk '{{print $0"\tGT\t0/1"}}') > {output} && [[ -s {output} ]] || cat {params.dummypath}/empty.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output}
		"""

rule vcf2gz:
	input: rules.correct_vcf.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full_unfiltered.vcf.gz")
	params: config['DUMMY_FILES']
	shell: "bgzip -c {input} > {output} ; tabix {output}"

rule bcftools_filter:
	"""	Filter with bcftools """
	input: rules.vcf2gz.output	
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full.vcf.gz"
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.filter.log"
	params: bed=config['BCFTOOLS_FILTER'],
			dummy=config['DUMMY_FILES']
	shell: "bcftools view {params.bed} {input} -o {output}  2> {log} && [[ -s {output} ]] || cat {params.dummy}/empty.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' | bgzip > {output} ; tabix {output}" 

rule AnnotSV:
	"""
	Annotate and rank Structural Variations from a vcf file
	Output a vcf file, if no annotation is done output an empty vcf
	- annotationMode can be split by exons/introns or full by genes
	- txtFile: path to a file containing a list of preferred genes transcripts for annotation
	- genomeBuild must be specified if not hg38
	"""
	input: rules.bcftools_filter.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full.tsv"
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.log"
	params:
		outputfile=f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full",
		genome=config['genomeBuild'],
		overlap=config['overlap'],
		mode=config['annotationMode'],
		dummypath = config['DUMMY_FILES'],
		annotation=config['annotationdir'],
		hpo=lambda wildcards: runDict[wildcards.sample]['hpo'],
		samplevcf=lambda wildcards: runDict[wildcards.sample][vcf_extension] # -snvIndelSamples to specify samples in a multisample vcf
	shell: 
		"""
		AnnotSV -SVinputFile {input} -outputFile {params.outputfile} -snvIndelFiles {params.samplevcf} -annotationMode {params.mode} -annotationsDir {params.annotation} -hpo {params.hpo} -txFile {annotation_file} -genomeBuild {params.genome} -overlap {params.overlap} > {log} && [[ -s {output} ]]  || cat {params.dummypath}/emptyAnnotSV.tsv | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output}
		"""

rule correct_tsv:
	input: rules.AnnotSV.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full_unfix.tsv")
	shell: 
		""" 
		sed '1s/[()/\\]/_/g' {input} > {output}
		"""

# AnnotSV drop chr from input, so we need to re-add chr to the SV_chrom column (column number 2)
rule correct_chr:
	"""
	Add chr to the SV_chrom column
	"""
	input: rules.correct_tsv.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full_chrfix.tsv")
	shell:
		"""
		awk '{{{{FS=OFS="\t"}};if(NR==1){{print; next}}; $2="chr"$2; print}}' {input} > {output}
		"""

rule variantconvert:
	input: rules.correct_chr.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full_fix.vcf")
	params:
			VARIANT_CONVERT_CONFIG=config['VARIANT_CONVERT_CONFIG'],
			dummypath=config['DUMMY_FILES']
	log: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.varianconvert.log"
	shell:
		"""
		variantconvert convert -i {input} -o {output} -c {params.VARIANT_CONVERT_CONFIG}/annotsv3_from_vcf.json 2> {log} && [[ -s {output} ]] || cat {params.dummypath}/emptyAnnotSV.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output}
		"""

rule correct_vcf_bis:
	input: rules.variantconvert.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full_correct.vcf")
	shell:
		""" 
		sed 's/SAMPLENAME/{wildcards.sample}/g' {input} > {output}
		"""

rule sort_vcf:
	input: rules.correct_vcf_bis.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full_sort.vcf")
	params: config['DUMMY_FILES']
	shell: " {{ grep \'^#\' {input} && grep -v \'^#\' {input} | sort -k1,1V -k2,2g; }} > {output} && [[ -s {output} ]] || cat {params}/empty.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output} "

rule fix_vcf_full:
	input: rules.sort_vcf.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full.vcf.gz"
	params:	config['PY_SCRIPTS']
	shell: " python {params}/multifix_vcf.py -i {input} -o {output} -z ; tabix {output}"

rule merge_vcf_full:
	"""	Copy or merge vcfs with bcftools """
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full.vcf.gz", aligner=aligner_list, sample=sample_list)
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Full.vcf.gz"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Full.bcftoolsmerge.log"
	shell: "if [ {sample_count} -eq 1 ]; then cp {input} {output}; else bcftools merge {input} -O z -o {output} 2> {log}; fi; tabix {output}"

rule filter_vcf_design:
	input: rules.merge_vcf_full.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Design.vcf.gz"
	params: bed=config['BED_FILE'],
			dummypath=config['DUMMY_FILES']
	shell: " bedtools intersect -header -a {input} -b {params.bed} | bgzip > {output} && [[ -s {output} ]] || cat {params.dummypath}/empty.vcf | bgzip > {output} ; tabix {output} "

use rule filter_vcf_design as filter_vcf_panel with:
	input: rules.filter_vcf_design.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel_unnorm.{{panel}}.vcf.gz")
	params: bed=lambda wildcards: f"{resultDir}/{wildcards.panel}",
			dummypath=config['DUMMY_FILES']

rule vcf_normalisation_panel:
	input: rules.filter_vcf_panel.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel.{{panel}}.vcf.gz"
	shell: "bcftools norm -d all -o {output} -Oz {input} ; tabix {output}"

rule vcf2tsv:
	""" vcf to tsv conversion """
	input: rules.merge_vcf_full.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Full.tsv"
	shell:
		"""
		if [ $(gunzip -c {input} | wc -c) -eq 0 ]; then
			echo -e "No data available in the input file" > {output}
		else
			vcf2tsvpy --keep_rejected_calls --input_vcf {input} --out_tsv {output}.tmp
			grep -v '^#' {output}.tmp > {output}
		fi
		"""

use rule vcf2tsv as vcf2tsv_full with:
	input: rules.filter_vcf_design.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Design.tsv"

use rule vcf2tsv as vcf2tsv_panel with:
	input: rules.vcf_normalisation_panel.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel.{{panel}}.tsv"

use rule vcf2tsv as vcf2tsv_full_sample with:
	input: rules.fix_vcf_full.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Full.tsv"

rule split_vcf_design:
	input: rules.filter_vcf_design.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.vcf.gz"
	params: config['DUMMY_FILES']
	shell: "mkdir -p {resultDir}/{wildcards.sample}/{serviceName} && bcftools view -c1 -Oz -s {wildcards.sample} -o {output} {input} && [[ -s {output} ]] || cat {params}/empty.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' | bgzip > {output} ; tabix {output}"

use rule vcf2tsv as vcf2tsv_design_sample with:
	input: rules.split_vcf_design.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.tsv"

use rule split_vcf_design as split_vcf_panel with:
	input: rules.split_vcf_design.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel_unnorm.{{panel}}.vcf.gz")

use rule vcf_normalisation_panel as vcf_normalisation_panel_sample with:
	input: rules.split_vcf_panel.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.vcf.gz"

use rule vcf2tsv as vcf2tsv_panel_sample with:
	input: rules.vcf_normalisation_panel_sample.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv"

use rule merge_vcf_full as merge_vcf_full_noannot with:
	input: expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Full.vcf.gz", aligner=aligner_list, sample=sample_list)
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Full.vcf.gz"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Full.bcftoolsmerge.log"

use rule filter_vcf_design as filter_vcf_design_noannot with:
	input: rules.merge_vcf_full_noannot.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design.vcf.gz"
	params: bed=config['BED_FILE'],
			dummypath=config['DUMMY_FILES']

use rule filter_vcf_design as filter_vcf_panel_noannot with:
	input: rules.filter_vcf_design_noannot.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel_unnorm.{{panel}}.vcf.gz")
	params: bed=lambda wildcards: f"{resultDir}/{wildcards.panel}",
			dummypath=config['DUMMY_FILES']

use rule vcf_normalisation_panel as vcf_normalisation_panel_all_noannot with:
	input: rules.filter_vcf_panel_noannot.output
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
	generate_html_report(resultDict, runName, serviceName, sample_list, f"{serviceName}.template.html" , f"{resultDir}/{serviceName}.{date_time}.report.html")
	copy2(config['TEMPLATE_DIR'] + '/' + serviceName + '.style.css', resultDir)

	# Clear existing output directories & copy results
	for sample in sample_list:
		shell(f"rm -f {outputDir}/{sample}/{serviceName}/* || true")
	shell("rsync -azvh --include={include} --exclude='*' {resultDir}/ {outputDir}")
	for sample in sample_list:
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