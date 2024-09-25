##########################################################################
# Snakemakefile Version:   1.0
# Description:             Snakemake file to run fusions module
##########################################################################

# DEV version 1.0 : 23/11/2022
# Authoring : Thomas LAVAUX

################## Context ##################
# launch snakemake -s  snakefile_fusions -c(numberofthreads) --use-conda --config run=absolutepathoftherundirectory without / at the end of the path
# to launch the snakemake file, use --config to replace variables that must be properly set for the pipeline to work ie run path directory
# every variable defined in the yaml file can be change
# separate multiple variable with a space (ex  --config run=runname var1=0.05 var2=12)
# also use option --configfile another.yaml to replace and merge existing config.yaml file variables

# use -p to display shell commands
# use --lt to display docstrings of rules

# input file = bam files indexed and .arriba.fusions.tsv from arriba or STAR Fusion
# output file = arriba pdf report
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
configfile: "/app/config/snakefile/fusions_default.yaml"

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

def generate_html_report(result_dict, run_name, service_name, sample_list, template_name, output_file='report.html', sample_list_added=None):
	env = Environment(loader=FileSystemLoader(config['TEMPLATE_DIR']))
	template = env.get_template(template_name)

	rendered_html = template.render(
		runDict=result_dict,
		runName=run_name,
		serviceName=service_name,
		sample_list=sample_list,
		sample_list_added=sample_list_added
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
				update_dictinnary[sample].pop(key, None)  # Remove the key if it exists

### END OF FUNCTIONS ###
serviceName = config['serviceName']
fusion_list=config['FUSION_CALLER']
date_time = config['DATE_TIME'] if config['DATE_TIME'] else datetime.now().strftime("%Y%m%d-%H%M%S")
runName = os.path.basename(os.path.normpath(config['run']))
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

runDict = populate_dictionary(sample_list, config['EXT_INDEX_LIST'], filtered_files, None, None)
print(dict(runDict))

# Find bed file (Design)
config['BED_FILE'] = config['BED_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.design.bed', '.genes.bed')
# Find genes file (Panel); we can't use .genes files because .list.genes and .genes are not distinctable from the indexing we made
config['GENES_FILE'] = config['GENES_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.genes.bed', '.list.genes')
# Find list.genes files 
config['LIST_GENES'] = config['LIST_GENES'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.list.genes', '.list.transcripts')
# Find transcripts files (NM)
config['TRANSCRIPTS_FILE'] = config['TRANSCRIPTS_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.transcripts', '.list.transcripts')

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

################################################## RULES ##################################################
# check the number of sample for copy or merge vcf rule
sample_count = len(sample_list) 

rule all:
	"""
	Rule will create a vcf for arriba and starfusion
	"""
	input:
		expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{fusion}}.vcf.gz", sample=sample_list, fusion=fusion_list),
		expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{fusion}}.tsv", fusion=fusion_list),
		expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{{sample}}.STAR.bam.bai", sample=sample_list),
		expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{fusion}}.pdf", sample=sample_list, fusion=fusion_list)

rule help:
	"""
	General help for drawfusions module
	Launch snakemake -s  snakefile_fusions -c(numberofthreads) --config DATA_DIR=absolutepathoftherundirectory (default is data) without / at the end of the path
	To launch the snakemake file, use --config to replace variables that must be properly set for the pipeline to work ie run path directory
	Every variable defined in the yaml file can be change
	Separate multiple variable with a space (ex  --config DATA_DIR=runname transProb=0.05 var1=0.05 var2=12)
	Also use option --configfile another.yaml to replace and merge existing config.yaml file variables
	Use -p to display shell commands
	Use --lt to display docstrings of rules
	Input file = bam files (if bai is needed, it will be generate), tsv output fusion list from arriba
	Output file = pdf report with arriba graphs
	"""


rule copy_fastq:
	output:
		fastqR1=temp(f"{resultDir}/{{sample}}.R1.fastq.gz"),
		fastqR2=temp(f"{resultDir}/{{sample}}.R2.fastq.gz")
	params:
		process = config['PROCESS_CMD'],
		download_link1 = lambda wildcards: runDict[wildcards.sample]['.R1.fastq.gz'],
		download_link2 = lambda wildcards: runDict[wildcards.sample]['.R2.fastq.gz']
	shell: "[ \"{params.process}\" = \"ln\" ] && ln -sfn {params.download_link1} {output.fastqR1} && ln -sfn {params.download_link2} {output.fastqR2} || rsync -azvh {params.download_link1} {output.fastqR1} && rsync -azvh {params.download_link2} {output.fastqR2}"


# indexing genome for STAR with the command :
# STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /STARK/databases/STAR/ --genomeFastaFiles /STARK/databases/STAR/hg19.fa --sjdbGTFfile /STARK/databases/gtf/current/gencode.v19.annotation.gtf_withproteinids --sjdbOverhang 99 --outFileNamePrefix /STARK/databases/STAR/log/ --genomeChrBinNbits 15

# STAR produces multiple output files. All files have standard name, however, you can change the file
# prefixes using --outFileNamePrefix /path/to/output/dir/prefix. By default, this parameter is
# ./, i.e. all output files are written in the current directory
# name.bam + Aligned.toTranscriptome.out.bam or Aligned.sortedByCoord.out.bam
# --limitBAMsortRAM 32000000000 is minimum
# output will be Sample.Chimeric.out.junction
# --genomeLoad LoadAndRemove or --genomeLoad NoSharedMemory if 2 pass
rule STAR:
	input:
		fastqR1=f"{resultDir}/{{sample}}.R1.fastq.gz", 
		fastqR2=f"{resultDir}/{{sample}}.R2.fastq.gz",
	output:
		bam = f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{{sample}}.STAR.bam",
		junction = f"{resultDir}/{{sample}}/{{sample}}.Chimeric.out.junction"
	threads: workflow.cores
	params:
		genomepath_indexed=config['GENOME_PATH_INDEXED'],
		RG="ID:1 PL:ILLUMINA PU:PU LB:001 SM:{sample}",
		tmpdir=f"{resultDir}/tmp/",
		juncfile=config['JUNCTION_FILE'],
		prefix=f"{resultDir}/{{sample}}/{{sample}}."
	shell:
		"""
		STAR --genomeDir {params.genomepath_indexed} --runThreadN {threads} --readFilesIn {input.fastqR1} {input.fastqR2} --genomeLoad NoSharedMemory --readFilesCommand zcat --outFileNamePrefix {params.prefix} --outSAMattrRGline {params.RG} --outSAMtype BAM SortedByCoordinate --chimOutJunctionFormat 1 --outSAMunmapped Within --outBAMcompression 6 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType {params.juncfile} HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --outTmpDir {params.tmpdir} && mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam} ; rm -rf {params.prefix}_STARgenome/ ; rm -rf {params.prefix}_STARpass1/
		"""

rule indexing:
	""" Indexing bam files with samtools """
	input: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{{sample}}.STAR.bam"
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{{sample}}.STAR.bam.bai"
	threads: workflow.cores
	shell: "samtools index -b -@ {threads} {input} {output}"

# -f will disable filters / -f known_fusions to reduce sensibility / -f isoforms to reduce false positive
# -S MIN_SUPPORTING_READS : The filter min_support discards all fusions with fewer than this many supporting reads (split reads and discordant mates combined). Default: 2
# -z MIN_ITD_ALLELE_FRACTION : Required fraction of supporting reads to report an internal tandem duplication. Default: 0.07
# -Z MIN_ITD_SUPPORTING_READS : Required absolute number of supporting reads to report an internal tandem duplication. Default: 10
# -k File containing known/recurrent fusions
# -t Tab-separated file containing fusions to annotate with tags in the 'tags' column
rule Arriba:
	input: rules.STAR.output.bam
	output:	f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.arriba.tsv"
	params:
		refgenome = config['REFGENOMEFA_PATH'],
		refgtfgencode = config['REFGTFGENCODE_PATH'],
		refprot = config['REFPROTDOMAIN_PATH'],
		blacklist = config['BLACKLIST_PATH'],
		reffusion = config['REFFUSION_PATH'],
		fusion = config ['REFFUSION_PATH']
	shell:
		"""
		arriba -x {input} -o {output} -a {params.refgenome} -g {params.refgtfgencode} -b {params.blacklist} -k {params.reffusion} -t {params.fusion} -p {params.refprot} -f isoforms
		"""

# Settings for STAR used by Starfusion
# STAR  --chimSegmentMin 12  --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 
rule STARFusion:
	input: rules.STAR.output.junction
	output: f"{resultDir}/{{sample}}/tmp/star-fusion.fusion_predictions.tsv"
	params:
		libgencode=config['CTAT_LIB_GENCODE'],
		fusiondir=f"{resultDir}/{{sample}}/tmp/"
	conda: "starfusion"
	shell:
		"""
		STAR-Fusion --genome_lib_dir {params.libgencode} -J {input} --output_dir {params.fusiondir}
		"""

rule rename_tsv:
	input: rules.STARFusion.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.starfusion.tsv"
	shell:
		"""
		mv {input} {output} ; rm -rf {resultDir}/{{sample}}/tmp/
		"""

rule variantconvert:
	input: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{fusion}}.tsv"
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{fusion}}.vcf"
	params: lambda wildcards: "/app/config/fileconversion/" + wildcards.fusion + ".json"
	shell: "variantconvert convert -i {input} -o {output} -c {params}"

rule vcf2gz:
	"""	Compress vcf with bgzip	"""
	input: rules.variantconvert.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{fusion}}.vcf.gz"
	shell: "bgzip -c {input} > {output} ; tabix {output}"

rule mergevcf:
	"""	Copy or merge vcfs with bcftools """
	input:
		lambda wildcards: expand(
			f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{wildcards.fusion}.vcf.gz", 
			sample=sample_list, 
			fusion=[wildcards.fusion])
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{fusion}}.vcf.gz"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{fusion}}.bcftoolsmerge.log"
	shell: "if [ {sample_count} -eq 1 ]; then cp {input} {output}; else bcftools merge {input} -O z -o {output} 2> {log}; fi; tabix {output}"

rule vcf2tsv:
	""" vcf to tsv conversion """
	input: rules.mergevcf.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{fusion}}.tsv"
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{fusion}}.vcf2tsv_converter.log"
	shell: "vcf2tsvpy --keep_rejected_calls --input_vcf {input} --out_tsv {output}.tmp && cat {output}.tmp | grep -v '^#' > {output}"

rule DrawR:
	input:
		bam=f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{{sample}}.STAR.bam",
		bai=f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{{sample}}.STAR.bam.bai",
		fusiontsv=f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{fusion}}.tsv"
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{fusion}}.pdf"
	params:
		refgtfgencode=config['REFGTFGENCODE_PATH'],
		refprot=config['REFPROTDOMAIN_PATH'],
		cytoband=config['REFCYTOBAND_PATH'],
		arriba_scripts=config['ARRIBA_SCRIPTS']
	shell:
		"""
		Rscript {params.arriba_scripts}/draw_fusions.R --annotation={params.refgtfgencode} --fusions={input.fusiontsv} --output={output} --alignments={input.bam} --cytobands={params.cytoband} --proteinDomains={params.refprot}
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