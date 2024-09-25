##########################################################################
# Snakemakefile Version:   3.0
# Description:             Snakemake file to run DECoN module (Detection of Exon Copy Number variants)
##########################################################################

# DEV version 0.1 : 10/11/2021
# INT version 0.1 : 17/03/2022
# PROD version 1 : 03/06/2022
# Authoring : Thomas LAVAUX

# PROD version 2 : 14/10/2022 changelog
	# add the possibility to analyse gender separately, extract gender from tag files
	# add a gene list restriction to limit the analysis to a certain list of gene's name (GENE_LIST_RESTRICT)
	# rewrite DECON bed / custom numbering formating
	# remove panel vcf filtering (output was essentially empty)
	# add a chr filter option to remove chromosome (essentially chrY) (CHR_LIST_RESTRICT)
	# add a chrX Male gender analysis only to avoid chrY failure (REMOVE_M_noY)
	# rename failures and results_all from txt to tsv
	# exclude samples list case insensitive
	# copying new analysis per sample in the root sample dir & removing the old ones
	# add vcf2tsv converter ((https://github.com/sigven/vcf2tsvpy) & corresponding rules : each vcf will be convert to tsv
	# correct run path by removing ending '/' if exist
	# keep "Running.txt" file if failed, avoiding multiple analysis launch by the listener

# PROD version 3 : 19/09/2023 changelog
	# AnnotSV version 3.4 (include the vcf converter)
	# more options to process DECON bed, k-merisation
	# R scripts rewrite, to speed up ReadInBam, separate plotting, call CNV with a reference bam list (see changlog in the R scripts for details)
	# re-arrange/simplify/compact code a lot, input files copying is within the rules, f-string, etc.
	# add a draft for an html report using jinja2
	# remove header from tsv converted
################## Import libraries ##################
import os
import glob
import pandas as pd
import json
import csv
import tempfile
import logging
from pypdf import PdfWriter
from shutil import copy2
from datetime import datetime
from itertools import product
from collections import defaultdict
from jinja2 import Environment, FileSystemLoader

################## Configuration file ##################
configfile: "/app/config/snakefile/decon_default.yaml"

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

def deconbedfromdesign(inputbed, outputbed, sepgenecolumn):
	""" Function to extract a DECON bed 4 columns from another bed, assuming that the 4th column contains a gene name with a unique separator """
	df = pd.read_csv(inputbed, sep='\t', names=["Chr", "Start", "End", "Gene"], index_col=False)
	df['Gene'] = df['Gene'].str.split(sepgenecolumn).str[0] # we extract the gene name == string before separator (ex TP73_chr1_3598878_3599018 sep is _chr to get TP73)
	df.to_csv(outputbed, sep='\t', index=False, header=False)

def create_list(txtlistoutput, file_list, ext_list, pattern):
	""" Function to create a txt file from a list of files filtered by extension and pattern """
	with open(txtlistoutput, 'a+') as f:
		f.writelines(f"{files}\n" for files in file_list if any(files.endswith(ext) and pattern in files for ext in ext_list))

def concatenatelist(inputlist, outputfile):
	""" Function to concatenate lines tab separated in a file """
	shell(f"xargs --delimiter='\\n' cat <{inputlist} >> {outputfile}")

def keep_unknown(df, columnname, pattern, prefix):
	""" Function to rename gene """
	for i, gene in enumerate(df[columnname]):
		if gene == pattern:
			df.at[i, columnname] = f"{prefix}{i}"

def extract_tag(tagfile, tag, tagsep, sep):
	""" Function to extract a tag """
	row = open(tagfile, 'r').readline().strip().split(tagsep)
	output_tag = next((items.split(sep)[-1] for items in row if tag in items and sep in items), "")
	return output_tag

def kmerisation(kmerSize, df=None, bedFile=None, kmerBedFile=None):
	""" Bed kmerisation using either DataFrame or file input """
	print(f"{kmerSize} kmerisation in progress")
	
	if df is not None:
		source = df.iterrows()
		source_type = 'DataFrame'
	elif bedFile is not None:
		source = open(bedFile, 'r')
		source_type = 'file'
	else:
		raise ValueError("Either df or bedFile must be provided")

	with open(kmerBedFile, 'w+') as writeKbed:
		if source_type == 'file':
			for line in source:
				if line.startswith('#'):
					continue
				chr, start, end, gene = line.split()[:4]
				start, end = int(start), int(end)
				write_kmers(writeKbed, chr, start, end, gene, kmerSize)
			source.close()
		else:
			for _, row in source:
				chr, start, end, gene = row["Chr"], int(row["Start"]), int(row["End"]), row["Gene"]
				write_kmers(writeKbed, chr, start, end, gene, kmerSize)
	
	print('Kmerisation done')

def write_kmers(writeKbed, chr, start, end, gene, kmerSize):
	""" Helper function to write kmerized intervals """
	diff = end - start
	
	while diff >= kmerSize:
		newEnd = start + kmerSize - 1
		writeKbed.write(f"{chr}\t{start}\t{newEnd}\t{gene}\n")
		start, diff = newEnd + 1, diff - kmerSize
	
	if diff > 0:
		writeKbed.write(f"{chr}\t{start}\t{end}\t{gene}\n")

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

def merge_pdfs(files, output_path):
	pdf_merger = PdfWriter()

	for file in files:
		with open(file, 'rb') as pdf_file:
			pdf_merger.append(pdf_file)

	with open(output_path, 'wb') as output_file:
		pdf_merger.write(output_file)

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

		if gender == '1':
			gender = 'M'
		elif gender == '2':
			gender = 'F'
		else:
			gender = None

		if individual_id in sample_list:
			dictionary.setdefault(individual_id, {})['hpo'] = hpo_list
			dictionary.setdefault(individual_id, {})['gender'] = gender


def intersectbed(input_df, ref_df):
	""" Function to intersect two DataFrames using bedtools and return the result as a DataFrame. """
	
	with tempfile.NamedTemporaryFile(delete=False, mode='w', newline='', suffix='.bed') as input_file, \
		tempfile.NamedTemporaryFile(delete=False, mode='w', newline='', suffix='.bed') as ref_file, \
		tempfile.NamedTemporaryFile(delete=False, mode='w', newline='', suffix='.bed') as output_file:
		input_path = input_file.name
		ref_path = ref_file.name
		output_path = output_file.name
		input_df.to_csv(input_path, sep='\t', index=False, header=False)
		ref_df.to_csv(ref_path, sep='\t', index=False, header=False)
		shell(f"intersectBed -a {input_path} -b {ref_path} -loj > {output_path}")
		result_df = pd.read_csv(output_path, sep='\t', header=None)
		print("Intersection complete.")
		os.remove(input_path)
		os.remove(ref_path)
		os.remove(output_path)
		return result_df

def process_bed_file(bed_file, inputbed_file, bed_process, refseqgene_path=None, transcripts_file=None, unknown_gene=False, gene_list_restrict=None, chr_list_restrict=None, old_bed=False, exon_sep=None, kmer=None, customexon=False, list_genes=None, genes_file=None):
	
	# Case: REGEN - Generate custom exon file and intersect with refseqgene
	if bed_process == 'REGEN':
		df = pd.read_csv(inputbed_file, sep='\t', names=["Chr", "Start", "End"], index_col=False)
		ref_df = pd.read_csv(refseqgene_path, sep='\t')  # Adjust as needed based on refseqgene file format
		intersected_df = intersectbed(df, ref_df)
		intersected_df['Custom.Exon'] = intersected_df['Tosplit'].str.split(',').str[-1].str.replace('exon', '')

		if transcripts_file and os.path.exists(transcripts_file) and os.path.getsize(transcripts_file) != 0:
			NM_list = pd.read_csv(transcripts_file, header=None).squeeze().tolist()
			intersected_df['NM'] = intersected_df['Tosplit'].str.split(',').str[0]
			intersected_df = intersected_df[intersected_df['NM'].isin(NM_list)].drop(columns=['NM'])

		if not unknown_gene:
			intersected_df = intersected_df[intersected_df.Gene != '-1']
		else:
			keep_unknown(intersected_df, 'Gene', '-1', 'Unknown')
		if gene_list_restrict:
			intersected_df = intersected_df[intersected_df['Gene'].isin(gene_list_restrict)]
		if chr_list_restrict:
			intersected_df = intersected_df[~intersected_df['Chr'].isin(chr_list_restrict)]

		df = intersected_df.drop(columns=['4', '5', '6', 'Tosplit', '9']).drop_duplicates()

	# Case: STANDARD - Intersect bed design with panel
	elif bed_process == 'STANDARD':
		if config['LIST_GENES']:
			basepath = os. path. dirname(config['LIST_GENES'])
			cat_panels_bed = f"/tmp/catpanel.bed"
			cat_list_genes = f"/tmp/catlist.genes"
			with open(config['LIST_GENES']) as f:
				cat_panel_list = f.read().splitlines()
				cat_panel_list = [basepath + "/" + s for s in cat_panel_list]
			with open(cat_list_genes, 'w+') as f:
				f.write('\t'.join(cat_panel_list))
			shell(f"xargs --delimiter='\\t' cat < {cat_list_genes} >> {cat_panels_bed}")
		else:
			shell(f"xargs --delimiter='\\t' cat < {config['GENES_FILE']} >> {cat_panels_bed}")
		concatenated_df = pd.read_csv(cat_panels_bed,  sep='\t')
		concatenated_df = concatenated_df.drop_duplicates()
		input_df = pd.read_csv(inputbed_file, sep='\t',header=None)
		intersected_df = intersectbed(input_df, concatenated_df)
		
		if not old_bed:
			headers = ["Chr", "Start", "End", "4", "5", "6", "7", "8", "9", "Gene", "11", "12", "13", "14", "15", "16", "17", "18", "Custom.Exon"]
			columns_to_drop = ['4', '5', '6', '7', '8', '9', '11', '12', '13', '14', '15', '16', '17', '18']
		else:
			headers = ["Chr", "Start", "End", "4", "5", "6", "7", "8", "9", "Gene", "11", "12"]
			columns_to_drop = ['4', '5', '6', '7', '8', '9', '11', '12']

		intersected_df.columns = headers
		df = intersected_df.drop(columns=columns_to_drop)

		if not unknown_gene:
			df = df[df.Gene != '.']
		else:
			keep_unknown(df, 'Gene', '.', 'Unknown')
		if exon_sep:
			df['Gene'] = df['Gene'].str.split(exon_sep).str[0]
		if gene_list_restrict:
			df = df[df['Gene'].isin(gene_list_restrict)]
		if chr_list_restrict:
			df = df[~df['Chr'].isin(chr_list_restrict)]

		df = df.drop_duplicates()

	# Case: NO - Directly convert bed file
	elif bed_process == 'NO':
		df = pd.read_csv(inputbed_file, sep='\t', names=["Chr", "Start", "End", "Gene"], index_col=False)

	# Final Processing: KMER or CustomExon or direct save
	if customexon:
		df.to_csv(bed_file, sep='\t', index=False, header=False)
	elif kmer:
		kmerisation(kmer, df, bed_file)
	else:
		df.to_csv(bed_file, sep='\t', index=False, header=False)

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
	- dictionnary: Dictionary containing sample results.
	- update_dictionnary: Dictionary to be updated with results.
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

# Populate dictionnary
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

# Extract the gender_list from dictionary with the key gender
gender_list = ['A'] + ['XX' if runDict[sample]['gender'] == 'F' else 'XY' for sample in sample_list if sample in runDict and 'gender' in runDict[sample]]
gender_list = list(set(gender_list)) # Removing duplicate

# Find bed file (Design)
config['BED_FILE'] = config['BED_FILE'] or find_item_in_dict(sample_list, config['EXT_INDEX_LIST'], runDict, '.design.bed', '.genes.bed')
if not config['BED_FILE']:
	print('No bed found, DECON cannot continue, exiting')
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

### DECON bed ###
# DECON need a 4 columns bed for the ReadInBams rule, if bed is not set correctly ReadInBams will crash
# Chromosome Start End Gene but no header, bed format (tsv)
# Chr must be formated chr1, Gene column must be a string, no negative numbers
# ex	:	chr17	4511	4889	BCRA1
# DECON can also use a bed with the exons number for an accurate formating of pdfs graph (makeCNVcalls) and the IdentifyFailures metrics
# ex	:	chr17	4511	4889	BCRA1	26
deconbed_file=f"{resultDir}/{serviceName}.{date_time}.bed"
process_bed_file(bed_file=deconbed_file,inputbed_file=config['BED_FILE'],bed_process=config['BED_PROCESS'],refseqgene_path=config['REFSEQGENE_PATH'],transcripts_file=config['TRANSCRIPTS_FILE'],unknown_gene=config['KEEP_UNKNOWN'],gene_list_restrict=config['GENE_LIST_RESTRICT'],chr_list_restrict=config['CHR_LIST_RESTRICT'],old_bed=config['OLD_BED'],exon_sep=config['EXON_SEP'],kmer=config['KMER'],customexon=config['customexon'],list_genes=config['LIST_GENES'],genes_file=config['GENES_FILE'])

# Create file_list of bam by gender with the runDict, depending on the gender_list
# A = all ; XX = Female only ; XY = Male only
# files_list_A contains the full path of all files (bam files for ex)
# Warning the key dictionnary for sexe is A/M/F but the gender_list is A/XY/XX
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
	""" Output a design vcf.gz with the bams list and the bed provided """
	input:
		expand(f"{resultDir}/{{sample}}.{{aligner}}.bam.bai", aligner=aligner_list, sample=sample_list) +
		expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Metrics.tsv", aligner=aligner_list) +
		expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.vcf.gz", aligner=aligner_list, sample=sample_list) +
		[f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Design.tsv"] +
		(expand(f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.plotSuccess", aligner=aligner_list, gender=gender_list) if config['DECON_PLOT'] else []) +
		(expand(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv", aligner=aligner_list, sample=sample_list, panel=panels_list) if config['GENES_FILE'] else []) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel.{{panel}}.tsv", panel=panels_list) if config['GENES_FILE'] else []) +
		(expand(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel.{{panel}}.vcf.gz", aligner=aligner_list, panel=panels_list) if config['GENES_FILE'] else [])


rule help:
	"""
	General help for DECON
	Launch snakemake -s snakefile_decon -c(numberofthreads) --config run=absolutepathoftherundirectory
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

rule ReadInBams:
	""" DECoN calculates FPKM for each exon in each samples BAM file, using a list of BAM files and a BED file """
	input:
		allbamlist=f"{resultDir}/{serviceName}.{date_time}.A.list.txt",
		bai=expand(f"{resultDir}/{{sample}}.{{aligner}}.bam.bai", sample=sample_list, aligner=aligner_list)
	output:
		f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.ReadInBams.RData"
	params:
		refbamlist=config['REF_BAM_LIST'],
		decondir=config['R_SCRIPTS'],
		refgene=config['REFGENEFA_PATH'],
		mcores=config['MCORES']
	log:
		log=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.ReadInBams.log",
		err=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.ReadInBams.err"
	shell:
		"""
		Rscript {params.decondir}/ReadInBams.R --maxcores {params.mcores} --bams {input.allbamlist} --bed {deconbed_file} --fasta {params.refgene} --rdata {output} {params.refbamlist} 1> {log.log} 2> {log.err}
		"""

rule IdentifyFailures:
	""" Evaluate samples and exons for suboptimal CNV calling using the summary RData file """
	input: rules.ReadInBams.output
	params:
		mincorr=config['mincorr'],
		mincov=config['mincov'],
		analysisfailure=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.MetricsFailed",
		decondir=config['R_SCRIPTS']
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Metrics.tsv"
	log:
		log=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.Metrics.log",
		err=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.Metrics.err"
	shell: 
		"""
		Rscript {params.decondir}/IdentifyFailures.R --rdata {input} --mincorr {params.mincorr} --mincov {params.mincov} --tsv {output} 1> {log.log} 2> {log.err} && [[ -s {output} ]] || touch {params.analysisfailure} ; [[ -s {output} ]] || touch {output}
		"""

rule makeCNVcalls:
	""" Call exon CNVs in each sample using reference samples and output correlation results """
	input: rules.ReadInBams.output
	params:
		prob=config['transProb'],
		output=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.Design_results",
		analysisfailure=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.CNVCalledFailed",
		bamlist=f"{resultDir}/{serviceName}.{date_time}.{{gender}}.list.txt",
		chromosome="{gender}",
		refbamlist=config['REF_BAM_LIST'],
		decondir=config['R_SCRIPTS']
	output:
		calltsv=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.Design_results_all.tsv",
		rdata=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.CNVcalls.RData"
	log:
		log=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.makeCNVcalls.log",
		err=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.makeCNVcalls.err"
	shell:
		"""
		Rscript {params.decondir}/makeCNVcalls.R --rdata {input} --samples {params.bamlist} --transProb {params.prob} --chromosome {params.chromosome} {params.refbamlist} --tsv {output.calltsv} --outrdata {output.rdata} 1> {log.log} 2> {log.err} && [[ -s {output.calltsv} ]] || touch {params.analysisfailure} ; [[ -s {output.calltsv} ]] || touch {output.calltsv}
		"""
 
rule merge_makeCNVcalls:
	""" Merge all CNV call results """
	input: expand(f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.Design_results_all.tsv", gender=gender_list, aligner=aligner_list)
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design_uncorr.tsv")
	shell: "cat {input} | sort -u | sort -r >> {output}"

rule fix_makeCNVcalls:
	input: rules.merge_makeCNVcalls.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design.tsv"
	shell: "awk -F'\\t' -v OFS='\\t' 'NR == 1 {{print; next}} {{ $12 = \"chr\" $12; print }}' {input} > {output}"

rule variantconvert:
	""" Convert tsv to a vcf """
	input: rules.fix_makeCNVcalls.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design_convert.vcf")
	params:
		variantconvertconfig=config['VARIANT_CONVERT_CONFIG'],
		dummypath=config['DUMMY_FILES']
	log: f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.variantconvert.log"
	shell: "variantconvert convert -i {input} -o {output} -c {params.variantconvertconfig}/config_decon.json 2> {log}"

rule correct_vcf_svtype:
	input: rules.variantconvert.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design_svtype.vcf")
	shell: 
		""" 
		sed -e 's/SVTYPE=deletion/SVTYPE=DEL/g' -e 's/SVTYPE=duplication/SVTYPE=DUP/g' {input} > {output}
		"""

rule sortvcf:
	""" Sort vcf """
	input: rules.correct_vcf_svtype.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design_sort.vcf")
	params: dummypath=config['DUMMY_FILES']
	shell: "{{ grep \'^#\' {input} && grep -v \'^#\' {input} | sort -k1,1V -k2,2g; }} > {output}"

rule vcf2gz:
	input: rules.sortvcf.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design_svlen.vcf.gz")
	params: config['DUMMY_FILES']
	shell: "bgzip -c {input} > {output} ; tabix {output}"


rule fix_svlen:
	input: rules.vcf2gz.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design_unnorm.vcf.gz")
	params:	config['PY_SCRIPTS']
	shell: " python {params}/multifix_vcf.py -i {input} -o {output} -z ; tabix {output} "

rule vcf_normalization:
	input: rules.fix_svlen.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Design.vcf.gz"
	shell: "bcftools norm -d all -o {output} -Oz {input} ; tabix {output}"

# --force-samples will output a vcf file with all annotations but without a sample name if the sample did not exist
# so we output a dummy vcf.gz
rule split_vcf:
	"""
	Split vcf with bcftools
	-s {sample}: comma-separated list of samples to include
	-Oz: output vcf compressed
	-c1: minimum allele count (INFO/AC) of sites to be printed
	"""
	input: rules.vcf_normalization.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.Design.Noannotation.vcf.gz"
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
	shell: 
		""" 
		sed '1s/[()/\\]/_/g' {input} > {output}
		"""

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

rule variantconvert_annotsv:
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

rule correct_vcf:
	input: rules.variantconvert_annotsv.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_fix.vcf")
	shell:
		""" 
		sed 's/SAMPLENAME/{wildcards.sample}/g' {input} > {output}
		"""

rule sortvcf_sample:
	input: rules.correct_vcf.output
	output: temp(f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design_sort.vcf")
	params: config['DUMMY_FILES']
	shell: " {{ grep \'^#\' {input} && grep -v \'^#\' {input} | sort -k1,1V -k2,2g; }} > {output} && [[ -s {output} ]] || cat {params}/empty.vcf | sed 's/SAMPLENAME/{wildcards.sample}/g' > {output} "

rule fix_vcf:
	input: rules.sortvcf_sample.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Design.vcf.gz"
	params:	config['PY_SCRIPTS']
	shell: " python {params}/multifix_vcf.py -i {input} -o {output} -z ; tabix {output} "

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

use rule vcf_normalization as vcf_normalisation_annot with:
	input: rules.filter_vcf.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel.{{panel}}.vcf.gz"

use rule vcf2tsv as vcf2tsv_panel with:
	input: rules.vcf_normalisation_annot.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.AnnotSV.Panel.{{panel}}.tsv"

use rule split_vcf as split_vcf_panel with:
	input: rules.fix_vcf.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.vcf.gz"

use rule vcf2tsv as vcf2tsv_panel_sample with:
	input: rules.split_vcf_panel.output
	output: f"{resultDir}/{{sample}}/{serviceName}/{{sample}}_{date_time}_{serviceName}/{serviceName}.{date_time}.{{sample}}.{{aligner}}.AnnotSV.Panel.{{panel}}.tsv"

use rule filter_vcf as filter_vcf_panel with:
	input: rules.vcf_normalization.output
	output: temp(f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel_unnorm.{{panel}}.vcf.gz")
	log: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel.{{panel}}.bedtoolsfilter.log"

use rule vcf_normalization as vcf_normalisation_all with:
	input: rules.filter_vcf_panel.output
	output: f"{resultDir}/{serviceName}.{date_time}.allsamples.{{aligner}}.Panel.{{panel}}.vcf.gz"


rule plot:
	""" This step inputs the CNVcall.Rdata file to output the plot files in pdf format """
	input:
		rules.makeCNVcalls.output.rdata
	params:
		folder=f"{resultDir}/pdfs/",
		decondir=config['R_SCRIPTS']
	output:
		f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.plotSuccess"
	log:
		log=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.plots.log",
		err=f"{resultDir}/{serviceName}.{date_time}.{{aligner}}.{{gender}}.plots.err"
	shell:
		"""
		mkdir -p {params.folder} &&
		Rscript {params.decondir}/DECONplot.R --rdata {input} --out {params.folder} --prefix {date_time} 1> {log.log} 2> {log.err} &&
		touch {output}
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
	
	# Move pdfs into sample's directories
	pdf_list = os.listdir(f"{resultDir}/pdfs")
	pdf_path_list = [f"{resultDir}/pdfs/{pdf}" for pdf in pdf_list]	
	
	output_pdf_list = []
	for sample in sample_list:
		matching_pdfs = [pdf for pdf in pdf_path_list if sample in pdf]
		
		if matching_pdfs:
			output_pdf = f"{resultDir}/pdfs/{serviceName}.{date_time}.{sample}.pdf"
			output_pdf_list.append(os.path.basename(output_pdf))
			merge_pdfs(matching_pdfs, output_pdf)

	pdf_list.extend(output_pdf_list)

	for sample in sample_list:
		for pdf in pdf_list:
			if sample in pdf:
				shell(f"mv {resultDir}/pdfs/{pdf} {resultDir}/{sample}/{serviceName}/{sample}_{date_time}_{serviceName}/")
	shell(f"rm -rf {resultDir}/pdfs")

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
	shell(f"rm -rf {resultDir}/pdfs")
	shell("rsync -azvh --include={include_log} --exclude='*' {resultDir}/ {outputDir}")