# -*- coding: utf-8 -*-
"""
@Goal: Expand Celine Besnard's script with infinite conversion abilities
@Author: Samuel Nicaise
@Date: 23/11/2021

Prerequisites: pandas, pyfasta

Usage examples:
#TSV (Decon) to VCF
python /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/file_converter.py -i /home1/BAS/nicaises/Tests/deconconverter/200514_NB551027_0724_AHTWHHAFXY.DECON_results_all.txt -o vcf_from_decon.vcf -fi tsv -fo vcf -c /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/fileconversion/config_decon.json

#AnnotSV3 to VCF
python /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/file_converter.py -i /home1/BAS/nicaises/Tests/deconconverter/DECoN.20211207-183955_results_both.tsv -o /home1/BAS/nicaises/Tests/deconconverter/decon__annotsv3.vcf -fi annotsv -fo vcf -c /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/fileconversion/config_annotsv3.json

----------------------------
Configfile guidelines (JSON)
----------------------------
1)[GENERAL] has 3 important fields
	#source format name: will show up in VCF meta fields
	#skip_rows: how many rows to skip before we reach indexes.
	This script cannot handle a tsv with unnamed columns (beds are fine)
	#unique_variant_id: useful in multisample files. List the
	columns that are needed to uniquely identify a variant.
2) [VCF_COLUMNS] describe the columns that will go in your VCF
	key: column name in VCF ; value: column name in source format
3) [COLUMNS_DESCRIPTION] describe the tsv columns
	Type and Description fields will be used in the VCF header
4) read HelperFunctions docstring

"""
from __future__ import division
from __future__ import print_function

import argparse
import json
import os
import sys
import time

from functools import lru_cache
from os.path import join as osj
import pandas as pd
#from pyfaidx import Fasta # bug so revert to pyfasta
from pyfasta import Fasta

class Converter:
	"""
	Main interface. Fetch the proper converter depending on input and output formats.
	"""
	def __init__(self, source_format, dest_format, config):
		factory = ConverterFactory()
		self.converter = factory.get_converter(source_format, dest_format)
		self.config = config
	def convert(self, file, output_path):
		self.converter.convert(file, self.config, output_path)

class ConverterFactory:
	"""
	Factory pattern implementation
	To add new converters, use register_converter() or add them directly in __init__()
	"""
	def __init__(self):
		self._converters = {}
		self._converters["varank>vcf"] = VcfFromVarank
		self._converters["annotsv>vcf"] = VcfFromAnnotsv
		self._converters["bed>vcf"] = VcfFromBed
		self._converters["tsv>vcf"] = VcfFromTsv

	def register_converter(self, source_format, dest_format, converter):
		self._converters[source_format + ">" + dest_format] = converter

	def get_converter(self, source_format, dest_format):
		converter = self._converters.get(source_format + ">" + dest_format)
		if not converter:
			raise ValueError("Unknown converter: " + source_format+">"+dest_format)
		return converter()

class VcfFromVarank:
	"""
	TODO: update vcffromvarank.py to fit the interface and import it instead
	"""
	def convert(self, tsv, config, output_path):
		print("Converting to vcf from Varank")
		raise ValueError("Not implemented yet")

class VcfFromAnnotsv:
	"""
	Specificities compared to TSV:
	- vcf-like INFO field  in addition to other annotation columns
	- vcf-like FORMAT and <sample> fields
	- full/split annotations. Each variant can have one "full" and 
	zero to many "split" annotations which result in additional lines in the file
	Very hard to deal with this with generic code --> it gets its own converter
	"""
	def _build_input_dataframe(self):
		df = pd.read_csv(self.filepath,
							skiprows=self.config["GENERAL"]["skip_rows"],
							sep = '\t',
							low_memory = False)
		df.sort_values([self.config["VCF_COLUMNS"]["#CHROM"],
							self.config["VCF_COLUMNS"]["POS"]],
							inplace = True)
		df.reset_index(drop = True, inplace = True)
		df.fillna('.', inplace = True)
		df = df.astype(str)
		# print(df)
		return df

	def _get_sample_list(self):
		samples_col = self.input_df[self.config["VCF_COLUMNS"]["SAMPLE"]]
		sample_list = []
		for cell in samples_col:
			if "," in cell:
				for sample in cell.split(","):
					sample_list.append(sample)
			else:
				sample_list.append(cell)
		# print(samples_col)
		sample_list = list(set(sample_list))
		# print("sample_list:", sample_list)
		if not set(sample_list).issubset(self.input_df.columns):
			raise ValueError("All samples in '" + samples_col + "' column are expected to "\
					"have their own column in the input AnnotSV file")
		return sample_list

	def _build_input_annot_df(self):
		"""
		remove FORMAT, Samples_ID, and each <sample> column
		TODO: remove vcf base cols ; INFO field
		"""
		columns_to_drop = [v for v in self.sample_list]
		columns_to_drop += [v for v in self.main_vcf_cols]
		columns_to_drop.append(self.config["VCF_COLUMNS"]["SAMPLE"])
		columns_to_drop.append(self.config["VCF_COLUMNS"]["FORMAT"])
		columns_to_drop.append(self.config["VCF_COLUMNS"]["INFO"]["INFO"])
		df = self.input_df.drop(columns_to_drop, axis=1)
		df = df.replace(';', ',', regex=True) #any ';' in annots will ruin the vcf INFO field
		return df

	def _merge_full_and_split(self, df):
		"""
		input: df of a single annotSV_ID ; containing only annotations (no sample/FORMAT data)
		it can contain full and/or split annotations
		
		returns a single line dataframe with all annotations merged properly
		"""
		if self.config["GENERAL"]["mode"] != "full&split":
			raise ValueError("Unexpected value in json config['GENERAL']['mode']: "\
							"only 'full&split' mode is implemented yet.")

		annots = {}
		dfs = {}
		for type, df_type in df.groupby(self.config["VCF_COLUMNS"]["INFO"]["Annotation_mode"]):
			if type not in ("full", "split"):
				raise ValueError("Annotation type is assumed to be only 'full' or 'split'")
			dfs[type] = df_type

		#deal with full
		if "full" not in dfs.keys():
			#still need to init columns
			if "split" not in dfs.keys():
				print("[WARNING] Input does not include AnnotSV's 'Annotation_mode' column. This is "\
						"necessary to know how to deal with annotations. The INFO field will be empty.")
				return {}
			for ann in dfs["split"].columns:
				annots[ann] = "."
		else:
			if len(dfs["full"].index) > 1 :
				raise ValueError("Each variant is assumed to only have one single line of 'full' annotation")
			for ann in dfs["full"].columns:
				annots[ann] = dfs["full"].loc[df.index[0], ann]

		#deal with split
		if "split" not in dfs.keys():
			return annots
		for ann in dfs["split"].columns:
			if ann == self.config["VCF_COLUMNS"]["INFO"]["Annotation_mode"]:
				annots[ann] = self.config["GENERAL"]["mode"]
				continue
			if annots[ann] != ".":
				continue #'full' annot is always prioritized
			annots[ann] = ",".join(dfs["split"][ann].tolist())

		#remove empty annots
		annots = {k:v for k,v in annots.items() if v != "."}
		return annots

	def _build_info_dic(self):
		"""
		Output: dictionary with key: annotsv_ID ; value: a key-value dictionary of all annotations
		This will be used to write the INFO field
		"""
		input_annot_df = self._build_input_annot_df()
		# print(input_annot_df)
		annots_dic = {}
		id_col = self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"]
		for variant_id, df_variant in input_annot_df.groupby(id_col):
			merged_annots = self._merge_full_and_split(df_variant)
			annots_dic[variant_id] = merged_annots
		return annots_dic

	#TODO: merge this with the other create_vcf_header method if possible
	#Making this method static is an attempt at making it possible to kick it out of the class
	@staticmethod 
	def _create_vcf_header(input_path, config, sample_list, input_df, info_keys):
		header = []
		header.append("##fileformat=VCFv4.3")
		header.append("##fileDate=%s" % time.strftime("%d/%m/%Y"))
		header.append("##source=" + config["GENERAL"]["origin"])
		header.append("##InputFile=%s" % os.path.abspath(input_path))

		if config["VCF_COLUMNS"]["FILTER"] in input_df.columns:
			for filter in set(input_df[config["VCF_COLUMNS"]["FILTER"]].to_list()):
				header.append('##FILTER=<ID=' + str(filter) + ',Description=".">')
		else:
			header.append('##FILTER=<ID=PASS,Description="Passed filter">')

		#identify existing values in header_dic...
		header_dic = {}
		header_dic["REF"] = set(input_df[config["VCF_COLUMNS"]["REF"]].to_list())
		header_dic["ALT"] = set(input_df[config["VCF_COLUMNS"]["ALT"]].to_list())
		#... then for each of them, check if a description was given in the config
		for section in header_dic:
			for key in header_dic[section]:
				if key in config["COLUMNS_DESCRIPTION"][section]:
					info_config = config["COLUMNS_DESCRIPTION"][section][key]
					header.append('##' + section + '=<ID=' + key + ',Description="' 
								+ info_config["Description"] + '">')
				else:
					header.append('##' + section + '=<ID=' + key 
									+ ',Description="Imported from ' 
									+ config["GENERAL"]["origin"] + '">')

		#same as before, but for header elements that also have a type...
		header_dic = {}
		header_dic["INFO"] = info_keys
		header_dic["FORMAT"] = set()
		for format_field in input_df[config["VCF_COLUMNS"]["FORMAT"]].to_list():
			for format in format_field.split(":"):
				header_dic["FORMAT"].add(format)
		#... then for each of them, check if a description was given in the config
		for section in header_dic:
			for key in header_dic[section]:
				if key in config["COLUMNS_DESCRIPTION"][section]:
					info_config = config["COLUMNS_DESCRIPTION"][section][key]
					header.append('##' + section + '=<ID=' + key + ',Number=.,Type='
								+ info_config["Type"] + ',Description="' 
								+ info_config["Description"] + '">')
				else:
					header.append('##' + section + '=<ID=' + key 
									+ ',Number=.,Type=String,Description="Imported from ' 
									+ config["GENERAL"]["origin"] + '">')

		header += config["GENOME"]["vcf_header"]
		header.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT",
									"QUAL", "FILTER", "INFO", "FORMAT"] + sample_list))
		return header

	def _get_main_vcf_cols(self):
		cols = [self.config["VCF_COLUMNS"]["#CHROM"]]
		cols.append(self.config["VCF_COLUMNS"]["POS"])
		cols.append(self.config["VCF_COLUMNS"]["ID"])
		cols.append(self.config["VCF_COLUMNS"]["REF"])
		cols.append(self.config["VCF_COLUMNS"]["ALT"])
		cols.append(self.config["VCF_COLUMNS"]["QUAL"])
		cols.append(self.config["VCF_COLUMNS"]["FILTER"])
		return cols

	def convert(self, tsv, json_config, output_path):
		"""
		Creates and fill the output file.
		
		For each annotSV_ID ; fetch all related lines of annotations in key value dics.
		A function identifies and merges the annotations. 
		then we build a dictionary with 1 dictionary per annotSV_ID containing all the annotations
		This is then used to make the header and fill the INFO field. 
		
		Note: the "INFO" field from annotSV is discarded for now, 
		because it only contains Decon annotations and they're useless.
		TODO: make an option to keep the "INFO" field in the annotations dictionary
		"""
		print("Converting to vcf from tsv using config: " + json_config)

		self.filepath = tsv
		with open(json_config, "r") as f:
			self.config = json.load(f)
		helper = HelperFunctions(self.config)

		self.input_df = self._build_input_dataframe()
		self.sample_list = self._get_sample_list()
		self.main_vcf_cols = self._get_main_vcf_cols()
		
		info_dic = self._build_info_dic()
		info_keys = set()
		for id, dic in info_dic.items():
			for k in dic:
				info_keys.add(k)

		# create the vcf
		with open(output_path, "w") as vcf:
			vcf_header = self._create_vcf_header(tsv, self.config, self.sample_list, self.input_df, info_keys)
			for l in vcf_header:
				vcf.write(l+"\n")

			id_col = self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"]
			for variant_id, df_variant in self.input_df.groupby(id_col):
				main_cols = "\t".join(df_variant[self.main_vcf_cols].iloc[0].to_list())
				vcf.write(main_cols + "\t")
				vcf.write(";".join([k+"="+v for k, v in info_dic[variant_id].items()]) + "\t")
				sample_cols = "\t".join(df_variant[[self.config["VCF_COLUMNS"]["FORMAT"]] 
													+ self.sample_list].iloc[0].to_list())
				vcf.write(sample_cols + "\t")
				vcf.write("\n")

class VcfFromBed:
	"""
	TODO: insert related celine's code here
	"""
	def convert(self, tsv, config, output_path):
		print("Converting to vcf from bed")
		raise ValueError("Not implemented yet")

class VcfFromTsv:
	def _init_dataframe(self):
		self.df = pd.read_csv(self.filepath,
							skiprows=self.config["GENERAL"]["skip_rows"],
							sep = '\t',
							low_memory = False)
		self.df.sort_values([self.config["VCF_COLUMNS"]["#CHROM"],
							self.config["VCF_COLUMNS"]["POS"]],
							inplace = True)
		self.df.reset_index(drop = True, inplace = True)
		self.df.fillna('.', inplace = True)
		print(self.df)
		self.df["__!UNIQUE_VARIANT_ID!__"] = self.df.apply(lambda row: self._get_unique_variant_id(row), axis=1)
		if self.config["VCF_COLUMNS"]["SAMPLE"] != "":
			self.df[self.config["VCF_COLUMNS"]["SAMPLE"]] = self.df.apply(lambda row: self._bwamem_name_bugfix(row), axis=1)
		print(self.df)

	def _get_sample_list(self):
		#is the file multisample?
		if self.config["VCF_COLUMNS"]["SAMPLE"] != "":
			sample_list = []
			for sample in self.df[self.config["VCF_COLUMNS"]["SAMPLE"]].unique():
				sample_list.append(sample)
			return sample_list
		else:
			return [os.path.basename(self.output_path)]

	def _bwamem_name_bugfix(self, row):
		"""remove .bwamem from the end of sample names if needed"""
		name = row[self.config["VCF_COLUMNS"]["SAMPLE"]]
		if name.endswith(".bwamem") and name != ".bwamem":
			return name[0:-7]
		else:
			return name

	def _get_unique_variant_id(self, row):
		id = []
		for col in self.config["GENERAL"]["unique_variant_id"]:
			id.append(str(row[col]))
		return "_".join(id)

	def _get_unique_id_to_index_list(self, data):
		id_dic = {}
		for k, v in data["__!UNIQUE_VARIANT_ID!__"].items():
			if v not in id_dic:
				id_dic[v] = [k]
			else:
				id_dic[v].append(k)
		return id_dic

	def convert(self, tsv, json_config, output_path):
		print("Converting to vcf from annotSV using config: " + json_config)

		self.filepath = tsv
		self.output_path = output_path
		with open(json_config, "r") as f:
			self.config = json.load(f)
		self._init_dataframe()
		sample_list = self._get_sample_list()
		helper = HelperFunctions(self.config)

		with open(output_path, "w") as vcf:
			vcf_header = create_vcf_header(tsv, self.config, sample_list)
			for l in vcf_header:
				vcf.write(l+"\n")

			data = self.df.astype(str).to_dict()
			#In Decon (and maybe others), TSV are given as a list of variant-sample associations
			#so the same variant can be on multiple TSV lines
			# __!UNIQUE_VARIANT_ID!__ allows to identify variants and only add them to the VCF once
			already_seen_variants = set()
			unique_id_to_index_list = self._get_unique_id_to_index_list(data)
			for i in range(len(data[self.config["VCF_COLUMNS"]["#CHROM"]])):
				if data["__!UNIQUE_VARIANT_ID!__"][i] in already_seen_variants:
					continue

				line = ""

				for vcf_col in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL"]:
					col = self.config["VCF_COLUMNS"][vcf_col]
					if is_helper_func(col):
						#col[1] is a function name, col[2] its list of args
						#the function named in col[1] has to be callable from this module
						func = helper.get(col[1])
						args = [data[c][i] for c in col[2:]]
						line += func(*args) + "\t"
					elif col == "":
						line += ".\t"
					else:
						line += data[col][i] + "\t"

				#Cutting-edge FILTER implementation
				line += "PASS\t"

				info_field = []
				for vcf_col, tsv_col in self.config["VCF_COLUMNS"]["INFO"].items():
					if is_helper_func(tsv_col):
						func = helper.get(tsv_col[1])
						args = [data[c][i] for c in tsv_col[2:]]
						s = vcf_col + "=" + func(*args)
					else:
						s = vcf_col + "=" + data[tsv_col][i]
					info_field.append(clean_string(s))
				line += ";".join(info_field) + "\t"

				vcf_format_fields = []
				tsv_format_fields = []
				for vcf_col, tsv_col in self.config["VCF_COLUMNS"]["FORMAT"].items():
					vcf_format_fields.append(vcf_col)
					tsv_format_fields.append(tsv_col)
				line += ":".join(vcf_format_fields) + "\t"

				#monosample input
				if len(sample_list) == 1:
					sample_field = []
					for key, val in self.config["VCF_COLUMNS"]["FORMAT"].items():
						if key == "GT" and val == "":
							sample_field.append("0/1")
							continue
						sample_field.append(data[val][index])
					line += ":".join(sample_field)

				#multisample input
				else:
					sample_field_dic = {}
					#If the variant exists in other lines in the source file, fetch their sample data now
					for index in unique_id_to_index_list[data["__!UNIQUE_VARIANT_ID!__"][i]]:
						sample_field = []
						for key, val in self.config["VCF_COLUMNS"]["FORMAT"].items():
							if key == "GT" and val == "":
								sample_field.append("0/1")
								continue
							sample_field.append(data[val][index])
						sample_field_dic[data[self.config["VCF_COLUMNS"]["SAMPLE"]][index]] = ":".join(sample_field)

					for sample in sample_list:
						if sample in sample_field_dic:
							line += sample_field_dic[sample] + "\t"
						else:
							if "GT" in self.config["VCF_COLUMNS"]["FORMAT"]:
								empty = "./.:" + ":".join(["." for i
															in range(len(self.config["VCF_COLUMNS"]["FORMAT"]) - 1)
														])
							else:
								empty = ":".join(["." for i
													in range(len(self.config["VCF_COLUMNS"]["FORMAT"]) - 1)
												])
							line += empty + "\t"
					already_seen_variants.add(data["__!UNIQUE_VARIANT_ID!__"][i])
					line.rstrip('\t')

				vcf.write(line+"\n")

class HelperFunctions:
	"""
	For when you can't just convert columns by changing column names
	Steps needed:
	- define the helper function
	- update self.dispatcher so this class can redirect to the proper function
	- tell the config file you want to use a HELPER_FUNCTION with the following pattern:
		[HELPER_FUNCTION, <name of your function>, <arg1>, <arg2>, ..., <argN>]

	Example: I need a LENGTH value in my destination format,
	but my source file only has START and END columns.
	You would need:
	# somewhere in the module
		def get_length(start, end):
			return str(end - start)
	#in this class __init__():
		self.dispatcher["get_length_from_special_format"]: get_length
	# in the JSON configfile
		LENGTH: [HELPER_FUNCTION, "get_length_from_special_format", START, END]
	"""
	def __init__(self, config):
		self.config = config
		self.dispatcher = {
			"get_ref_from_decon": self.get_ref_from_decon,
			"get_alt_from_decon": self.get_alt_from_decon,
			"get_svlen_from_decon": self.get_svlen_from_decon,
			"get_info_from_annotsv": self.get_info_from_annotsv
		}

	def get(self, func_name):
		return self.dispatcher[func_name]

	def get_ref_from_decon(self, chr, start):
		f = get_genome(self.config["GENOME"]["path"])
		return f[chr][int(start) -1]

	@staticmethod
	def get_alt_from_decon(cnv_type_field):
		if cnv_type_field == "deletion":
			return "<DEL>"
		if cnv_type_field == "duplication":
			return "<DUP>"
		raise ValueError("Unexpected CNV.type value:" + str(cnv_type_field))

	@staticmethod
	def get_svlen_from_decon(start, end):
		return str(int(end) - int(start))

	@staticmethod
	def get_info_from_annotsv(info):
		"""
		only used in attempts to convert annotsv files 
		as if they were generic TSV (not recommended)
		"""
		return "."

def is_helper_func(arg):
	if isinstance(arg, list):
		if arg[0] == "HELPER_FUNCTION":
			return True
		else:
			raise ValueError("This config file value should be \
							a String or a HELPER_FUNCTION pattern:" + arg)
	return False

@lru_cache
def get_genome(fasta_path):
	return Fasta(fasta_path)

def clean_string(s):
	"""
	replace characters that will crash bcftools and/or cutevariant
	those in particular come from Varank files
	"""
	replace = {";": ",",
				"“": '"',
				"”": '"',
				"‘": "'",
				"’": "'"}
	for k, v in replace.items():
		s = s.replace(k, v)
	return s

def create_vcf_header(input_path, config, sample_list):
	header = []
	header.append("##fileformat=VCFv4.3")
	header.append("##fileDate=%s" % time.strftime("%d/%m/%Y"))
	header.append("##source=" + config["GENERAL"]["origin"])
	header.append("##InputFile=%s" % os.path.abspath(input_path))

	#TODO: FILTER is not present in any tool implemented yet
	#so all variants are set to PASS
	if config["VCF_COLUMNS"]["FILTER"] != "":
		raise ValueError('Filters are not implemented yet. '\
						'Leave config["COLUMNS_DESCRIPTION"]["FILTER"] empty '\
						'or whip the developer until he does it.'\
						'If you are trying to convert an annotSV file, '\
						'use "annotsv" in the input file format argument')
	header.append('##FILTER=<ID=PASS,Description="Passed filter">')

	if "ALT" in config["COLUMNS_DESCRIPTION"]:
		for key, desc in config["COLUMNS_DESCRIPTION"]["ALT"].items():
			header.append('##ALT=<ID=' + key
							+ ',Description="' + desc + '">')
	if "INFO" in config["COLUMNS_DESCRIPTION"]:
		for key, dic in config["COLUMNS_DESCRIPTION"]["INFO"].items():
			header.append('##INFO=<ID=' + key + ',Number=1,Type='
							+ dic["Type"] + ',Description="' + dic["Description"] + '">')
	if "FORMAT" in config["COLUMNS_DESCRIPTION"]:
		for key, dic in config["COLUMNS_DESCRIPTION"]["FORMAT"].items():
			header.append('##FORMAT=<ID=' + key + ',Number=1,Type='
							+ dic["Type"] + ',Description="' + dic["Description"] + '">')

	header += config["GENOME"]["vcf_header"]
	header.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT",
								"QUAL", "FILTER", "INFO", "FORMAT"] + sample_list))
	return header

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='python fileconverter.py')
	parser.add_argument("-i", "--inputFile", type=str, required=True, help="Input file")
	parser.add_argument("-o", "--outputFile", type=str, required=True, help="Output file")
	parser.add_argument("-fi", "--inputFormat", type=str, required=True, help="Input file format")
	parser.add_argument("-fo", "--outputFormat", type=str, required=True, help="Output file format")
	parser.add_argument("-c", "--configFile", type=str, required=True,
							help="JSON config file describing columns. See script's docstring."
						)
	args = parser.parse_args()

	if args.inputFormat.lower() == "decon":
		print("[ERROR] DECON is handled as a TSV conversion. Use 'tsv' as input format")
		sys.exit()
	c = Converter(args.inputFormat.lower(), args.outputFormat.lower(), args.configFile)
	c.convert(args.inputFile, args.outputFile)
