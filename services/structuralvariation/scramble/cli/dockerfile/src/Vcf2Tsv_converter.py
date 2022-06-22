##########################################################################
# Vcf2Tsv Version:			1.0
# Description:				Script to convert vcf to a tsv file
##########################################################################

# DEV version 0.1 : 21/06/2022

# Authoring : Jean-baptiste LAMOUCHE

# TL : Add a input and an output option, argparse


################## Context ##################
# type python Vcf2Tsv_converter.py -h for help
# work only for mono sample vcfs
####################################

import pandas as pd
from tqdm import tqdm
import subprocess
import sys
import numpy as np

import argparse
import doctest

def parse_sample_field(dfVar):
	#############
	### Parsing Sample Field in VCF
	#############

	dico = []
	bad_annotation = []
	sample_list = []

	# Parsing FORMAT field in VCF
	# print("[#INFO] Parsing FORMAT field")
	isample = list(dfVar.columns).index("FORMAT") + 1
	# index: line where the caller identify an event somethings
	for col in dfVar.columns[isample:]:
		# print("#[INFO] " + col + "\n")
		sample_list.append(col)
		for i, row in tqdm(
			dfVar.iterrows(),
			total=dfVar.shape[0],
			leave=True,
			desc="#[INFO] SAMPLE column split",
		):
			# print_progress_bar(i, len(dfVar.index)-1)
			# print('\n')
			# print(row['FORMAT'].split(':'), row['bwamem.VarScan_HUSTUMSOL.howard'].split(':'))
			if len(row["FORMAT"].split(":")) != len(row[col].split(":")):
				bad_annotation.append(pd.Series(row[:], index=dfVar.columns))
				continue
			else:
				toadd = pd.Series(
					row[col].split(":"), index=row["FORMAT"].split(":")
				).to_dict()
				# toadd.update({"caller":col})
				dico.append(toadd)

	dfSample = pd.DataFrame(dico)
	df_bad_anno = pd.DataFrame(bad_annotation)
	try:
		df_final = dfVar.join(dfSample, how="inner")
	except ValueError:
		df_final = dfVar.join(dfSample, how="inner", lsuffix="_INFO", rsuffix="_SAMPLE")
		print(
			"WARNING some columns are present in both INFO and SAMPLE field, add _INFO and _SAMPLE suffix"
		)
	df_final.drop(columns=sample_list, inplace=True)
	# Remove FORMAT col
	df_final.drop(columns="FORMAT", inplace=True)
	return df_final


def parse_info_field(dfVar):
	"""
	input: take a dataframe (from vcf)

	output: return a dataframe of the vcf when the info field is parsed
	"""

	############
	# Parsing INFO field from dfVar dataframe containing all informations from vcf
	############

	# print(dfVar.head())
	infoList = []
	dicoInfo = []
	headers = []

	# print("#[INFO] Parsing INFO field")
	for i, elems in tqdm(
		dfVar.iterrows(), total=dfVar.shape[0], desc="#[INFO] INFO column split"
	):
		# print_progress_bar(i, len(dfVar.index)-1)
		infoList.append([x.split("=") for x in elems["INFO"].split(";")])

	# print("\n")

	[headers.append(elems[0]) for ite in infoList for elems in ite]
	dfInfo = pd.DataFrame(columns=np.unique(np.array(headers)))
	# print(np.unique(np.array(headers)))

	# print("#[INFO] From INFO field to Dataframe")
	for j, elems in enumerate(infoList):
		# print_progress_bar(j, len(infoList)-1)
		add = {}
		for fields in elems:
			if len(fields) <= 1:
				f = {fields[0]: "TRUE"}
				add.update(f)
			else:
				f = dict([fields])
				add.update(f)

		dicoInfo.append(add)

	# print("\n")
	# print(dicoInfo.keys())
	# print(dict(list(dicoInfo.items())[0:2]))

	df_final = pd.DataFrame(dicoInfo, columns=np.unique(np.array(headers)))

	dfInfo = dfVar.iloc[:, :7].join(df_final, how="inner")
	# Drop old columns info
	# dfInfo.drop(columns="INFO", inplace=True)
	df = dfInfo.join(dfVar.iloc[:, 8:], how="inner")
	# drop possible no annotation, stack like a '.' col
	if "." in df.columns:
		df.drop(columns=".", inplace=True)
	return df


def systemcall(command):
	"""
	Passing command to the shell, return first item in list containing stdout lines
	"""
	try:
		print("#[INFO] " + command)
		p = subprocess.check_output([command], shell=True)
	except subprocess.CalledProcessError as err:
		sys.exit("ERROR " + str(err.returncode))
	return (
		subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
		.stdout.read()
		.decode("utf8")
		.strip()
		.split("\n")[0]
	)


def preprocess_vcf(file):
	"""
	from vcf file as input
	return dico with 2 lists, header with each row from header and field which contain the vcf header field,
			and number of row of header (without field)
	"""
	data = {"header": []}
	skip = []
	with open(file, "r") as f:
		for i, lines in enumerate(f):
			if lines.split("\t")[0] != "#CHROM":
				data["header"].append(lines.strip())
			else:
				# data["fields"] = lines.strip().split("\t")
				data["fields"] = lines.strip()
				skip.append(i)
				break
	if len(systemcall("grep -v '^#' " + file + " | wc -l")) > 10000:
		df = pd.read_csv(
			file, skiprows=skip[0], sep="\t", header=0, chunksize=1000, low_memory=False
		)
	else:
		df = pd.read_csv(file, skiprows=skip[0], sep="\t", header=0)
	return df


def main(inputfile, outputfile):
	df_variants = parse_sample_field(parse_info_field(preprocess_vcf(inputfile)))
	df_variants.to_csv(outputfile, sep="\t", header=True, index=False)


def myoptions():
	'''
	*arg parser*
	*return options*
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", type = str, default = "", help = "Input vcf file", dest = 'inputfile')
	parser.add_argument("-o", "--output", type = str, default = "", help = "Output tsv file", dest = 'outputfile')
	return parser.parse_args()

if __name__ == "__main__":
	doctest.testmod()
	args = myoptions()
	main(args.inputfile, args.outputfile)