##########################################################################
# merge_tsv.py		       	Version: 1.0
# Description:          	Python script to merge tsv files, combining headers and optionally filling blank columns with 'NA'.
##########################################################################

########## Note ########################################################################################
# DEV v1 27/06/2024
# Changelog
#   - v1.0: Initial version
########################################################################################################

# usage python merge_tsv.py -i file1.tsv file2.tsv file3.tsv -o merged_output.tsv

import pandas as pd
import argparse
import os

def merge_tsv(input_files, output_file, fill_na=False):
	# Read and merge TSV files
	dfs = [pd.read_csv(file, sep='\t') for file in input_files]
	merged_df = pd.concat(dfs, ignore_index=True, sort=False)

	# Optionally fill missing values with 'NA' for non-numeric columns
	if fill_na:
		for col in merged_df.columns:
			if merged_df[col].dtype == 'float64':
				merged_df[col] = merged_df[col].fillna(pd.NA)
			else:
				merged_df[col] = merged_df[col].fillna('NA')
	
	# Write the merged DataFrame to the output file
	merged_df.to_csv(output_file, sep='\t', index=False)

def main():
	# Argument parser
	parser = argparse.ArgumentParser(description="Merge multiple TSV files, combining headers and optionally filling blank columns with 'NA'.")
	parser.add_argument('-i', '--input', nargs='+', required=True, help="Input TSV files to merge.")
	parser.add_argument('-o', '--output', required=True, help="Output TSV file after merging.")
	parser.add_argument('--fill-na', action='store_true', help="Fill missing values with 'NA' in non-numeric columns and NaN in numeric columns.")

	args = parser.parse_args()

	# Check if input files exist
	for file in args.input:
		if not os.path.isfile(file):
			print(f"Error: File '{file}' does not exist.")
			return

	# Merge the TSV files with or without filling NA based on the argument
	merge_tsv(args.input, args.output, fill_na=args.fill_na)
	print(f"Successfully merged files into {args.output}")

if __name__ == "__main__":
	main()