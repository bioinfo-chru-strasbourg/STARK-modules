##########################################################################
# split_tsv.py		Version: 1.3
# Description:		Python script to split a TSV based on a specific column name
##########################################################################

########## Note ########################################################################################
# DEV v1 19/10/2024
# Need testing EXPERIMENTAL!
########################################################################################################

# usage
# python split_tsv.py -i path_to_your_input_file.tsv --sample_column Column_Name --output path_to_output_file.tsv --specific_sample Sample_Name

import argparse
import os
import pandas as pd

def split_tsv(input_file, sample_column, output_file, specific_sample):
	"""
	Extract a specific sample from a TSV file and save it to a TSV file.
	
	Args:
		input_file (str): Path to the input TSV file.
		sample_column (str): Column name in the TSV file to filter the data by.
		output_file (str): Output file path where the specified sample will be saved.
		specific_sample (str): The specific sample to extract from the input file.
	"""
	
	# Load the TSV file
	try:
		df = pd.read_csv(input_file, sep='\t')
	except FileNotFoundError:
		print(f"Error: The input file {input_file} does not exist.")
		return
	
	if sample_column not in df.columns:
		print(f"Error: The column '{sample_column}' does not exist in the input file.")
		return
	
	# Filter the DataFrame for the specific sample
	filtered_group = df[df[sample_column] == specific_sample]

	# Check if the filtered group is not empty
	if not filtered_group.empty:
		# Ensure the directory for the output file exists
		os.makedirs(os.path.dirname(output_file), exist_ok=True)

		# Check if the file already exists
		if os.path.exists(output_file):
			print(f"Warning: The file {output_file} already exists and will be overwritten.")

		# Save the filtered group to a TSV file
		filtered_group.to_csv(output_file, sep='\t', index=False)
		print(f"Saved: {output_file}")
	else:
		print(f"No data found for sample: {specific_sample}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Extract a specific sample from a TSV file.")
	parser.add_argument('-i', '--input', required=True, help='Input TSV file path.')
	parser.add_argument('--sample_column', default='Sample', help='Column name to filter by (default: Sample).')
	parser.add_argument('--output', required=True, help='Output file path for the specified sample.')
	parser.add_argument('--specific_sample', required=True, help='The specific sample to extract.')

	args = parser.parse_args()

	split_tsv(args.input, args.sample_column, args.output, args.specific_sample)
