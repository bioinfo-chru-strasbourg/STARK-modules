##########################################################################
# split_tsv.py		Version: 1.3
# Description:		Python script to split a TSV based on a specific column name
##########################################################################

########## Note ########################################################################################
# DEV v1 19/10/2024
# Need testing EXPERIMENTAL!
########################################################################################################

# usage
# python split_tsv.py -i path_to_your_input_file.tsv --sample_column Column_Name --output_template path_to_output_file_template.tsv

import argparse
import os
import pandas as pd

def split_tsv(input_file, sample_column, output_template):
	
	"""
	Split a TSV file based on unique values in the specified sample column and save each group to a separate TSV file.
	
	Args:
		input_file (str): Path to the input TSV file.
		sample_column (str): Column name in the TSV file to group the data by.
		output_template (str): Template for the output file path where {sample_id} is replaced by the actual sample ID.
	
	The function splits the data into separate files by the unique values in the specified column.
	The file names are constructed based on the provided template, where "{sample_id}" will be replaced by
	the actual sample ID.
	"""
	
	# Load the TSV file
	df = pd.read_csv(input_file, sep='\t')

	# Group the DataFrame by the specified sample column
	grouped = df.groupby(sample_column)

	# Iterate over each group and save it to a separate TSV file
	for sample_id, group in grouped:
		# Replace {sample_id} in the output template with the actual sample ID
		output_file = output_template.format(sample_id=sample_id)
		
		# Ensure the directory for the output file exists
		os.makedirs(os.path.dirname(output_file), exist_ok=True)
		
		# Save the group to a TSV file
		group.to_csv(output_file, sep='\t', index=False)
		print(f"Saved: {output_file}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Split a TSV file by specified sample column.")
	parser.add_argument('-i', '--input', required=True, help='Input TSV file path.')
	parser.add_argument('--sample_column', default='Sample', help='Column name to group by (default: Sample).')
	parser.add_argument('--output_template', required=True, help='Output file path template (use {sample_id} for the sample identifier).')

	args = parser.parse_args()

	split_tsv(args.input, args.sample_column, args.output_template)