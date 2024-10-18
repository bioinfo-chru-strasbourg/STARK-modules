##########################################################################
# split_tsv.py		Version: 1.0
# Description:		Python script split a tsv based on a specific column name
##########################################################################

########## Note ########################################################################################
# DEV v1 19/10/2024
# Need testing EXPERIMENTAL!
########################################################################################################

# usage
# python split_tsv.py -i path_to_your_input_file.tsv -o path_to_output_directory -s path_to_sample_list.txt


import argparse
import pandas as pd

def split_tsv(input_file, output_dir, sample_list_file, sample_column):
	
	"""
	Split a TSV file based on unique values in the specified sample column and save each group to a separate TSV file.
	
	Args:
		input_file (str): Path to the input TSV file.
		output_dir (str): Path to the output directory where split files will be saved.
		sample_list_file (str): Path to a file containing the list of valid sample IDs to filter.
		sample_column (str): Column name in the TSV file to group the data by (default is 'Samples_ID').
	
	The function filters the input TSV file based on the provided sample list, splits the data
	into separate files by the unique values in the specified column, and saves them to the 
	output directory. Each output file is named after the unique sample ID and retains the 
	header from the input TSV file.
	"""
	
	# Load the TSV file
	df = pd.read_csv(input_file, sep='\t')

	# Load the sample list
	with open(sample_list_file, 'r') as f:
		sample_list = [line.strip() for line in f]

	# Filter the DataFrame based on the sample list
	df_filtered = df[df[sample_column].isin(sample_list)]

	# Group the DataFrame by the specified sample column
	grouped = df_filtered.groupby(sample_column)

	# Iterate over each group and save it to a separate TSV file
	for sample_id, group in grouped:
		# Construct the output file path
		filename = f"{output_dir}/{sample_id}.tsv"
		# Save the group to a TSV file
		group.to_csv(filename, sep='\t', index=False)
		print(f"Saved: {filename}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Split a TSV file by specified sample column.")
	parser.add_argument('-i', '--input', required=True, help='Input TSV file path.')
	parser.add_argument('-o', '--output', required=True, help='Output directory path.')
	parser.add_argument('-s', '--samples', required=True, help='File with the list of valid Samples_IDs.')
	parser.add_argument('--sample_column', default='Samples_ID', help='Column name to group by (default: Samples_ID).')

	args = parser.parse_args()

	split_tsv(args.input, args.output, args.samples, args.sample_column)