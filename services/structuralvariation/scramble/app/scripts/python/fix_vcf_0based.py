##########################################################################
# fix_vcf_0based.py       	Version: 1.0
# Description:          	Python script to adjust positions in a VCF file based on variant type and handle .gz compression
##########################################################################

########## Note ########################################################################################
# DEV v1 27/06/2024
# Changelog
#   - v1.0: Initial version
########################################################################################################

import argparse
import pysam
import gzip

def increment_positions(input_file, output_file, gzipped_output):
	# Open the VCF file using pysam
	if input_file.endswith('.gz'):
		vcf_in = pysam.VariantFile(gzip.open(input_file, 'rt'), 'r')
	else:
		vcf_in = pysam.VariantFile(input_file, 'r')
	
	# Open the output file for writing
	if gzipped_output:
		vcf_out = pysam.VariantFile(gzip.open(output_file, 'wt'), 'w', header=vcf_in.header)
	else:
		vcf_out = pysam.VariantFile(output_file, 'w', header=vcf_in.header)
	
	# Process each record
	for rec in vcf_in.fetch():
		# Get the variant type from INFO field
		var_type = rec.info.get('TYPE', 'SNP')  # Default to 'SNP' if TYPE is not present
		
		if var_type == 'DEL':
			# Increment start position only for deletions
			rec.pos += 1
			if 'END' in rec.info:
				rec.info['END'] = rec.info['END']  # Keep END as is
		elif var_type == 'INS':
			# Increment end position only for insertions
			if 'END' in rec.info:
				rec.info['END'] = rec.info['END'] + 1  # Increment end position for insertions
		else:
			# For other types, no changes
			pass
		
		# Write the modified record to the output VCF file
		vcf_out.write(rec)
	
	# Close the input and output files
	vcf_in.close()
	vcf_out.close()

def main():
	parser = argparse.ArgumentParser(description='Adjust positions in a VCF file based on variant type and handle .gz compression.')
	parser.add_argument('-i', '--input', required=True, help='Path to the input VCF file (can be gzipped)')
	parser.add_argument('-o', '--output', required=True, help='Path to the output VCF file')
	parser.add_argument('-z', '--gzip', action='store_true', help='Compress the output VCF file with gzip')

	args = parser.parse_args()

	increment_positions(args.input, args.output, args.gzip)

if __name__ == '__main__':
	main()