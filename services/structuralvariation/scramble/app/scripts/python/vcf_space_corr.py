import argparse
import gzip
import re
import os
import subprocess

def remove_whitespace(vcf_file, output_file):
	# Detect if the input file is gzipped
	is_gzipped = vcf_file.endswith('.gz')
	
	# Open the input file based on its type
	open_func = gzip.open if is_gzipped else open
	
	with open_func(vcf_file, 'rt') as infile, open(output_file, 'wt') as outfile:
		for line in infile:
			if line.startswith('#'):
				outfile.write(line)
			else:
				fields = line.strip().split('\t')
				info_field = fields[7]

				# Remove all types of whitespace characters, including non-breaking spaces
				info_field = re.sub(r'\s+', '', info_field)
				info_field = info_field.replace('\xa0', '')  # Remove non-breaking spaces
				
				fields[7] = info_field
				outfile.write('\t'.join(fields) + '\n')
	
	# Compress the output file using bgzip
	subprocess.run(['bgzip', '-f', output_file])

def main():
	parser = argparse.ArgumentParser(description='Remove whitespace in the INFO field of a VCF file (plain or gzipped).')
	parser.add_argument('input_vcf', help='Input VCF file (.vcf or .vcf.gz)')
	parser.add_argument('output_vcf', help='Output VCF file with removed whitespace (will be bgzipped)')
	
	args = parser.parse_args()
	
	output_vcf_temp = args.output_vcf.rstrip('.gz')
	remove_whitespace(args.input_vcf, output_vcf_temp)
	
	# Rename the temp output file to the desired output file name
	os.rename(output_vcf_temp + '.gz', args.output_vcf)

if __name__ == '__main__':
	main()