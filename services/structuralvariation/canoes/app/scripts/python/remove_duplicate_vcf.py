import re
import argparse
import gzip
import shutil
import os
import tempfile
import subprocess

def clean_info_field(info):
	fields = info.split(';')
	unique_fields = {}
	
	for field in fields:
		key_value = field.split('=')
		if len(key_value) == 2:
			key, value = key_value
			unique_fields[key] = value
	
	cleaned_info = ';'.join([f"{key}={value}" for key, value in unique_fields.items()])
	return cleaned_info

def process_vcf(input_file, output_file, is_gz):
	open_func = gzip.open if is_gz else open
	with open_func(input_file, 'rt') as infile, tempfile.NamedTemporaryFile(delete=False, mode='w') as tmpfile:
		for line in infile:
			if line.startswith('#'):
				tmpfile.write(line)
			else:
				columns = line.strip().split('\t')
				columns[7] = clean_info_field(columns[7])
				tmpfile.write('\t'.join(columns) + '\n')
		tmpfile_path = tmpfile.name

	temp_output_file = output_file + ".tmp.vcf"
	shutil.move(tmpfile_path, temp_output_file)

	# Compress the output file using bgzip
	subprocess.run(['bgzip', '-f', temp_output_file])
	shutil.move(temp_output_file + ".gz", output_file)

def main():
	parser = argparse.ArgumentParser(description='Clean INFO fields in a VCF file and output a bgzipped VCF file.')
	parser.add_argument('input_vcf', help='Input VCF or VCF.GZ file')
	parser.add_argument('output_vcf', help='Output VCF.GZ file')
	args = parser.parse_args()
	
	input_file = args.input_vcf
	output_file = args.output_vcf
	
	is_gz = input_file.endswith('.gz')
	
	process_vcf(input_file, output_file, is_gz)

if __name__ == "__main__":
	main()