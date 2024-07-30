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

def correct_svlen(info, fields):
	if 'SVTYPE=DEL' in info:
		match = re.search(r'END=(\d+)', info)
		if match:
			end_pos = int(match.group(1))
			start_pos = int(fields[1])
			svlen = start_pos - end_pos  # Calculate the length, should be negative for deletions

			# Ensure SVLEN is negative
			if svlen > 0:
				svlen = -svlen

			# Extract the existing SVLEN value if it exists
			svlen_match = re.search(r'SVLEN=(-?\d+)', info)
			if svlen_match:
				existing_svlen = int(svlen_match.group(1))
				if existing_svlen != svlen:
					# Update the SVLEN value
					info = re.sub(r'SVLEN=[-]?\d+', f'SVLEN={svlen}', info)
			else:
				# If SVLEN does not exist, add it
				info += f';SVLEN={svlen}'
	return info

def process_vcf(input_file, output_file, is_gz, compress_output):
	open_func = gzip.open if is_gz else open
	with open_func(input_file, 'rt') as infile, tempfile.NamedTemporaryFile(delete=False, mode='w') as tmpfile:
		for line in infile:
			if line.startswith('#'):
				tmpfile.write(line)
			else:
				columns = line.strip().split('\t')
				info_field = re.sub(r'\s+', '', columns[7])  # Remove all types of whitespace characters
				info_field = info_field.replace('\xa0', '')  # Remove non-breaking spaces
				info_field = clean_info_field(info_field)
				info_field = correct_svlen(info_field, columns)
				columns[7] = info_field
				tmpfile.write('\t'.join(columns) + '\n')
		tmpfile_path = tmpfile.name

	if compress_output:
		temp_output_file = output_file + ".tmp.vcf"
		shutil.move(tmpfile_path, temp_output_file)
		# Compress the output file using bgzip
		subprocess.run(['bgzip', '-f', temp_output_file])
		shutil.move(temp_output_file + ".gz", output_file)
	else:
		shutil.move(tmpfile_path, output_file)

def main():
	parser = argparse.ArgumentParser(description='Clean INFO fields and correct SVLEN for deletions in a VCF file and output a bgzipped VCF file.')
	parser.add_argument('-i', '--input', required=True, help='Input VCF or VCF.GZ file')
	parser.add_argument('-o', '--output', required=True, help='Output VCF or VCF.GZ file')
	parser.add_argument('-z', '--gzip', action='store_true', help='Compress the output file using bgzip')
	args = parser.parse_args()

	input_file = args.input
	output_file = args.output
	compress_output = args.gzip or output_file.endswith('.gz')

	is_gz = input_file.endswith('.gz')

	process_vcf(input_file, output_file, is_gz, compress_output)

if __name__ == "__main__":
	main()