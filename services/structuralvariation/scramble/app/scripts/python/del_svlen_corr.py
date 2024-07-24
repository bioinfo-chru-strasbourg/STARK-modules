import argparse
import gzip
import re
import os
import subprocess

def correct_svlen(vcf_file, output_file):
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

				# Ensure there are no extra spaces in the INFO field
				info_field = info_field.replace(" ", "")

				if 'SVTYPE=DEL' in info_field:
					match = re.search(r'END=(\d+)', info_field)
					if match:
						end_pos = int(match.group(1))
						start_pos = int(fields[1])
						svlen = start_pos - end_pos  # Calculate the negative length
						# If SVLEN exists, replace it; otherwise, add it
						if 'SVLEN=' in info_field:
							info_field = re.sub(r'SVLEN=[^;]+', f'SVLEN={svlen}', info_field)
						else:
							info_field += f';SVLEN={svlen}'
						fields[7] = info_field
				
				# Ensure the INFO field is properly formatted
				fields[7] = re.sub(r'\s+', '', fields[7])
				
				outfile.write('\t'.join(fields) + '\n')
	
	# Compress the output file using bgzip
	subprocess.run(['bgzip', '-f', output_file])

def main():
	parser = argparse.ArgumentParser(description='Correct SVLEN for deletions in a VCF file (plain or gzipped).')
	parser.add_argument('input_vcf', help='Input VCF file (.vcf or .vcf.gz)')
	parser.add_argument('output_vcf', help='Output VCF file with corrected SVLEN (will be bgzipped)')
	
	args = parser.parse_args()
	
	output_vcf_temp = args.output_vcf.rstrip('.gz')
	correct_svlen(args.input_vcf, output_vcf_temp)
	
	# Rename the temp output file to the desired output file name
	os.rename(output_vcf_temp + '.gz', args.output_vcf)

if __name__ == '__main__':
	main()