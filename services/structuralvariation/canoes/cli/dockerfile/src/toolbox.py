import pandas as pd
import inspect
import os
import subprocess
from collections import OrderedDict

def systemcall(command):
	'''
	*passing command to the shell*
	*return list containing stdout lines*
	command - shell command (string)
	'''
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read().decode('utf8').strip().split('\n')

def parse_samplesheet(samplesheet_path):
	'''
	taking a samplesheet path, illumina format
	return an ordered dict with data field informations
	'''
	samplesheet_data = OrderedDict()
	samplesheet_header = []
	if os.path.exists(samplesheet_path):
		print("#[INFO] "+os.path.abspath(samplesheet_path+" to ordered dict"))
		#Open samplesheet file
		with open(samplesheet_path, 'r') as f:
			v = False
			#iterate over each row in samplesheet
			for lines in f:
				lines = lines.strip()
				#within data field
				if v:
					#fill each key with samplesheet values
					for i, keys in enumerate(samplesheet_data):
						samplesheet_data[keys].append(lines.split(',')[i])
				#catch data field header	
				if 'Sample_ID' in lines:
					v = True
					#add a list for each header field key
					samplesheet_header.append(lines.split(','))
					for fields in lines.split(','):
						samplesheet_data[fields] = []
	else:
		print("ERROR "+os.path.abspath(samplesheet_path)+" does not exist !!")
		exit()
	return samplesheet_data

def cramTobam(cramfile, bamfile, samtools, genome):
	'''
	Used samtools to convert cram into bamfile
	samtools should be in path or provide abspath of samtools binary
	'''
	if os.path.exists(bamfile):
		print("ERROR "+bamfile+" already exists ")
		exit()
	else:
		print("#[INFO] From "+cramfile+" to "+bamfile)
		systemcall(samtools+" view -C -T "+genome+" -o "+bamfile+" "+cramfile)

cramTobam("/home1/BAS/lamouchj/scripts/B9173.archive.cram", "/home1/BAS/lamouchj/scripts/B9173.bam", "/home1/TOOLS/tools/samtools/1.9/bin/samtools", "/home1/data/STARK/databases/genomes/current/hg19.fa")
#print(parse_samplesheet("/home1/BAS/lamouchj/scripts/TEST/SampleSheet_tumsol_210219.csv"))