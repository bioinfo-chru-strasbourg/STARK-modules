#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import pandas as pd
import subprocess
import time
import re
import shutil
import time
from datetime import datetime, timedelta
import argparse
import json
from os.path import join as osj



def systemcall(command):
	'''
	*passing command to the shell*
	*return list containing stdout lines*
	command - shell command (string)
	'''
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read().decode('utf8').strip().split('\n')

def nbCPUs():
	command = "lscpu | grep -E '^CPU\\(s\\)' | awk '{print $NF}'"
	return systemcall(command)

def find_samplesheet(runDir, fromResDir = False):
	"""

	"""
	ss = subprocess.Popen("find -L "+runDir+" -maxdepth 3 -name *SampleSheet.csv -print -quit", stdout=subprocess.PIPE, shell=True)
	if ss:
		return ss
	else:
		return "NO_SAMPLESHEET_FOUND"

def parse_samplesheet(samplesheet_path):
	## Sample_ID extraction ##
	samplesheet_data = []
	samplesheet_header = []
	with open(samplesheet_path, 'r') as f:
		v = False
		for lines in f:
			lines = lines.strip()
			if v:
				samplesheet_data.append(lines.split(','))
			if 'Sample_ID' in lines:
				v = True
				samplesheet_header.append(lines.split(','))

	df = pd.DataFrame(samplesheet_data, columns=samplesheet_header)
	
	#sample ID
	sample_list = df.iloc[:, 0].tolist()

	#return dataframe with full SS's samples informations
	return df

#def createCompleteFile(output):
#	file = open(osj(output,"CANOESComplete.txt"), "w+")
#	file.write("# ["+(datetime.now() + timedelta(hours = 2)).strftime("%d/%m/%Y %H:%M:%S")+"] "+os.path.basename(output)+" analyzed with CANOES\n")
#	file.close()
#	if os.path.exists(osj(output,"CANOESRunning.txt")):
#		os.remove(osj(output,"CANOESRunning.txt"))

