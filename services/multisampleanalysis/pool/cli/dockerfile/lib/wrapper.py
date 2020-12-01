#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Goal: Start the pools pipeline and include its result properly in a STARK repository

@Author: Samuel Nicaise (2020)

#Build the Image
bioinfo@int0663/home1/TOOLS/tools/pool/1.2$ docker build . -t pool:1.2

/home1/TOOLS/tools/pool/dev/bin/pool sample \
	-o /home1/data/WORK_DIR_VINCENT/tmp/test01 \
	-s "/home1/L_PROD/DIAG/DIAG/DI/180112_NB551027_0230_AHWFMMAFXX/ASG160332/ASG160332.final.vcf,/home1/L_PROD/DIAG/DIAG/DI/180112_NB551027_0230_AHWFMMAFXX/ADN160249/ADN160249.final.vcf" \
	-p "POOL_F:/home1/L_PROD/DIAG/DIAG/DI/180112_NB551027_0230_AHWFMMAFXX/POOL_AHWFMMAFXX_F_13/POOL_AHWFMMAFXX_F_13.final.vcf:/home1/L_PROD/DIAG/DIAG/DI/180112_NB551027_0230_AHWFMMAFXX/POOL_AHWFMMAFXX_F_13/DATA/POOL_AHWFMMAFXX_F_13.bwamem.bam,POOL_M:/home1/L_PROD/DIAG/DIAG/DI/180112_NB551027_0230_AHWFMMAFXX/POOL_AHWFMMAFXX_M_12/POOL_AHWFMMAFXX_M_12.final.vcf:/home1/L_PROD/DIAG/DIAG/DI/180112_NB551027_0230_AHWFMMAFXX/POOL_AHWFMMAFXX_M_12/DATA/POOL_AHWFMMAFXX_M_12.bwamem.bam" \
	-b /home1/L_PROD/DIAG/DIAG/DI/180112_NB551027_0230_AHWFMMAFXX/ASG160332/DATA/ASG160332.bed \
	-g /home1/TOOLS/genomes/hg19/hg19.fa

##The following should be done with the results:
#copier le configfile VaRank:
cp /home1/TOOLS/tools/varank/VaRank_1.4.3/configfile.pool configfile
#lancer VaRank:
/home1/TOOLS/tools/stark/current/bin/launch.VaRank.sh -a /home1/data/WORK_DIR_VINCENT/dev/TEST_tristan
#suivre le log
tail -f /home1/data/WORK_DIR_VINCENT/dev/TEST_tristan/VaRank.log
#copier les fichiers tsv sur EMCLABO
cp /home1/data/WORK_DIR_VINCENT/dev/TEST_tristan/*.tsv /home1/L_PROD/DIAG/DIAG/GENODENT/
"""

from __future__ import division
from __future__ import print_function

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import time
from os.path import join as osj

def find_any_samplesheet(runDir, fromResDir = False):
	"""
	Input:
		runDir: String ; path/to/run where you want to find a samplesheet
		fromResDir: boolean ; if True the run dir is organized as a STARK result dir,
							else as a STARK repository dir (i.e. patient dirs contains a STARK dir)
	Output:
		String with samplesheet path or "NO_SAMPLESHEET_FOUND"

	Method:
	1) look up recursively all files named SampleSheet.csv in the runDir
	2) check if file path follows an expected samplesheet name and location
		(the latter depends on if we're in a STARK result or repository dir,
		defined by the bool fromResDir)
	3) first correct file path is returned
	"""
	p = subprocess.Popen("find -L "+runDir+" -maxdepth 3 -name *SampleSheet.csv", stdout=subprocess.PIPE, shell=True)
	out = p.stdout.readlines()
	for ss in out:
		ss = ss.decode("utf-8").strip()
		if fromResDir:
			r = re.match(runDir.rstrip("/")+"/(.*)/(.*).SampleSheet.csv", ss)
		else:
			r = re.match(runDir.rstrip("/")+"/(.*)/STARK/(.*).SampleSheet.csv", ss)
		if r is None:
			continue
		elif r.group(1) == r.group(2): #checks if (.*) == (.*)
			return ss
	return "NO_SAMPLESHEET_FOUND"

def get_sample_list_from_samplesheet(samplesheetPath):
	"""
	Returns a python list containing all sample names in a samplesheet.
	"""
	assert samplesheetPath != "NO_SAMPLESHEET_FOUND", \
			"[ERROR] find_any_samplesheet() couldn't find any samplesheet. Check if the --fromResultDir argument is set correctly."
	inDataTable = False
	sampleList = []
	with open(samplesheetPath, "r") as f:
			for l in f:
				if not inDataTable:
					if l.startswith("Sample_ID,Sample_Name,"):
						inDataTable = True
				else:
					if "," in l:
						sampleList.append(l.strip().split(",")[0])
	#if there are spaces in samplesheet names, change them to "_" because that's what demultiplexing.sh will do
	#otherwise the fastq won't be found when looking in the DEM dir
	for i in range(len(sampleList)):
		if " " in sampleList[i]:
			sampleList[i] = sampleList[i].replace(" ", "_")
	return sampleList

def main(args):
	dockerOutputDir = "/app/res"
	# dockerOutputDir = osj(args.runDir,"res")
	if not os.path.exists(dockerOutputDir):
		os.mkdir(dockerOutputDir)
	vcfList = []
	poolFStr = "init"
	poolMStr = "init"
	bed = "init"
	
	sampleList = get_sample_list_from_samplesheet(find_any_samplesheet(args.runDir))
	removedSamples = []

	for s in sampleList:
		if osj(args.runDir,s) not in glob.glob(osj(args.runDir, "*")):
			print("[WARNING] Sample "+s+" from samplesheet not in repository, ignoring it. Could be normal if it belongs to a different application")
			removedSamples.append(s)
		else:
			vcf = osj(args.runDir, s, "STARK", s+".reports", s+".final.vcf")
			bam = osj(args.runDir, s, "STARK", s+".bwamem.bam")
			# tagFile = osj(args.runDir, s, s+".tag")
			vcfList.append(vcf)
			if re.match("(POOL_[A-Z0-9]*_M_[0-9]*)", s):
				assert poolMStr == "init", "[ERROR] More than one sample is named POOL_([A-Z]*)_M_([0-9]*) in the samplesheet"
				poolMStr = "POOL_M:"+vcf+":"+bam
			if re.match("(POOL_[A-Z0-9]*_F_[0-9]*)", s):
				assert poolFStr == "init", "[ERROR] More than one sample is named POOL_([A-Z]*)_F_([0-9]*) in the samplesheet"
				poolFStr = "POOL_F:"+vcf+":"+bam
				bed = osj(args.runDir, s, "STARK", s+".bed")
	
	for removed in removedSamples:
		sampleList.remove(removed)
	
	assert poolFStr != "init"
	assert poolMStr != "init"
	assert bed != "init"
	#TODO: be able to fetch a pool from a different run ?
	cmd = '/app/bin/pool sample -o '+dockerOutputDir+' -s "'+','.join(vcfList)+'" -p "'+poolFStr+","+poolMStr+'" -b '+bed+' -g '+args.genome
	# cmd = '/home1/TOOLS/tools/pool/1.2/bin/pool sample -o '+dockerOutputDir+' -s "'+','.join(vcfList)+'" -p "'+poolFStr+","+poolMStr+'" -b '+bed+' -g '+args.genome
	print(cmd)
	subprocess.call(cmd, shell=True)
	for s in sampleList:
		resDir = osj(args.runDir, s, "POOL")
		logDir = osj(resDir, "logs")
		if not os.path.exists(resDir):
			os.mkdir(resDir)
		if not os.path.exists(logDir):
			os.mkdir(logDir)
		shutil.copyfile(osj(dockerOutputDir, s+".final.vcf.gz"), osj(resDir, s+".final.vcf.gz"))
		shutil.copyfile(osj(dockerOutputDir, s+".final.vcf.gz.tbi"), osj(resDir, s+".final.vcf.gz.tbi"))
		shutil.copyfile(osj(dockerOutputDir, ".log", s+".bcftools.log"), osj(logDir, s+".bcftools.log"))
		shutil.copyfile(osj(dockerOutputDir, ".log", s+".howard.log"), osj(logDir, s+".howard.log"))
		for f in glob.glob(osj(dockerOutputDir, ".snakemake", "log", "*.log")):
			shutil.copyfile(f, osj(logDir, os.path.basename(f)))
	shutil.rmtree(dockerOutputDir)
	with open(osj(args.runDir, "POOLComplete.txt"), "w") as f:
		f.write(time.ctime())
	if  os.path.exists(osj(args.runDir, "POOLRunning.txt")):
		os.remove(osj(args.runDir, "POOLRunning.txt"))

if __name__=="__main__":
	parser = argparse.ArgumentParser(prog='wrapper to launch routine analysis')
	parser.set_defaults(mode=main)
	parser.add_argument("-i", "--runDir", type=str, help="path to run in a STARK 0.9.18 repository",required=True)
	parser.add_argument("-g","--genome", help="genome file", type=str, dest='genome',required=True)
	
	args = parser.parse_args()
	args.mode(args)