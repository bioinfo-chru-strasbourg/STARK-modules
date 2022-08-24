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
import doctest
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
					if l.startswith("Sample_ID,"):
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

# def getDescriptionFromSample(sample, runDir):
	# samplesheetPath = find_any_samplesheet(runDir)
	# assert samplesheetPath != "NO_SAMPLESHEET_FOUND", \
			# "[ERROR] find_any_samplesheet() couldn't find any samplesheet. Check if the --fromResultDir argument is set correctly."
	# inDataTable = False
	# description = []
	# with open(samplesheetPath, "r") as f:
			# for l in f:
				# if not inDataTable:
					# if l.startswith("Sample_ID,Sample_Name,"):
						# inDataTable = True
						# DescriptionIndex = l.strip().split(",").index("Description")
				# else:
					# if sample in l:
						# description.append(l.strip().split(",")[DescriptionIndex])
	# return description

# def getPool(description, sexTag):
	# if re.search("APP#[A-Z0-9]*.[A-Z0-9]*#POOL", description) and sexTag in description:
		# return True
	# else:
		# return False

def cleanList(str):
	return [ x.strip() for x in str.split(",") ]

def isValidSample(runDir, sample, exclude):
	with open(osj(runDir, sample,"STARK", sample+".tag"), "r") as f:
		for l in f: 
			for item in exclude:
				if item in l:
					return False
	return True

def tag_field_to_dict(tag_string):
	"""
	>>> tag_field_to_dict("APP#HUS")
	{'APP': ['HUS']}
	>>> tag_field_to_dict("APP#HUSTUMSOL.XTHS!#test!ceci est une description")
	{'': ['test'], 'APP': ['HUSTUMSOL.XTHS']}
	>>> tag_field_to_dict("SEX#F!APP#DIAG.DI#POOL!")
	{'APP': ['DIAG.DI', 'POOL'], 'SEX': ['F']}
	
	Trailing "!" are authorized and must be dealt with.
	Fields without a # are not returned. 
	"""
	tag_dict = {}
	tag_string = tag_string.rstrip("\r\n") # as it will usually come from a samplesheet
	if "#" not in tag_string:
		return {}
	if "!" not in tag_string:
		v = tag_string.split("#")
		if len(v) == 2:
			tag_dict[v[0]] = [v[1]]
		elif len(v) > 2:
			tag_dict[v[0]] = v[1:]
	else:
		for ts in tag_string.split("!"):
			if ts == "": #deals with trailing "!"
				continue
			v = ts.split("#")
			if len(v) == 2:
				tag_dict[v[0]] = [v[1]]
			elif len(v) > 2:
				tag_dict[v[0]] = v[1:]
	return tag_dict

def is_pool(runDir, sample, sexTag):
	"""
	checks if sample is a pool of samples of given sex according to sample tag file
	"""
	isPool = False
	correctSex = False
	with open(osj(runDir, sample, "STARK", sample+".tag"), "r") as f:
		for l in f:
			tagsDict = tag_field_to_dict(l)
			for key in tagsDict.keys():
				if key == 'PLUGAPP':
					if 'POOL' in tagsDict[key]:
						isPool = True
				elif key == 'APP':
					if 'POOL' in tagsDict[key]:
						isPool = True
				if key == 'SEX':
					if tagsDict[key] == tag_field_to_dict(sexTag)[key]:
						correctSex = True
	if isPool and correctSex:
		return True
	else:
		return False

def getPoolDict(sampleList, runDir):
	poolDict = {}
	for s in sampleList:
		with open(osj(runDir, s, "STARK", s+".tag"), "r") as f:
			for l in f:
				if 'POOL' in tag_field_to_dict(l).keys():
					poolIndex = tag_field_to_dict(l)['POOL']
					poolIndex.sort()
					if not "#".join(poolIndex) in poolDict.keys():
						poolDict["#".join(poolIndex)] = []
					poolDict["#".join(poolIndex)].append(s)
	if not poolDict:
		poolList = []
		for s in sampleList:
			with open(osj(runDir, s, "STARK", s+".tag"), "r") as f:
				for l in f:
					if is_pool(runDir, s, "SEX#F"):
						poolList.append(s)
					if is_pool(runDir, s, "SEX#M"):
						poolList.append(s)
		poolList.sort()
		poolDict["#".join(poolList)] = sampleList
	return poolDict

def createSampleRepository(runDir, sample):
	resDir = osj(runDir, sample, "POOL")
	logDir = osj(resDir, "logs")
	if not os.path.exists(resDir):
		os.mkdir(resDir)
	if not os.path.exists(logDir):
		os.mkdir(logDir)
	return (resDir, logDir)

def writeErrorLog(sample, runDir, errorMessage):
	resDir, logDir = createSampleRepository(runDir, sample)
	with open(osj(resDir, "ERROR.log"), "w") as f:
		f.write(errorMessage)

def bcftoolsCompress(vcf, folder):
	vcfCompressed = osj(folder,os.path.splitext(os.path.basename(vcf))[0]+".vcf.gz")
	cmd = "/STARK/tools/htslib/current/bin/bgzip -c "+vcf+" > "+vcfCompressed
	subprocess.call(cmd, shell = True)
	cmd = "/STARK/tools/htslib/current/bin/tabix -f "+vcfCompressed
	subprocess.call(cmd, shell = True)
	return vcfCompressed

def createEmptyVcfAndBam(pool, dockerOutputDir, sampleList, runDir, bcftools, samtools):
	if not os.path.exists(osj(dockerOutputDir, pool)):
		os.mkdir(osj(dockerOutputDir, pool))
	vcf = osj(dockerOutputDir, pool, pool+".vcf")
	bam = osj(dockerOutputDir, pool, pool+".bam")
	vcfCompressed = bcftoolsCompress(osj(runDir, sampleList[0], "STARK", sampleList[0]+".reports", sampleList[0]+".final.vcf"), dockerOutputDir)
	cmd = bcftools+' view -h -o '+vcf+' -O v '+osj(runDir,sampleList[0],sampleList[0]+'.final.vcf.gz')
	subprocess.call(cmd, shell=True)
	cmd = samtools+' view -H -b -o '+bam+' '+osj(runDir,sampleList[0],"STARK",sampleList[0]+'.bwamem.bam')
	subprocess.call(cmd, shell=True)
	cmd = samtools+' index -b '+bam+' '+osj(bam+'.bai')
	subprocess.call(cmd, shell=True)
	poolStr = pool+":"+vcf+":"+bam
	return poolStr

def copyResults(sampleList, runDir, dockerOutputDir):
	for s in sampleList:
		if os.path.exists(osj(dockerOutputDir, s+".final.vcf.gz")):
			resDir, logDir = createSampleRepository(runDir, s)
			shutil.copyfile(osj(dockerOutputDir, s+".final.vcf.gz"), osj(resDir, s+".final.vcf.gz"))
			shutil.copyfile(osj(dockerOutputDir, s+".final.vcf.gz.tbi"), osj(resDir, s+".final.vcf.gz.tbi"))
			shutil.copyfile(osj(dockerOutputDir, ".log", s+".bcftools.log"), osj(logDir, s+".bcftools.log"))
			shutil.copyfile(osj(dockerOutputDir, ".log", s+".howard.log"), osj(logDir, s+".howard.log"))
			for f in glob.glob(osj(dockerOutputDir, ".snakemake", "log", "*.log")):
				shutil.copyfile(f, osj(logDir, os.path.basename(f)))
		else:
			writeErrorLog(s, runDir,  "[ERROR] Missing final VCF file. Are the input BAM & VCF correct?")

def launchAnalysis(sampleList, key, runDir, dockerOutputDir, genome, bcftools, samtools):
	vcfList = []
	poolFStr = "init"
	poolMStr = "init"
	bed = "init"
	for s in sampleList:
		vcf = osj(runDir, s, "STARK", s+".reports", s+".final.vcf")
		# bam = osj(runDir, s, "STARK", s+".bwamem.bam")
		vcfList.append(vcf)
	for pool in key.split('#'):
		if os.path.exists(osj(runDir, pool, "STARK",pool+".tag")):
			if is_pool(runDir, pool, "SEX#M"):
				assert poolMStr == "init", "[ERROR] More than one sample is named POOL_([A-Z]*)_M_([0-9]*) in the samplesheet"
				vcf = osj(runDir, pool, "STARK", pool+".reports", pool+".final.vcf")
				bam = osj(runDir, pool, "STARK", pool+".bwamem.bam")
				vcfList.append(vcf)
				poolMStr = "POOL_M:"+vcf+":"+bam
			elif is_pool(runDir, pool, "SEX#F"):
				assert poolFStr == "init", "[ERROR] More than one sample is named POOL_([A-Z]*)_F_([0-9]*) in the samplesheet"
				vcf = osj(runDir, pool, "STARK", pool+".reports", pool+".final.vcf")
				bam = osj(runDir, pool, "STARK", pool+".bwamem.bam")
				vcfList.append(vcf)
				poolFStr = "POOL_F:"+vcf+":"+bam
			else:
				for s in sampleList:
					writeErrorLog(s, runDir, "[ERROR] Can't find "+pool+" sample's sex.")
				return "[ERROR] Can't find "+pool+" sample's sex."
		else:
			for s in sampleList:
				writeErrorLog(s, runDir, "[ERROR] Specified "+pool+" sample not found.")
			return "[ERROR] Specified "+pool+" sample not found."
	bed = osj(runDir, pool, "STARK", pool+".bed")
	# assert poolFStr != "init"
	# assert poolMStr != "init"
	# assert bed != "init"
	if poolFStr == "init":
		poolFStr = createEmptyVcfAndBam("POOL_F", dockerOutputDir, sampleList, runDir, bcftools, samtools)
		# for s in sampleList:
			# writeErrorLog(s, runDir, "[ERROR] Missing POOL_F for POOL analysis.")
		# return "[ERROR] Missing POOL_F for POOL analysis."
	if poolMStr == "init":
		poolMStr = createEmptyVcfAndBam("POOL_M", dockerOutputDir, sampleList, runDir, bcftools, samtools)
		# for s in sampleList:
			# writeErrorLog(s, runDir, "[ERROR] Missing POOL_M for POOL analysis.")
		# return "[ERROR] Missing POOL_M for POOL analysis."
	if bed == "init":
		for s in sampleList:
			writeErrorLog(s, runDir, "[ERROR] Missing bed for POOL analysis.")
		return "[ERROR] Missing bed for POOL analysis."
	#TODO: be able to fetch a pool from a different run ?
	cmd = 'python /app/lib/pool/pool.py sample -o '+dockerOutputDir+' -s "'+','.join(vcfList)+'" -p "'+poolFStr+","+poolMStr+'" -b '+bed+' -g '+genome
	# cmd = 'python /home1/TOOLS/tools/pool/1.2/lib/pool/pool.py sample -o '+dockerOutputDir+' -s "'+','.join(vcfList)+'" -p "'+poolFStr+","+poolMStr+'" -b '+bed+' -g '+genome
	print(cmd)
	subprocess.call(cmd, shell=True)
	copyResults(sampleList, runDir, dockerOutputDir)

def main(args):
	dockerOutputDir = "/app/res"
	# dockerOutputDir = osj(args.runDir,"res")
	if not os.path.exists(dockerOutputDir):
		os.mkdir(dockerOutputDir)
	
	#get only samples that will be analysed
	sampleList = get_sample_list_from_samplesheet(find_any_samplesheet(args.runDir))
	removedSamples = []
	for s in sampleList:
		if not isValidSample(args.runDir, s, cleanList(args.exclude)):
			removedSamples.append(s)
		elif osj(args.runDir,s) not in glob.glob(osj(args.runDir, "*")):
			print("[WARNING] Sample "+s+" from samplesheet not in repository, ignoring it. Could be normal if it belongs to a different application")
			removedSamples.append(s)
	for removed in removedSamples:
		sampleList.remove(removed)
	#create a dictionary with POOL_ID1#POOL_ID2 as key, and samples as values, depending on tag POOL#POOL_ID1#POOL_ID2
	poolDict = getPoolDict(sampleList, args.runDir)
	for key in poolDict:
		launchAnalysis(poolDict[key], key, args.runDir, dockerOutputDir, args.genome, args.bcftools, args.samtools)
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
	parser.add_argument("-e","--exclude", help="list of tags for which samples are removed if have them", type=str, dest='exclude', default="CQI#")
	parser.add_argument("-b","--bcftools", help="path to bcftools bin", type=str, dest='bcftools', default="/STARK/tools/bcftools/current/bin/bcftools")
	parser.add_argument("-s","--samtools", help="path to samtools bin", type=str, dest='samtools', default="/STARK/tools/samtools/current/bin/samtools")
	
	args = parser.parse_args()
	args.mode(args)