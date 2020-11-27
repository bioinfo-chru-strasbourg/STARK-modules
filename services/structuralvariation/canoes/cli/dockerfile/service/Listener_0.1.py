#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Goal:
- Listener for CANOES module for STARK18

@Author: Victor Grentzinger (2020)
"""

from __future__ import division
from __future__ import print_function

import glob
import hashlib
import json
import os
import re
import subprocess
import sys
from os.path import join as osj

def find_any_samplesheet(runDir, fromResDir = False):
	"""
	Adapted from runmetrics.py

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

def getMd5(run):
	runMd5 = hashlib.md5()
	with open(run, "rb") as f:
		for chunk in iter(lambda: f.read(4096), b""):
			runMd5.update(chunk)
	return runMd5.hexdigest()

def startCanoes(run, serviceName, serviceDockerImage):
	samplesheet = find_any_samplesheet(run)
	with open(osj(run,"analysis.json"),'r') as jsonAnalysisFile:
			jsonData = json.load(jsonAnalysisFile)
	bed = jsonData["analysis"]["bed"]
	genome = jsonData["analysis"]["genome"]
	# md5 =
	cmd = "docker run -dti --name=service-"+serviceName+"-"+md5+" -v "+run+":"+run+" -v "+genome+":"+genome+" "+serviceDockerImage+" /bin/bash /app/bin/canoes run -r "+run+" -s "+samplesheet+" -l "+bed+" -e A,F,M -g "+genome+" -o "+run+"/CANOES"
	subprocess.call(cmd, shell=True)
	cmd = [ "$( docker container inspect -f '{{.State.Status}}' service-canoes )" == "running" ]
	if subprocess.call(cmd, shell=True) == True:
			return True
	return False

def complete(run,serviceName):
	if glob.glob(osj(run,serviceName+"Complete.txt")):
			return True
	return False

def checkTags(serviceTag, sampleTagsList, analysisTagsList):
	for t in sampleTagsList:
			if serviceTag in t:
					return True
	for t in analysisTagsList:
			if serviceTag in t:
					return True
	return False

def getAnalysisTagsFromSampleList(sampleList,run):
	tagsList = []
	for s in sampleList:
			with open (osj(run,s,"STARK",s+".analysis.tag"), "r") as tagsFile:
					for l in tagsFile:
							tagsList.append(l.strip())
	return tagsList

def getSampleTagsFromSampleList(sampleList,run):
	tagsList = []
	for s in sampleList:
			with open (osj(run,s,"STARK",s+".tag"), "r") as tagsFile:
					for l in tagsFile:
							tagsList.append(l.strip())
	return tagsList

def getSampleListFromRunPath(run):
	"""
	Only returns folders with a SampleSheet
	"""
	sampleList = []
	for patient in glob.glob(osj(run,"*/")):
			pName = os.path.basename(os.path.normpath(patient))
			if os.path.exists(osj(patient, "STARK", pName+".SampleSheet.csv")) or os.path.exists(osj(patient, pName+".SampleSheet.csv")):
					sampleList.append(pName)
	return sampleList

def getRunName(starkComplete):
	assert starkComplete.endswith('/STARKComplete.txt'), "[ERROR]"
	if starkComplete.endswith('/STARKComplete.txt'):
			return starkComplete[:-18]

def main(starkComplete, serviceTag, serviceName):
	run = getRunName(starkComplete)
	sampleList = getSampleListFromRunPath(run)
	sampleTagsList = getSampleTagsFromSampleList(sampleList, run)
	analysisTagsList = getAnalysisTagsFromSampleList(sampleList, run)
	if not checkTags(serviceTag, sampleTagsList, analysisTagsList):
			sys.exit()
	if complete(run,serviceName):
			sys.exit()
	print(startCanoes(run, serviceName, serviceDockerImage))

if __name__ == "__main__":
	#Path for BAS use
	starkComplete = "/home1/BAS/grentziv/Runs/input/200114_NB551027_0671_AHK2C5AFXY/STARKComplete.txt"
	serviceTag = ""
	serviceName = "CANOES"
	serviceDockerImage = "canoes:1.1"
	main(starkComplete, serviceTag, serviceName)