#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Goal: 
- Listener for CANOES module for STARK18
- Launch CANOES module in a Docker container

@Author: Victor Grentzinger (2020)
"""

from __future__ import division
from __future__ import print_function

import argparse
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
import time

from datetime import datetime
from os.path import join as osj

def assert_file_exists_and_is_readable(filePath):
	assert os.path.isfile(filePath) and os.access(filePath, os.R_OK), \
			"[ERROR] File "+filePath+" doesn't exist or isn't readable"

def checkTriggersFromConfig(config, serviceName, run):
	jconfig = getDataFromJson(config)['services'][serviceName]['triggers']
	return checkTriggers(jconfig, serviceName, run)

def checkTriggers(jconfig, serviceName, run):
	# print(run)
	if len(jconfig.keys()) == 0:
		return True
	sampleProject = []
	for key in jconfig.keys():
		if key == "AND":
			if all([checkTriggers({andKey : jconfig[key][andKey]}, serviceName, run) for andKey in jconfig[key].keys()]):
				return True
			else: 
				return False
		elif key.startswith("OR_"):
			if any([checkTriggers({andKey : jconfig[key][andKey]}, serviceName, run) for andKey in jconfig[key].keys()]):
				return True
			else: 
				return False
		elif key == "file":
			listNotFile = []
			listFile = []
			for file in jconfig[key]:
				if file.startswith("!"):
					if "Running.txt" in file:
						listNotFile.append(running(run, file[1:-11]))
					elif "Complete.txt" in file:
						listNotFile.append(complete(run, file[1:-12]))
					elif "STARKCopyComplete.txt" in file:
						listNotFile.append(complete(run, "STARK"))
				else:
					if "Running.txt" in file:
						listFile.append(not running(run, file[:-11]))
					elif "Complete.txt" in file:
						listFile.append(not complete(run, file[:-12]))
					elif "STARKCopyComplete.txt" in file:
						listFile.append(not complete(run, "STARK"))
			if all(listNotFile + listFile):
				# print("Files ok")
				return True
			else:
				# print("Files not ok")
				return False
		elif key == "tags":
			listTag = []
			# if all([checkTags(tag, getSampleTagsFromSampleList(getSampleListFromRunPath(run), run), getAnalysisTagsFromSampleList(getSampleListFromRunPath(run), run)) for tag in jconfig[key]]):
			for tag in jconfig[key]:
				listTag.append(checkTags(tag, getSampleTagsFromSampleList(getSampleListFromRunPath(run), run), getAnalysisTagsFromSampleList(getSampleListFromRunPath(run), run)))
			if all(listTag):
				# print("Tags ok")
				return True
			else:
				# print("Tags not ok")
				return False
		elif key == "group" or key == "project":
			if not sampleProject:
				sampleProject = get_sample_project_from_samplesheet(find_any_samplesheet(run))
			if key == "group":
				listGroup = []
				for group in jconfig[key]:
					listGroup.append(checkGroup(group, sampleProject))
				if any(listGroup):
					return True
				else:
					return False
			elif key == "project":
				listProject = []
				for project in jconfig[key]:
					listProject.append(checkProjects(project, sampleProject))
				if any(listProject):
					return True
				else:
					return False

def get_sample_project_from_samplesheet(samplesheetPath):
	"""
	Adapted from functions.py
	
	Returns a python list containing all tags from samplesheet.
	"""
	assert samplesheetPath != "NO_SAMPLESHEET_FOUND", \
			"[ERROR] find_any_samplesheet() couldn't find any samplesheet. Check if the --fromResultDir argument is set correctly."
	assert_file_exists_and_is_readable(samplesheetPath)
	inDataTable = False
	sampleProject = []
	with open(samplesheetPath, "r") as f:
		for l in f:
			if not inDataTable:
				if l.startswith("Sample_ID,Sample_Name,"):
					inDataTable = True
					sampleProjectIndex = l.strip().split(",").index("Sample_Project")
			else:
				if "," in l:
					sampleProject.append(l.strip().split(",")[sampleProjectIndex])
		if not sampleProject:
			for l in f:
				if l.startswith("Investigator Name"):
					sampleProject.append(l.strip().split(",")[1])
	return sampleProject

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

def findAnyBed(run):
	p = subprocess.Popen("find -L "+run+" -maxdepth 3 -name *.bed", stdout=subprocess.PIPE, shell=True)
	out = p.stdout.readlines()
	for bed in out:
		bed = bed.decode("utf-8").strip()
		r = re.match(run.rstrip("/")+"/(.*)/STARK/(.*).bed", bed)
		if r is None:
			continue
		elif r.group(1) == r.group(2): #checks if (.*) == (.*)
			return bed
	return "NO_BED_FOUND"

def createContainerFile(config, run, serviceName):
	if os.path.exists(osj(os.path.dirname(config),"ContainerList"+serviceName+".txt")):
		file = open(osj(os.path.dirname(config),"ContainerList"+serviceName+".txt"), "a")
		file.write(run+"\n")
		file.close()
	else:
		file = open(osj(os.path.dirname(config),"ContainerList"+serviceName+".txt"), "w+")
		file.write(run+"\n")
		file.close()

def createRunningFile(run, serviceName):
	file = open(osj(run,serviceName+"Running.txt"), "w+")
	file.write("# ["+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+"] "+os.path.basename(run)+" running with "+serviceName+"\n")
	file.close()

def getDataFromJson(file):
	with open(file,'r') as jsonAnalysisFile:
		jsonData = json.load(jsonAnalysisFile)
	return jsonData

def getMd5(run):
	runMd5 = hashlib.md5()
	runMd5.update(hashlib.md5(run).hexdigest())
	return runMd5.hexdigest()

def startService(run, serviceName, serviceDockerImage, annotsvServer, annotsvContainer, genome, config):
	samplesheet = find_any_samplesheet(run)
	assert samplesheet != "NO_SAMPLESHEET_FOUND",\
		"[ERROR] find_any_samplesheet() couldn't find any samplesheet in run"+run+"."
	bed = findAnyBed(run)
	assert bed != "NO_BED_FOUND",\
		"[ERROR] findAnyBed() couldn't find any bed in run"+run+"."
	md5 = getMd5(run)
	containerName = "service-"+serviceName+"-"+md5
	if serviceName == "CANOES":
		cmd = "docker run --rm -dti --name="+containerName+" -v "+run+":"+run+" -v "+genome+":"+genome+" -v "+annotsvServer+":"+annotsvContainer+" "+serviceDockerImage+" /bin/bash /app/bin/canoes run -r "+run+" -s "+samplesheet+" -l "+bed+" -g "+genome+"/hg19.fa -o "+run
	subprocess.call(cmd, shell = True)
	createRunningFile(run, serviceName)
	createContainerFile(config, run, serviceName)

def running(run, serviceName):
	if glob.glob(osj(run,serviceName+"Running.txt")):
		return False
	return True

def complete(run, serviceName):
	if glob.glob(osj(run,serviceName+"Complete.txt")):
		return False
	return True

def checkGroup(serviceGroups, sampleGroupsList):
	listGroup = []
	for group in sampleGroupsList:
		listGroup.append(group.split("-")[0])
	for serviceGroup in serviceGroups:
		for g in listGroup:
			if serviceGroup in g:
				return True
	return False

def checkProjects(serviceProjects, sampleProjectsList):
	listProject = []
	for project in sampleProjectsList:
		listProject.append(project.split("-")[0])
	for serviceProject in serviceProjects:
		for p in listProject:
			if serviceProject in p:
				return True
	return False

def checkTags(tag, sampleTagsList, analysisTagsList):
	for serviceTag in tag:
		if serviceTag.startswith("!"):
			for t in sampleTagsList:
				if not serviceTag in t:
					return True
			for t in analysisTagsList:
				if not serviceTag in t:
					return True
			return False
		else:
			for t in sampleTagsList:
				if serviceTag in t:
					return True
			for t in analysisTagsList:
				if serviceTag in t:
					return True
			return False

def getAnalysisTagsFromSampleList(sampleList,run):
	analysisTagsList = []
	for s in sampleList:
		with open (osj(run,s,s+".analysis.tag"), "r") as tagsFile:
			for l in tagsFile:
				analysisTagsList.append(l.strip())
	return analysisTagsList

def getSampleTagsFromSampleList(sampleList,run):
	sampleTagsList = []
	for s in sampleList:
		with open (osj(run,s,s+".tag"), "r") as tagsFile:
			for l in tagsFile:
				sampleTagsList.append(l.strip())
	return sampleTagsList

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
	assert starkComplete.endswith('/STARKCopyComplete.txt'), "[ERROR]"
	if starkComplete.endswith('/STARKCopyComplete.txt'):
		return starkComplete[:-22]

def getStarkCopyComplete(groupInput, days):
	starkCompleteList = []
	p = subprocess.Popen("find "+groupInput+"/*/*/*/STARKCopyComplete.txt -mtime -"+str(days)+" 2>/dev/null", stdout=subprocess.PIPE, shell=True)
	out = p.stdout.readlines()
	for complete in out:
		complete = complete.decode("utf-8").strip()
		starkCompleteList.append(complete)
	return starkCompleteList

def main(groupInputList, serviceName, config, serviceDockerImage, annotsvServer, annotsvContainer, minDelay, days, genome):
	while True:
		for groupInput in groupInputList:
			starkCompleteList = getStarkCopyComplete(groupInput, days)
			for starkComplete in starkCompleteList:
				run = getRunName(starkComplete)
				if checkTriggersFromConfig(config, serviceName, run):
					startService(run, serviceName, serviceDockerImage, annotsvServer, annotsvContainer, genome, config)
		time.sleep(60.0*minDelay)

def myoptions():
	'''
	*arg parser*
	*return options*
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", type = str, help = "list of group to listen to : <PATH_RUN_1>,<PATH_RUN_2>...", dest = 'groupListeningInput', required = True)
	parser.add_argument("-n", "--servicename", type = str, help = "name of the service to use", dest = 'serviceName', required = True)
	parser.add_argument("-d", "--dockerimage", type = str, help = "docker image name of the service", dest = 'serviceDockerImage', required = True)
	parser.add_argument("-p", "--config", type = str, default = "/home1/data/STARK_09181/config/services.json", help = "path to the config file services.json", dest = 'config')
	parser.add_argument("-S", "--annotsvserver", type = str, default = "/home1/TOOLS/tools/AnnotSV/AnnotSV_2.1", help = "path to the AnnotSV_2.1 folder on the server", dest = 'annotsvServer')
	parser.add_argument("-c", "--annotsvcontainer", type = str, default = "/src/AnnotSV/2.1", help = "path to the AnnotSV_2.1 folder in the container", dest = 'annotsvContainer')
	parser.add_argument("-D", "--minDelay", type = int, default = 5, help = "daemon refresh delay", dest = 'minDelay')
	parser.add_argument("-t", "--nbDaysBack", type = int, default = 30, help = "folder older than x days", dest = 'days')
	parser.add_argument("-g", "--genome", type = str, default = "/home1/data/STARK_09181/databases/genomes/current", help = "path to genome folder", dest = 'genome')
	return parser.parse_args()

if __name__ == "__main__":
	#Path for HUX234 use
	args = myoptions()
	groupInputList = args.groupListeningInput.split(",")
	main(groupInputList, args.serviceName, args.config, args.serviceDockerImage, args.annotsvServer, args.annotsvContainer, args.minDelay, args.days, args.genome)