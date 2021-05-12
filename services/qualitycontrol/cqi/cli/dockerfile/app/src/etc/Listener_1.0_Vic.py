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

def checkTriggersFromConfig(services, serviceName, run):
	config = getDataFromJson(services)['services'][serviceName]['triggers']
	return checkTriggers(config, serviceName, run)

def checkTriggers(config, serviceName, run):
	if len(config.keys()) == 0:
		return True
	sampleProject = []
	for key in config.keys():
		if key == "AND":
			if all([checkTriggers({andKey : config[key][andKey]}, serviceName, run) for andKey in config[key].keys()]):
				return True
			else: 
				return False
		elif key.startswith("OR_"):
			if any([checkTriggers({andKey : config[key][andKey]}, serviceName, run) for andKey in config[key].keys()]):
				return True
			else: 
				return False
		elif key == "file":
			listNotFile = []
			listFile = []
			for file in config[key]:
				if file.startswith("!"):
					if "Running.txt" in file:
						listNotFile.append(running(run, file[1:-11]))
					elif "Complete.txt" in file:
						listNotFile.append(complete(run, file[1:-12]))
					elif "STARKComplete.txt" in file:
						listNotFile.append(complete(run, "STARK"))
				else:
					if "Running.txt" in file:
						listFile.append(not running(run, file[:-11]))
					elif "Complete.txt" in file:
						listFile.append(not complete(run, file[:-12]))
					elif "STARKComplete.txt" in file:
						listFile.append(not complete(run, "STARK"))
			print(listNotFile + listFile)
			print(key,config[key])
			if all(listNotFile + listFile):
				return True
			else:
				return False
		elif key == "tags":
			print([checkTags(tag, getSampleTagsFromSampleList(getSampleListFromRunPath(run), run), getAnalysisTagsFromSampleList(getSampleListFromRunPath(run), run)) for tag in config[key]])
			print(key,config[key])
			if all([checkTags(tag, getSampleTagsFromSampleList(getSampleListFromRunPath(run), run), getAnalysisTagsFromSampleList(getSampleListFromRunPath(run), run)) for tag in config[key]]):
				return True
			else:
				return False
		elif key == "group" or key == "project":
			if not sampleProject:
				sampleProject = get_sample_project_from_samplesheet(find_any_samplesheet(run))
			if key == "group":
				print([checkGroup(group, sampleProject) for group in config[key]])
				print(key,config[key])
				if any([checkGroup(group, sampleProject) for group in config[key]]):
					return True
				else:
					return False
			elif key == "project":
				print([checkProjects(project, sampleProject) for project in config[key]])
				print(key,config[key])
				if any([checkProjects(project, sampleProject) for project in config[key]]):
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



def removeContainer(containerName, run, serviceName):
	for name in containerName:
		if subprocess.Popen("docker container inspect -f '{{.State.Status}}' "+name, stdout=subprocess.PIPE, shell=True).communicate()[0][:-1] == 'running':
			continue
		else:
			subprocess.Popen("docker rm "+name, stdout=subprocess.PIPE, shell = True)
			createCompleteFile(run, serviceName)
			del(containerName[containerName.index(name)])

def createCompleteFile(run, serviceName):
	file = open(osj(run,serviceName+"Complete.txt"), "w+")
	file.write("# ["+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+"] "+os.path.basename(run)+" analyzed with "+serviceName)
	file.close()

def createRunningFile(run, serviceName):
	file = open(osj(run,serviceName+"Running.txt"), "w+")
	file.write("# ["+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+"] "+os.path.basename(run)+" running with "+serviceName)
	file.close()

def getDataFromJson(file):
	with open(file,'r') as jsonAnalysisFile:
		jsonData = json.load(jsonAnalysisFile)
	return jsonData

def getMd5(run):
	runMd5 = hashlib.md5()
	runMd5.update(hashlib.md5(run).hexdigest())
	return runMd5.hexdigest()

def startService(run, serviceName, serviceDockerImage, annotsvServer, annotsvContainer):
	samplesheet = find_any_samplesheet(run)
	assert samplesheet != "NO_SAMPLESHEET_FOUND",\
		"[ERROR] find_any_samplesheet() couldn't find any samplesheet."
	bed = getDataFromJson(osj(run,"analysis.json"))["analysis"]["bed"]
	genome = getDataFromJson(osj(run,"analysis.json"))["analysis"]["genome"][:-7]
	md5 = getMd5(run)
	containerName = "service-"+serviceName+"-"+md5
	if serviceName == "CANOES":
		cmd = "docker run -dti --name="+containerName+" -v "+run+":"+run+" -v "+genome+":"+genome+" -v "+annotsvServer+":"+annotsvContainer+" "+serviceDockerImage+" /bin/bash /app/bin/canoes run -r "+run+" -s "+samplesheet+" -l "+bed+" -e A -g "+genome+"/hg19.fa -o "+run
	subprocess.call(cmd, shell = True)
	createRunningFile(run, serviceName)
	return containerName

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

def checkTags(serviceTags, sampleTagsList, analysisTagsList):
	for serviceTag in serviceTags:
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

def getStarkComplete(groupInput, days):
	starkCompleteList = []
	p = subprocess.Popen("find "+groupInput+"/*/STARKComplete.txt -mtime -"+str(days)+" 2>/dev/null", stdout=subprocess.PIPE, shell=True)
	out = p.stdout.readlines()
	for complete in out:
		complete = complete.decode("utf-8").strip()
		starkCompleteList.append(complete)
	return starkCompleteList

def main(groupInputList, serviceTags, serviceName, services, serviceDockerImage, annotsvServer, annotsvContainer, minDelay, days):
	containerName = []
	while True:
		for groupInput in groupInputList:
			starkCompleteList = getStarkComplete(groupInput, days)
			for starkComplete in starkCompleteList:
				run = getRunName(starkComplete)
				if checkTriggersFromConfig(services, serviceName, run):
					containerName.append(startService(run, serviceName, serviceDockerImage, annotsvServer, annotsvContainer))
		# Launch loop if list "containerName" is not empty
		if containerName:
			removeContainer(containerName, run, serviceName)
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
	parser.add_argument("-T", "--tags", type = str, help = "tags of services", dest = 'serviceTagsInput', required = True)
	parser.add_argument("-s", "--servicesfile", type = str, default = "/home1/BAS/grentziv/CANOES/1.1.0/services.json", help = "path to the services.json configfile", dest = 'services')
	parser.add_argument("-S", "--annotsvserver", type = str, default = "/home1/TOOLS/tools/AnnotSV/AnnotSV_2.1", help = "path to the AnnotSV_2.1 folder on the server", dest = 'annotsvServer')
	parser.add_argument("-c", "--annotsvcontainer", type = str, default = "/src/AnnotSV/2.1", help = "path to the AnnotSV_2.1 folder in the container", dest = 'annotsvContainer')
	parser.add_argument("-D", "--minDelay", type = int, default = 5, help = "daemon refresh delay", dest = 'minDelay')
	parser.add_argument("-t", "--nbDaysBack", type = int, default = 30, help = "folder older than x days", dest = 'days')
	return parser.parse_args()

if __name__ == "__main__":
	#Path for BAS use
	args = myoptions()
	groupInputList = args.groupListeningInput.split(",")
	serviceTags = args.serviceTagsInput.split("!")
	main(groupInputList, serviceTags, args.serviceName, args.services, args.serviceDockerImage, args.annotsvServer, args.annotsvContainer, args.minDelay, args.days)