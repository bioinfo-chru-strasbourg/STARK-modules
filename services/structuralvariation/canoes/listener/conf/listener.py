#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Goal: 
- Listener for STARK18 modules
- Launch analysis if triggers are passed

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
import time

from datetime import datetime
from launcher import launch
from os.path import join as osj

def createContainerFile(containersFile, run, containerName):
	file = open(osj(containersFile,containerName+".log"), "w+")
	file.write("RUN: "+os.path.basename(run)+"\n")
	file.write("FOLDER: "+run+"\n")
	file.write("EXEC_DATE: "+datetime.now().strftime("%d%m%Y-%H%M%S")+"\n")
	file.write("ID: "+containerName+"\n")
	file.close()

def getMd5(run):
	runMd5 = hashlib.md5()
	runMd5.update(hashlib.md5(run).hexdigest())
	return runMd5.hexdigest()

def createRunningFile(run, serviceName):
	file = open(osj(run,serviceName+"Running.txt"), "w+")
	file.write("# ["+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+"] "+os.path.basename(run)+" running with "+serviceName+"\n")
	file.close()

def checkProjects(serviceProject, sampleProjectsList):
	listProject = []
	for project in sampleProjectsList:
		listProject.append(project.split("-")[0])
	for p in listProject:
		if serviceProject in p:
			return True
	return False

def checkGroup(serviceGroup, sampleGroupsList):
	# print(serviceGroup)
	listGroup = []
	for group in sampleGroupsList:
		listGroup.append(group.split("-")[0])
	# print(listGroup)
	for g in listGroup:
		if serviceGroup in g:
			return True
	return False

def assert_file_exists_and_is_readable(filePath):
	assert os.path.isfile(filePath) and os.access(filePath, os.R_OK), \
			"[ERROR] File "+filePath+" doesn't exist or isn't readable"

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
	i=0
	with open(samplesheetPath, "r") as f:
		for l in f:
			if not inDataTable:
				if l.startswith("Sample_ID,Sample_Name,"):
					inDataTable = True
					sampleProjectIndex = l.strip().split(",").index("Sample_Project")
			else:
				if "," in l:
					sampleProject.append(l.strip().split(",")[sampleProjectIndex])
	while i < len(sampleProject):
		if sampleProject[i] == "":
			del sampleProject[i]
		else:
			i+=1
	if not sampleProject:
		with open(samplesheetPath, "r") as f:
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

def complete(run, serviceName):
	if glob.glob(osj(run,serviceName+"Complete.txt")):
		return False
	return True

def running(run, serviceName):
	if glob.glob(osj(run,serviceName+"Running.txt")):
		return False
	return True

def checkTriggers(jconfig, serviceName, run):
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
			# print(listNotFile)
			# print(listFile)
			if all(listNotFile + listFile):
				return True
			else:
				return False
		elif key == "tags":
			listTag = []
			for tag in jconfig[key]:
				listTag.append(checkTags(tag, getSampleTagsFromSampleList(getSampleListFromRunPath(run), run), getAnalysisTagsFromSampleList(getSampleListFromRunPath(run), run)))
			# print(listTag)
			if all(listTag):
				return True
			else:
				return False
		elif key == "group" or key == "project":
			if not sampleProject:
				sampleProject = get_sample_project_from_samplesheet(find_any_samplesheet(run))
				# print(sampleProject)
			if key == "group":
				listNotGroup = []
				listGroup = []
				for group in jconfig[key]:
					if group.startswith("!"):
						listNotGroup.append(not checkGroup(group[1:], sampleProject))
					else:
						listGroup.append(checkGroup(group, sampleProject))
				# print(listNotGroup)
				# print(listGroup)
				if any(listGroup) and all(listNotGroup):
					return True
				else:
					return False
			elif key == "project":
				listProject = []
				for project in jconfig[key]:
					listProject.append(checkProjects(project, sampleProject))
				# print(listProject)
				if any(listProject):
					return True
				else:
					return False

def getDataFromJson(jsonFile):
	with open(jsonFile,'r') as jsonAnalysisFile:
		jsonData = json.load(jsonAnalysisFile)
	return jsonData

def checkTriggersFromConfig(jsonFile, serviceName, run):
	jconfig = getDataFromJson(jsonFile)['services'][serviceName]['triggers']
	return checkTriggers(jconfig, serviceName, run)

def verifyTriggers(run, serviceName, jsonFile):
	return checkTriggersFromConfig(jsonFile, serviceName, run)

def getRunName(starkComplete):
	assert starkComplete.endswith('/STARKCopyComplete.txt'), "[ERROR]"
	if starkComplete.endswith('/STARKCopyComplete.txt'):
		return starkComplete[:-22]

def getStarkCopyComplete(groupInput, days):
	starkCompleteList = []
	p = subprocess.Popen("find "+groupInput+"/*/STARKCopyComplete.txt -mtime -"+str(days)+" 2>/dev/null", stdout=subprocess.PIPE, shell=True)
	out = p.stdout.readlines()
	if not out:
		p = subprocess.Popen("find "+groupInput+"/*/*/STARKCopyComplete.txt -mtime -"+str(days)+" 2>/dev/null", stdout=subprocess.PIPE, shell=True)
		out = p.stdout.readlines()
		if not out:
			p = subprocess.Popen("find "+groupInput+"/*/*/*/STARKCopyComplete.txt -mtime -"+str(days)+" 2>/dev/null", stdout=subprocess.PIPE, shell=True)
			out = p.stdout.readlines()
	for complete in out:
		complete = complete.decode("utf-8").strip()
		starkCompleteList.append(complete)
	return starkCompleteList

def main(groupInputList, serviceName, jsonFile, days, delay, containersFile, configFile):
	while True:
		for groupInput in groupInputList:
			starkCompleteList = getStarkCopyComplete(groupInput, days)
			for starkComplete in starkCompleteList:
				run = getRunName(starkComplete)
				if verifyTriggers(run, serviceName, jsonFile):
					print("Launching "+serviceName+" analysis for run "+run)
					launch(run, serviceName, containersFile, os.getenv('MICROSERVICE_MONTAGE'), os.getenv('MICROSERVICE_IMAGE'), os.getenv('MICROSERVICE_LAUNCH'), configFile)
		time.sleep(60.0*delay)

def myoptions():
	'''
	*arg parser*
	*return options*
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", type = str, default = " ", help = "list of runs to listen to : <PATH_RUN_1>,<PATH_RUN_2>...", dest = 'listInput')
	parser.add_argument("-n", "--servicename", type = str, default = "TEST", help = "name of the service to use", dest = 'serviceName')
	parser.add_argument("-j", "--json", type = str, default = "listener.json", help = "path to the triggers file listener.json", dest = 'jsonFile')
	parser.add_argument("-c", "--config", type = str, default = "listener.conf", help = "path to the config file listener.conf", dest = 'configFile')
	parser.add_argument("-f", "--containersfile", type = str, default = "", help = "path to the container's file folder", dest = 'containersFile')
	parser.add_argument("-t", "--nbDaysBack", type = int, default = 30, help = "folder older than x days", dest = 'days')
	parser.add_argument("-d", "--minDelay", type = int, default = 5, help = "delay to sleep between two listening", dest = 'delay')
	return parser.parse_args()

if __name__ == "__main__":
	args = myoptions()
	groupInputList = args.listInput.split(",")
	main(groupInputList, args.serviceName, args.jsonFile, args.days, args.delay, args.containersFile, args.configFile)

