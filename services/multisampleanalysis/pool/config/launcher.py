#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import glob
import json
import os
import re
import subprocess

from datetime import datetime
from os.path import join as osj

def createContainerFile(containersFile, run, containerName):
	file = open(osj(containersFile,containerName+".log"), "w+")
	file.write("RUN: "+os.path.basename(run)+"\n")
	file.write("FOLDER: "+run+"\n")
	file.write("EXEC_DATE: "+datetime.now().strftime("%d%m%Y-%H%M%S")+"\n")
	file.write("ID: "+containerName+"\n")
	file.close()

def createRunningFile(run, serviceName):
	file = open(osj(run,serviceName+"Running.txt"), "w+")
	file.write("# ["+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+"] "+os.path.basename(run)+" running with "+serviceName+"\n")
	file.close()

def checkAllSamples(sampleList, run):
	sampleFound = []
	for sample in sampleList:
		if glob.glob(osj(run,sample,"STARKCopyComplete.txt")):
			sampleFound.append(True)
		else:
			sampleFound.append(False)
	if all(sampleFound):
		return True
	else:
		return False

def get_sample_id_from_samplesheet(samplesheetPath):
	"""
	Adapted from functions.py
	
	Returns a python list containing all sample ID from samplesheet.
	"""
	if samplesheetPath == "NO_SAMPLESHEET_FOUND" or not samplesheetPath:
		return []
	inDataTable = False
	sampleID = []
	analysis_json = samplesheetPath.replace("SampleSheet.csv", "analysis.json")
	with open(analysis_json, "r") as jf:
		application = json.load(jf)["application"][0].split("+")[0]
	with open(samplesheetPath, "r") as f:
		for l in f:
			if not inDataTable:
				if l.startswith("Sample_ID,"):
					inDataTable = True
					sampleIDIndex = l.strip().split(",").index("Sample_ID")
			else:
				if "," in l:
					if (
                        "APP" in l.strip().split(",")[-1]
                        and application not in l.strip().split(",")[-1]
                    ):
						continue
					else:
						sampleID.append(l.strip().split(",")[sampleIDIndex])
	return sampleID

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

def launch(run, serviceName, containersFile, montage, image, launchCommand, configFile, repository):
	samplesheet = find_any_samplesheet(run)
	assert samplesheet != "NO_SAMPLESHEET_FOUND",\
		"[ERROR] find_any_samplesheet() couldn't find any samplesheet in run"+run+"."
	if not checkAllSamples(get_sample_id_from_samplesheet(samplesheet), run):
		return "Not all sample in directory"
	createRunningFile(run, serviceName)
	containerName = serviceName+"-NAME-"+os.path.basename(run)
	if run.startswith("/STARK/"):
		outerRunPath = run.replace("/STARK/output/repository", repository)
		cmd = "docker run --rm --name="+containerName+" -v "+outerRunPath+":"+run+" "+montage+" "+image+" "+launchCommand+" -i "+run
		print(cmd)
		subprocess.call(cmd, shell = True)
	else:
		cmd = "docker run --rm --name="+containerName+" -v "+run+":"+run+" "+montage+" "+image+" "+launchCommand+" -i "+run
		print(cmd)
		subprocess.call(cmd, shell = True)
	createContainerFile(containersFile, run, containerName)