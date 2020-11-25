#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import hashlib
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

def getMd5(run):
	runMd5 = hashlib.md5()
	runMd5.update(hashlib.md5(run).hexdigest())
	return runMd5.hexdigest()

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

def createRunningFile(run, serviceName):
	file = open(osj(run,serviceName+"Running.txt"), "w+")
	file.write("# ["+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+"] "+os.path.basename(run)+" running with "+serviceName+"\n")
	file.close()

def launch(run, serviceName, containersFile, montage, image, launchCommand, configFile, repository):
	createRunningFile(run, serviceName)
	samplesheet = find_any_samplesheet(run)
	assert samplesheet != "NO_SAMPLESHEET_FOUND",\
		"[ERROR] find_any_samplesheet() couldn't find any samplesheet in run"+run+"."
	bed = findAnyBed(run)
	md5 = getMd5(run)
	containerName = serviceName+"-"+md5+"-NAME-"+os.path.basename(run)
	#cmd = "docker run --rm --name="+containerName+" -v "+run+":"+run+" "+montage+" "+image+" "+launchCommand+" --run="+run+"'"
	cmd = "docker run --rm --name="+containerName+" "+montage+" "+image+" "+launchCommand+" --run="+run+"'"
	subprocess.call(cmd, shell = True)
	createContainerFile(containersFile, run, containerName)
