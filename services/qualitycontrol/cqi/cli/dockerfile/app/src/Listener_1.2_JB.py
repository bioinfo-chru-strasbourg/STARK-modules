#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Goal: 
- Listener for CQI module for STARK18
- Launch CQI module in a Docker container

@Author: Victor Grentzinger (2020)
@Addapted by Jean-Baptiste Lamouche (2020)
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
	#underTag(sampleList, run)

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



def removeContainer(md5, run, serviceName, groupInput):
	if subprocess.Popen("docker container inspect -f '{{.State.Status}}' service-"+serviceName+"-"+md5, stdout=subprocess.PIPE, shell=True).communicate()[0][:-1] == 'exited':
		subprocess.Popen("docker rm service-"+serviceName+"-"+md5, stdout=subprocess.PIPE, shell = True)
		subprocess.Popen("sed -i '/"+os.path.basename(os.path.normpath(run))+"/d' "+osj(groupInput,"ContainerList"+serviceName+".txt"), stdout=subprocess.PIPE, shell = True)

def createContainerFile(groupInput, run, serviceName):
	if os.path.exists(osj(groupInput,"ContainerList"+serviceName+".txt")):
		file = open(osj(groupInput,"ContainerList"+serviceName+".txt"), "a")
		file.write(run+"\n")
		file.close()
	else:
		file = open(osj(groupInput,"ContainerList"+serviceName+".txt"), "w+")
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
	runMd5 = run.encode("utf-8") 
	return hashlib.md5(runMd5).hexdigest()

#def underTag(sampleList, run):
#	print(getSampleTagsFromSampleList(sampleList, run))

def parseTag(run):
	samplesheet = find_any_samplesheet(run)
	tag = {}
	tmp = []
	with open(samplesheet) as origin_file:
		for line in origin_file:
			if re.match("(.*)Description,(.*)", line):
				print(line)
				if line:
					line = line.split(',')
					for elem in line:
						if elem != "Description":
							tmp.append(elem.strip())
	
	for item in tmp:
		item = item.split('!')
		if item:
			for x in item:
				x = x.split('#')
				tag[x[0]] = x[1:]
	
	return tag


def startService(run, serviceName, serviceDockerImage, JSON, env, cqiFolder, groupInput, genome, serviceDockerImageRtg, cqiBin):
	samplesheet = find_any_samplesheet(run)
	assert samplesheet != "NO_SAMPLESHEET_FOUND", \
				"[ERROR] find_any_samplesheet() couldn't find any samplesheet."
	bed = findAnyBed(run)
	tag = parseTag(run)
	assert bed != "NO_BED_FOUND", \
				"[ERROR] findAnyBed() couldn't find any bed."
	md5 = getMd5(run)
	containerName = "service-"+serviceName+"-"+md5
	
	#tmp = []
	#tag = []
	#with open(samplesheet) as origin_file:
	#	for line in origin_file:
	#		if re.match("(.*)#CQI(.*)", line):
	#			if line:
	#				line = line.split('!')
	#				for elem in line:
	#					tmp.append(elem.strip())
	#for item in tmp:
	#	item = item.split('#')
	#	if 'CQI' in item:
	#		tag.append(item)
	#
	#print(tag)
	
	#V1
	#ref = run+"/list_vcf_ref.txt"
	#vcfRef = []
	#if len(tag) <= 1:
	#	for elem in tag:
	#		vcfRef.append(parseJsonCQI(JSON, tag[0][2]).strip())
	#	vcfStark = tag[0][2] #it will be a CQI sample name ex: "HORIZON" "TeCoriell" etc
	#	print(vcfRef)
	#	print(vcfStark)
	#	with open(ref, 'w') as r:
	#		for item in vcfRef:
	#			r.write("%s\n" % item)
	#elif len(tag) > 1 and tag[0][1] =="CQI":
	#	for item in tag:
	#		vcfRef.append(parseJsonCQI(JSON, tag[0][2]).strip()) 
	#	vcfStark = tag[0][2]
	#	with open(ref, 'w') as r:
	#		for item in vcfRef:
	#			r.write("%s\n" % item)
	##useless cuz in nonregression we directly call compare scripts we do not use Listener
	#elif len(tag) == 0:
	#	ref = None
	#	vcfStark = None
		
	#V2
	vcfRef = []
	ref = run+"/list_vcf_ref.txt"
	for key, value in tag.items():
		if key == 'CQI':
			for sample in value:
				vcfRef.append(parseJsonCQI(JSON, sample))
	print(vcfRef)
	with open(ref, 'w') as r:
		for item in vcfRef:
			r.write("%s\n" % item)		
	
		
	#with open(run+'/list_vcf_stark.txt', 'w') as stark:
	#	for v in vcfStark:
	#		stark.write("%s\n" % v)

	#if serviceName == "CQI":
	cmd = "docker run -ti --rm --name="+containerName+" --env-file="+env+" -v "+genome+":"+genome+" -v "+cqiFolder+":"+cqiFolder+" -v "+run+":"+run+" -v "+cqiBin+":"+cqiBin+" "+serviceDockerImage+" bash -c 'source activate variant && /home1/BAS/lamouchj/CQI/bin/compare_vcf_1.1.sh --run="+run+" --genes="+bed+" --cqi="+vcfStark+" --cqigz="+ref+"'"
	print("[#INFO] Start of CQI analysis; Container launch")
	subprocess.call(cmd, shell = True)
	createRunningFile(run, serviceName)
	createContainerFile(groupInput, run, serviceName)
	while True:
		if os.path.isfile(run+'/'+'COMPAREComplete.txt'):
			name = "service-"+serviceName+"-RTG-"+md5
			rtg = "docker run -ti --rm --name="+name+" --entrypoint"+" /bin/bash"+" --env-file="+env+" -v "+run+":"+run+" -v "+genome+":"+genome+" -v "+cqiBin+":"+cqiBin+" "+serviceDockerImageRtg+" /home1/BAS/lamouchj/CQI/bin/rtg_vcf.sh --cqi="+tag[0][2]+" --run="+run
			subprocess.call(rtg, shell=True)
			print("[#INFO] Launch RTG container")
			break
		else:
			time.sleep(60.0)
		
#def startServiceSpe(run, serviceName, serviceDockerImage, env, genome):
#	md5 = getMd5(run)
#	containerName = "service-"+serviceName+"-RTG-"+md5
#	if serviceName == "CQI":
#		cmd = "docker run -dti --rm --name="+containerName+" --env-file="+env+" -v "+run+":"+run+" -v "+genome+":"+genome+" "+serviceDockerImage+" /bin/bash -c '/home1/BAS/lamouchj/CQI/bin/rtg_vcf.sh --cqi='"
#		print(cmd)
#	subprocess.call(cmd, shell=True)
#	createRunningFile(run, serviceName)
#
#	return containerName


def parseJsonCQI(JSON, CQI_SAMPLE):
	with open (JSON) as f:
		jload = json.load(f)
		for item in jload['CQI']:
			#print(item['name'])
			if item['name'] == CQI_SAMPLE:
				#print(item['VCF'])
				print("[#INFO] CQI " + item['name'] + " detected !")
				CQI_VCF = item['VCF']
				print("CQI_VCF="+CQI_VCF)
				return CQI_VCF
			else:
				continue
		print("[ERROR] Wrong CQI provided! got "+"{"+CQI_SAMPLE+ "}")
		print("[#INFO] Looking into other run")

#def nonRegression():
	#List of VCF REF for --cqi arg and list of VCF STARK for --cqigz
	


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

def main(groupInputList, serviceName, services, serviceDockerImage, JSON, env, minDelay, days, genome, cqiFolder, serviceDockerImageRtg, cqiBin):
	while True:
		for groupInput in groupInputList:
			starkCompleteList = getStarkComplete(groupInput, days)
			#print("groupInput="+groupInput)
			#print("serviceName="+serviceName)
			for starkComplete in starkCompleteList:
				run = getRunName(starkComplete)
				if checkTriggersFromConfig(services, serviceName, run):
					startService(run, serviceName, serviceDockerImage, JSON, env, cqiFolder, groupInput, genome, serviceDockerImageRtg, cqiBin)
		# Launch loop if list "containerName" is not empty
			with open (osj(groupInput,"ContainerList"+serviceName+".txt"), "r") as containerList:
							for runName in containerList:
								removeContainer(getMd5(runName.strip()), runName, serviceName, groupInput)
			time.sleep(60.0*minDelay)

def test(serviceName):
	mess ="[#INFO] Listener correctly start"
	print("serviceName="+serviceName)
	return mess

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
	parser.add_argument("-s", "--servicesfile", type = str, default = "/home1/BAS/lamouchj/CQI/bin/services.json", help = "path to the services.json configfile", dest = 'services')
	parser.add_argument("-e", "--env", type = str, help = "ENV file which contain envrionement variables", dest = 'env', required = True)
	parser.add_argument("-j", "--json", type = str, help = "JSON file to associate the CQI tag and the vcf", dest = 'JSON', required = True)
	parser.add_argument("-D", "--minDelay", type = int, default = 5, help = "daemon refresh delay", dest = 'minDelay')
	parser.add_argument("-t", "--nbDaysBack", type = int, default = 30, help = "folder older than x days", dest = 'days')
	parser.add_argument("-g", "--genome", type = str, default = "/home1/TOOLS/genomes/hg19/hg19.fa", help = "path to genome folder", dest = 'genome')
	parser.add_argument("-c", "--cqiFolder", type = str, default = "/home1/BAS/lamouchj/input/CQI_REF", help = "path to CQI ref vcf", dest = 'cqiFolder')
	parser.add_argument("-r", "--rtg", type = str, help = "docker image of rtg tools", dest = 'serviceDockerImageRtg', required = True)
	parser.add_argument("-b", "--bin", type = str, help = "Path of the CQI bin and", dest = 'cqiBin', required = True)
	return parser.parse_args()

if __name__ == "__main__":
	#Path for BAS use
	args = myoptions()
	print(args)
	print(args.groupListeningInput)
	groupInputList = args.groupListeningInput.split(",")
	serviceTags = args.serviceTagsInput.split("!")
	print(args.serviceName)
	main(groupInputList, args.serviceName, args.services, args.serviceDockerImage, args.JSON, args.env, args.minDelay, args.days, args.genome, args.cqiFolder, args.serviceDockerImageRtg, args.cqiBin)