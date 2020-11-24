#!/usr/bin/python
# -*- coding: utf-8 -*-

import time
import json
import argparse
import os, sys, subprocess
from datetime import datetime


def systemcall(command,log=False):
	"""
	passing command to the shell
	return list containing stdout lines
	"""
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	if log:
		with open(log, "a") as f:
			f.write("\n>>>[START] - "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+"\n")
			f.write(p.stdout.read().decode('utf8'))
			f.write("<<<[END] - "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+"\n")
	return p.stdout.read().decode('utf8').strip().split('\n')

"""
def launch():
    
    Once container is launched, infinite loop which identify Starkcomplete files lower than 30 days in the repository.
    return list containing 
    
    while True:
"""

def listen(app,input,output,repository,day,cmd):
    list_path=[]
	for group in app:
		path=os.path.join(input, group)
        # -mtime -n   means less than args.day
		folder_to_add=systemcall("find "+path+"/*/*/STARKCopyComplete.txt -mtime -"+str(day)+" 2>/dev/null")
		if not folder_to_add == ['']:
			list_path += folder_to_add

	run_list = [ "/".join(x.split('/')[:-1]) for x in list_path ]
	for run in run_list:
		folder_output = run.replace(input, output, 1)
		folder_repository = run.replace(input, repository, 1)
		if systemcall("find "+folder_output+" -maxdepth 1 -type d -name CQI 2>/dev/null") == [""]:
			print("[INFO] - "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+": Analysis launched: "+run)
			command = cmd+" run --run "+run+" --genome "+genome+" --output "+folder_output+"/CQI"+" --repository "+folder_repository
			systemcall(command)
			print("[INFO] - "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+": Analysis ended: "+run)
			return True

def get_sample_list_from_samplesheet(samplesheetPath):
        """
        Returns a python list containing all sample names in a samplesheet.
        """
        assert samplesheetPath != "NO_SAMPLESHEET_FOUND", \
                        "[ERROR] find_any_samplesheet() couldn't find any samplesheet. Check if the --fromResultDir argument is set correctly."
        assert_file_exists_and_is_readable(samplesheetPath)
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


def find_any_samplesheet(runDir, fromResDir):
        """
        1) look up recursively all files named SampleSheet.csv in the runDir
        2) check if file path follows an expected samplesheet name and location
                (the latter depends on if we're in a STARK result or repository dir,
                defined by the bool fromResDir)
        3) first correct file path is returned
        """
        p = subprocess.Popen("find -L "+runDir+" -name *SampleSheet.csv", stdout=subprocess.PIPE, shell=True)
        out = p.stdout.readlines()
        for ss in out:
                ss = ss.decode("utf-8").strip()
                if fromResDir:
                        r = re.match(runDir.rstrip("/")+"/(.*)/(.*).SampleSheet.csv", ss)
                else:
                        r = re.match(runDir.rstrip("/")+"/(.*)/DATA/(.*).SampleSheet.csv", ss)
                if r is None:
                        continue
                elif r.group(1) == r.group(2): #checks if (.*) == (.*)
                        return ss
        return "NO_SAMPLESHEET_FOUND"


def get_tags_list(sampleList, sampleDirList, samplesheet=""):
        """
        Returns a list of tags string in the order corresponding to sampleList

        Tries to fetch tags through STARK analysis and sample tags files.
        If it fails and a samplesheet is provided, fetches them from there instead.
        """
        try:
                tagsList = []
                for sample, sampleDir in zip(sampleList, sampleDirList):
                        tags = ""
                        with open(osj(sampleDir, sample+".analysis.tag"), "r") as f:
                                for l in f:
                                        tags += l.strip()
                        with open(osj(sampleDir, sample+".tag"), "r") as f:
                                for l in f:
                                        if tags != "":
                                                tags += "!"+l.strip()
                                        else:
                                                tags = l.strip()
                        tagsList.append(tags)
                print(tagsList)
                return tagsList
                # return ["APP#SWAG!#TUMCELL" for v in sampleList]
        except IOError as e:
                print(e)
                print("[WARNING] Couldn't use <sample>.tag and <sample>.analysis.tag to determine tags, trying to use samplesheet instead.")
                assert samplesheet != "", "[ERROR] Tried to use samplesheet to determine tags, but no samplesheet was provided."
                if os.path.getsize(samplesheet) == 0:
                        #then no tags to be found
                        return ["" for v in sampleList]
                        # return ["APP#SWAG!#TUMCELL" for v in sampleList]
                else:
                        #get both samples and tags to make sure the tags returned
                        #are in the same order as the sampleList provided in input
                        samples = get_sample_list_from_samplesheet(samplesheet)
                        tags = get_descriptions_from_samplesheet(samplesheet)
                        return [tags[samples.index(s)] for s in sampleList]

def main():
    args = myoptions()
    print("[INFO] - "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+": canoes-service is listening to: "+",".join([ "/data/input/"+x for x in args.group ]))
    binary = os.path.dirname(os.path.realpath(__file__))+"/../bin/canoes"
    while True:
        time.sleep(args.minutes)
        listen(args.group,"/data/input","/data/output","/data/repository",args.genome,args.days,binary)

def myoptions():
	'''
	*arg parser*
	*return options*
	'''
	parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--lenghtToSleep" , type=int, default=5, help="time between each starkcomplete search", dest='minutes') 
	parser.add_argument("-t", "--nbDaysBack", type=int, default=30, help="folder older than x days", dest='days')
	parser.add_argument("-l", "--genome", type=str, help="path to genome.fa file", dest='genome')
	parser.add_argument("-g", "--group", action='append', help="STARK group", dest='group', required=True)
	return parser.parse_args()

if __name__ == '__main__':
	main()