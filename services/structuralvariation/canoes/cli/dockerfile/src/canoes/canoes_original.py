#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os, subprocess, sys
import argparse
import json
import glob
import re
import shutil
import time
from datetime import datetime, timedelta
from os.path import join as osj

class Run():
	def __init__(self, path, samplesheet, bed):
		"""
		initialize the Run object
		"""
		self.path    = path
		self.info    = self.extract(path, samplesheet, bed)
		self.sample  = self.info.keys()
		self.bam     = [ self.info[x]["bam"] for x in self.sample ]
		self.sex     = [ self.info[x]["sex"] for x in self.sample ]
		self.bed     = [ self.info[x]["bed"] for x in self.sample ][0]

	def extract(self, path, samplesheet, bed):
		"""
		Extract run parameters from run folder
		Based on the folder configuration of EMCLABO (11/09/2019)
		"""
		info = {}
		ana_samplesheet = False
		command = 'find '+self.path+' -maxdepth 1 -name "*SampleSheet.csv"'
		sample_list = get_sample_id_from_samplesheet(systemcall(command)[0])
		if not sample_list:
			command = 'find '+path+' -maxdepth 3 -name "*SampleSheet.csv"'
			sample_list = get_sample_id_from_samplesheet(systemcall(command)[0])
		for sample in sample_list:
			if sample not in info:
				info[sample]={}
			command = 'find -L '+path+' -maxdepth 3 -name "'+sample+'*.bam" | paste -sd, - | awk -F "," \'{for(i=1;i<=NF;i++){if($i !~ /validation/){print $i}}}\''
			info[sample]["bam"] = systemcall(command)[0]
			if not bed:
				command = 'find -L '+path+' -maxdepth 3 -name "*bed" | paste -sd, - | awk -F "," \'{for(i=1;i<=NF;i++){if($i !~ /genes/){print $i}}}\''
				info[sample]["bed"] = systemcall(command)[0]
			else:
				info[sample]["bed"] = bed
			if not samplesheet:
				command = 'find '+path+' -maxdepth 1 -name "*SampleSheet.csv"'
				info[sample]["SampleSheet"] = systemcall(command)[0]
				if not info[sample]["SampleSheet"]:
					command = 'find '+osj(path,sample)+' -maxdepth 4 -name "*SampleSheet.csv"'
					info[sample]["SampleSheet"] = systemcall(command)[0]
			else:
				info[sample]["SampleSheet"] = samplesheet
			if not ana_samplesheet:
				ana_samplesheet = info[sample]["SampleSheet"]
		
		sex = {}
		nextline = True
		description = -1
		if not ana_samplesheet:
			print("[ERROR] - unable to find a SampleSheet in this run: "+path+" (EXIT)",file=sys.stderr)
			sys.exit()
		with open(ana_samplesheet, "r") as f:
			for line in f:
				if nextline:
					if line.startswith("Sample_ID"):
						nextline = False
						description = line.strip().split(",").index("Description")
				if not nextline and line.strip():
					field = line.strip().split(",")
					tags = field[description].split("!")
					ID = field[0]
					for tag in tags:
						if "SEX#" in tag:
							sex[ID] = tag.split("#")[-1]
						elif "SEX_" in tag:
							sex[ID] = tag.split("_")[-1]
						else:
							sex[ID] = ""
		for sample in sample_list:
			if sample in sex:
				info[sample]["sex"]=sex[sample]
			else:
				print("[WARNING] - Sample ID:\""+sample+"\" removed from analysis (sample was not find in SampleSheet)",file=sys.stderr)
				print("[INFO] - Could it be that you try to use a wrong SampleSheet?")
				del info[sample]
		return info

def write_json(output,bed,bam,sample,sex,ana,exclude,gatk,bedtools,annotsv,java,R,genome,window,version):
	config_file = os.path.dirname(os.path.realpath(__file__))+"/../../config/config.json"
	with open(config_file, 'r') as f:
		datastore = json.load(f)

	data = {}
	data['tools'] = {}
	data['tools']["canoes"] = os.path.dirname(os.path.realpath(__file__))
	data['tools']["canoes_version"] = str(version)
	if gatk:
		data['tools']["gatk"] = os.path.abspath(gatk)
	elif "gatk" in datastore:
		data['tools']["gatk"] = os.path.abspath(datastore["gatk"])
	else:
		print("[ERROR] - please define path to gatk binary: -k --gatk <PATH> (EXIT)",file=sys.stderr)
		sys.exit()
	if bedtools:
		data['tools']["bedtools"] = os.path.abspath(bedtools)
	elif "bedtools" in datastore:
		data['tools']["bedtools"] = os.path.abspath(datastore["bedtools"])
	else:
		print("[ERROR] - please define path to gatk binary: -t --bedtools <PATH> (EXIT)",file=sys.stderr)
		sys.exit()
	if annotsv:
		data['tools']["annotsv"] = os.path.abspath(annotsv)
	elif "annotsv" in datastore:
		data['tools']["annotsv"] = os.path.abspath(datastore["annotsv"])
	else:
		print("[ERROR] - please define path to gatk binary: -a --annotsv <PATH> (EXIT)",file=sys.stderr)
		sys.exit()
	if java:
		data['tools']["java"] = os.path.abspath(java)
	elif "java" in datastore:
		data['tools']["java"] = os.path.abspath(datastore["java"])
	else:
		print("[ERROR] - please define path to gatk binary: -j --java <PATH> (EXIT)",file=sys.stderr)
		sys.exit()
	if R:
		data['tools']["R"] = os.path.abspath(R)
	elif "R" in datastore:
		data['tools']["R"] = os.path.abspath(datastore["R"])
	else:
		print("[ERROR] - please define path to gatk binary: -r --R <PATH> (EXIT)",file=sys.stderr)
		sys.exit()
	
	data['analysis'] = {}
	data['analysis']["type"] = cleanList(ana)
	if genome:
		data['analysis']["genome"] = os.path.abspath(genome)
	elif "genome" in datastore:
		data['analysis']["genome"] = os.path.abspath(datastore["genome"])
	else:
		print("[ERROR] - please defie path to genome.fa: -g --genome <PATH> (EXIT)",file=sys.stderr)
		sys.exit()
	data['analysis']["output"] = osj("/app",os.path.basename(os.path.normpath(output)))
	data['analysis']["repository"] = os.path.abspath(output)
	data['analysis']["bed"] = os.path.abspath(bed)
	if window:
		data['analysis']["window"] = window

	data['analysis']['sample'] = []
	bamlist = cleanList(bam)
	for i in range(0,len(bamlist)):
		sample_data = {}
		
		sample_data['bam'] = os.path.abspath(bamlist[i])
		if sample:
			sample_list = sample.split(',')
			sample_ctrl = sample_list[i]
		else:
			sample_ctrl = bamlist[i].split('/')[-1].split('.')[0]

		if exclude:
			if not isValidSample(sample_ctrl,cleanList(exclude)):
				continue

		sample_data['name'] = sample_ctrl
		if sex:
			sex_list = cleanList(sex)
			if sex_list[i]:
				sample_data['sex'] = sex_list[i]
			else:
				sample_data['sex'] = 'UNKNOWN'
		else:
			sample_data['sex'] = 'UNKNOWN'
		data['analysis']['sample'].append(sample_data)

	with open(output+'/analysis.json', 'w') as outfile:
		json.dump(data, outfile, indent=4)

def systemcall(command):
	'''
	*passing command to the shell*
	*return list containing stdout lines*
	command - shell command (string)
	'''
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read().decode('utf8').strip().split('\n')

def get_sample_id_from_samplesheet(samplesheetPath):
	"""
	Adapted from functions.py
	
	Returns a python list containing all sample ID from samplesheet.
	"""
	if samplesheetPath == "NO_SAMPLESHEET_FOUND" or not samplesheetPath:
		return []
	inDataTable = False
	sampleID = []
	with open(samplesheetPath, "r") as f:
		for l in f:
			if not inDataTable:
				if l.startswith("Sample_ID,"):
					inDataTable = True
					sampleIDIndex = l.strip().split(",").index("Sample_ID")
			else:
				if "," in l:
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
		if ss:
			return ss
	return "NO_SAMPLESHEET_FOUND"

def nbCPUs():
	command = "lscpu | grep -E '^CPU\\(s\\)' | awk '{print $NF}'"
	return systemcall(command)

def isValidSample(sample,exclude):
	for item in exclude:
		if sample.startswith(item) or sample.endswith(item):
			return False
	return True

def createDir(path):
	if not os.path.isdir(path):
		os.makedirs(path)
		os.chmod(path, 0o777)
		print("Directory ",path," Created")

def cleanList(str):
	return [ x.strip() for x in str.split(",") ]

def cleanString(lst):
	return ",".join([ x.strip() for x in lst ])

def canoes(config,output,jobs,dryrun,dag):
	if not jobs:
		jobs = str(int(nbCPUs()[0]) -1)
	snk_opt = ""
	if dryrun:
		snk_opt = " -np "
	if dag:
		snk_opt += " --dag | dot -Tsvg > " + output + "/dag.svg "
	current_dir = os.path.dirname(os.path.realpath(__file__))
	command = "/usr/local/lib/miniconda3/bin/snakemake --snakefile " + current_dir + "/../snakemake/Snakefile --configfile " + config + " --directory " + output + " --jobs "+ jobs + snk_opt
	systemcall(command)
	coverageTsvFile = osj(output,"CANOES/all.canoes.coverage.tsv")
	if os.path.exists(coverageTsvFile):
		generateGraphs(coverageTsvFile, output)

def generateGraphs(coverageTsvFile, output):
	createStatsFiles(coverageTsvFile)
	writePlotScript(coverageTsvFile, output)
	command = "R CMD BATCH "+osj(output,"CANOES/script_plot.R")+" "+osj(output,"logs/R_plotting.log")
	systemcall(command)

def createStatsFiles(coverageTsvFile):
	command = "grep -vE \"chrX|chrY\" "+coverageTsvFile+" > "+os.path.dirname(coverageTsvFile)+"/all.CANOES.stats.tsv.tmp;\\\
head -1 "+os.path.dirname(coverageTsvFile)+"/all.CANOES.stats.tsv.tmp | awk '{printf $1\":\"$2\"-\"$3\"\\t\"} {for(i=4;i<=NF;i++){printf \"%s\\t\", $i}; printf \"\\n\"}' > "+os.path.dirname(coverageTsvFile)+"/all.canoes.stats.tsv;\\\
cat "+os.path.dirname(coverageTsvFile)+"/all.CANOES.stats.tsv.tmp | sed -e \"1d\"| awk '{printf $1\":\"$2\"-\"$3\"\\t\"} {for(i=4;i<=NF;i++){printf \"%s\\t\", $i/($3-$2)}; printf \"\\n\"}' | sed 's/\,/./g' >> "+os.path.dirname(coverageTsvFile)+"/all.canoes.stats.tsv;\\\
sed -i 's/.$//' "+os.path.dirname(coverageTsvFile)+"/all.canoes.stats.tsv;\\\
head -1 "+os.path.dirname(coverageTsvFile)+"/all.CANOES.stats.tsv.tmp | awk '{printf $1\":\"$2\"-\"$3\"\\t\"} {for(i=4;i<=NF;i++){printf \"%s\\t\", $i}; printf \"\\n\"}' > "+os.path.dirname(coverageTsvFile)+"/all.canoes.stats2.tsv;\\\
cat "+os.path.dirname(coverageTsvFile)+"/all.CANOES.stats.tsv.tmp | sed -e \"1d\"| awk '{printf $1\":\"$2\"-\"$3\"\\t\"} {for(i=4;i<=NF;i++){printf \"%s\\t\", $i/($3-$2)}; printf \"\\n\"}' | sed 's/\,/./g' >> "+os.path.dirname(coverageTsvFile)+"/all.canoes.stats2.tsv;\\\
sed -i 's/.$//' "+os.path.dirname(coverageTsvFile)+"/all.canoes.stats2.tsv;"
	systemcall(command)

def writePlotScript(coverageTsvFile, output):
	rScript = osj(os.path.dirname(os.path.realpath(__file__)),"plotCoverage.R")
	with open(osj(os.path.dirname(coverageTsvFile),"all.canoes.A.coverage.tsv.conv"), "r") as f:
		sample_list = ",".join([ "\""+x+"\"" for x in f.read().split('\n')[0].split("\t") ][3:] )
	with open(osj(output,"CANOES/script_plot.R"), "w") as f:
		f.write("setwd(\""+ output +"\")\n")
		f.write("source(\""+ rScript +"\")\n")
		f.write("sample.names <- c("+ sample_list +")\n")
		f.write("for (i in 1:length(sample.names)){\n")
		f.write('PlotCoverage(paste(sample.names[i]),\"'+os.path.dirname(coverageTsvFile)+'/all.canoes.stats.tsv\",paste(\"'+output+'/\",sample.names[i],\"/CANOES/\",sample.names[i],\".canoes.boxplot.coverage.pdf\", sep=\"\"))\n')
		f.write('PlotCoverage2(paste(sample.names[i]),\"'+os.path.dirname(coverageTsvFile)+'/all.canoes.stats2.tsv\",paste(\"'+output+'/\",sample.names[i],\"/CANOES/\",sample.names[i],\".canoes.barplot.coverage.pdf\", sep=\"\"))\n')
		f.write("}\n")
		f.close()

def copyAllResults(output, exclude, depository, samplesheet):
	sample_list = []
	if not samplesheet:
		sample_list = get_sample_id_from_samplesheet(find_any_samplesheet(output))
	else:
		sample_list = get_sample_id_from_samplesheet(samplesheet)
	results = osj("/app",os.path.basename(os.path.normpath(output)))
	if os.path.exists(osj(output,".snakemake/")):
		shutil.rmtree(osj(output,".snakemake/"))
	if not sample_list:
		command = 'find '+osj(results,"*/CANOES/*.canoes.annotsv.tsv")
		tsv_list = systemcall(command)
		for sample in tsv_list:
			sample_list.append(os.path.basename(sample)[:-19])
	for sample in sample_list:
		if isValidSample(sample, cleanList(exclude)):
			if not os.path.isdir(osj(output,sample,"CANOES/")):
				createDir(osj(output,sample,"CANOES/"))
			if not os.path.exists(osj(depository,sample,"CANOES/")):
				createDir(osj(depository,sample,"CANOES/"))
			shutil.copy2(osj(results,"CANOES/all.canoes.annotsv.sorted.tsv"),osj(output,sample,"CANOES/"))
			shutil.copy2(osj(results,"CANOES/all.canoes.annotsv.sorted.tsv"),osj(depository,sample,"CANOES/"))
			shutil.copy2(osj(results,sample,"CANOES",sample+".canoes.annotsv.tsv"),osj(output,sample,"CANOES/"))
			shutil.copy2(osj(results,sample,"CANOES",sample+".canoes.annotsv.tsv"),osj(depository,sample,"CANOES/"))
			if not os.path.exists(osj(output,sample,"CANOES/logs")):
				shutil.copytree(osj(results,"logs/"),osj(output,sample,"CANOES/logs"))
			if not os.path.exists(osj(depository,sample,"CANOES/logs")):
				shutil.copytree(osj(results,"logs/"),osj(depository,sample,"CANOES/logs"))

def runDepositoryGenerator(output):
	if (output).endswith('/'):
		runRepository = output[-1]
	else:
		runRepository = output
	runDepository = re.sub(os.environ['DOCKER_MODULE_STRUCTURALVARIATION_SUBMODULE_CANOES_FOLDER_REPOSITORY'], os.environ['DOCKER_MODULE_STRUCTURALVARIATION_SUBMODULE_CANOES_FOLDER_DEPOSITORY'], runRepository)
	return runDepository

def createCompleteFile(output):
	file = open(osj(output,"CANOESComplete.txt"), "w+")
	file.write("# ["+(datetime.now() + timedelta(hours = 2)).strftime("%d/%m/%Y %H:%M:%S")+"] "+os.path.basename(output)+" analyzed with CANOES\n")
	file.close()
	if os.path.exists(osj(output,"CANOESRunning.txt")):
		os.remove(osj(output,"CANOESRunning.txt"))

def main_sample(args):
	createDir(args.output)
	depository = runDepositoryGenerator(args.output)
	configfile = args.output+'/analysis.json'
	write_json(args.output,args.bed,args.bam,args.sample,args.sex,args.ana,args.exclude,args.gatk,args.bedtools,args.annotsv,args.java,args.R,args.genome,args.window,args.version)
	canoes(configfile,args.output,args.jobs,args.dryrun,args.dag)
	copyAllResults(args.output, args.exclude, depository)
	createCompleteFile(args.output)
	if args.rep:
		systemcall("cp -r "+args.output+"/* "+args.rep)

def main_run(args):
	createDir(args.output)
	depository = runDepositoryGenerator(args.output)
	configfile = args.output+'/analysis.json'
	run    = Run(args.run,args.samplesheet,args.bed)
	bed    = run.bed
	bam    = cleanString(run.bam)
	sample = cleanString(run.sample)
	sex    = cleanString(run.sex)
	write_json(args.output,bed,bam,sample,sex,args.ana,args.exclude,args.gatk,args.bedtools,args.annotsv,args.java,args.R,args.genome,args.window,args.version)
	canoes(configfile,args.output,args.jobs,args.dryrun,args.dag)
	copyAllResults(args.output, args.exclude, depository, args.samplesheet)
	createCompleteFile(args.output)
	if args.rep:
		createDir(args.rep)
		systemcall("cp -r "+args.output+"/* "+args.rep)

def main():
	'''
	*arg parser*
	*return options*
	'''
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(help='sub-command help')

	parser_a = subparsers.add_parser("sample", help="CANOES sample analysis")
	group0 = parser_a.add_argument_group('Analysis options')
	group0.add_argument("-b","--bam",type=str,help="bam files - \"file1,file2,...\"",dest="bam",required=True)
	group0.add_argument("-s","--sample",type=str,help="sample names \"file1,file2,...\"",dest="sample")
	group0.add_argument("-x","--sex",type=str,help="sample sex \"sex1,sex2,...\"",dest="sex")
	group0.add_argument("-l","--bed",type=str,help="the bed file",dest="bed",required=True)
	group0.add_argument("-e","--ana",type=str,help="list of analysis to run (A:autosomal, F:sex F, M:sex M) - Default: \"A,F,M\"",dest="ana",default="A,M,F")
	group0.add_argument("-c","--exclude",type=str,help="list of patern to exclude - Default: \"POOL_,BlcADN,blanc,BlcPCR,blcPCR,Z_NTC\"",dest="exclude",default="POOL_,BlcADN,blanc,BlcPCR,blcPCR,Z_NTC")
	group0.add_argument("-w","--window",type=int,help="cut bed region in windows of n bp",dest="window")
	group0.add_argument("-v","--version",type=int,help="CANOES version - 1(default) or 2 available",default=1,dest="version")
	group0.add_argument("-o","--output",type=str,help="output folder",dest="output",required=True)
	group0.add_argument("-u","--repository", help="Repository folder", type=str, dest='rep')

	group1 = parser_a.add_argument_group('Ressources')
	group1.add_argument("-g","--genome",type=str,help="path to the genome.fa file",dest="genome",)
	group1.add_argument("-i","--R",type=str,help="path to the R executable",dest="R")
	group1.add_argument("-j","--java",type=str,help="path to the java executable",dest="java")
	group1.add_argument("-t","--bedtools",type=str,help="path to the bedtools executable",dest="bedtools")
	group1.add_argument("-k","--gatk",type=str,help="path to the gatk executable",dest="gatk")
	group1.add_argument("-a","--annotsv",type=str,help="path to the annotsv executable",dest="annotsv")

	group2 = parser_a.add_argument_group('Snakemake options')
	group2.add_argument("-z","--jobs",type=int,help="number of parallel cores to use (default: nbCPUs -1)",dest="jobs")
	group2.add_argument("-d","--dag",action="store_true",help="write dag",dest="dag")
	group2.add_argument("-y","--dryrun",action="store_true",help="dry run",dest="dryrun")
	parser_a.set_defaults(mode=main_sample)

	parser_b = subparsers.add_parser("run", help="CANOES run analysis")
	group0 = parser_b.add_argument_group('Analysis options')
	group0.add_argument("-r","--run", help="The run folder", type=str, required=True, dest="run")
	group0.add_argument("-s","--samplesheet", help="path to SampleSheet.csv", type=str, dest='samplesheet')
	group0.add_argument("-l","--bed",type=str,help="the bed file",dest="bed")
	group0.add_argument("-e","--ana",type=str,help="list of analysis to run (A:autosomal, F:sex F, M:sex M) - Default: \"A,F,M\"",dest="ana",default="A,M,F")
	group0.add_argument("-c","--exclude",type=str,help="list of patern to exclude - Default: \"POOL_,BlcADN,blanc,BlcPCR,blcPCR,Z_NTC\"",dest="exclude",default="POOL_,BlcADN,blanc,BlcPCR,blcPCR,Z_NTC")
	group0.add_argument("-w","--window",type=int,help="cut bed region in windows of n bp",dest="window")
	group0.add_argument("-v","--version",type=int,help="CANOES version - 1(default) or 2 available",default=1,dest="version")
	group0.add_argument("-o","--output",type=str,help="output folder",dest="output",required=True)
	group0.add_argument("-u", "--repository", help="Repository folder", type=str, dest='rep')

	group1 = parser_b.add_argument_group('Ressources')
	group1.add_argument("-g","--genome", help="Path to the genome.fa file", type=str, dest='genome')
	group1.add_argument("-i","--R",type=str,help="path to the R executable",dest="R")
	group1.add_argument("-j","--java",type=str,help="path to the java executable",dest="java")
	group1.add_argument("-t","--bedtools",type=str,help="path to the bedtools executable",dest="bedtools")
	group1.add_argument("-k","--gatk",type=str,help="path to the gatk executable",dest="gatk")
	group1.add_argument("-a","--annotsv",type=str,help="path to the annotsv executable",dest="annotsv")

	group2 = parser_b.add_argument_group('Snakemake options')
	group2.add_argument("-z","--jobs",type=int,help="number of parallel cores to use (default: nbCPUs -1)",dest="jobs")
	group2.add_argument("-d","--dag",action="store_true",help="write dag",dest="dag")
	group2.add_argument("-y","--dryrun",action="store_true",help="dry run",dest="dryrun")
	parser_b.set_defaults(mode=main_run)
	
	args = parser.parse_args()
	args.mode(args)

if __name__ == '__main__':
	main()
