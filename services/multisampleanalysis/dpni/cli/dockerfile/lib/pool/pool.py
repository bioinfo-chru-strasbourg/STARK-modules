#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function
import argparse
import json
import os
import random
import string
import subprocess
import sys

def systemcall(command):
	"""
	Passing command to the shell, return list containing stdout lines
	"""
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read().decode('utf8').strip().split('\n')

def randomStringDigits(stringLength=10):
	"""
	Generate a random string of letters and digits
	"""
	lettersAndDigits = string.ascii_letters + string.digits
	return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))

def basename(n):
	"""
	Extract basename from path
	input:  "/foo/bar/name.file.ext.txt"
	output: "name"
	"""
	return n.split('/')[-1].split('.')[0]

def current_dir():
	return os.path.dirname(os.path.realpath(__file__))

def add_option(config,option,name,group,option_type,required=False):
	"""
	Add options to the config file
	option_type f: file
	option_type d: directory
	option_type n: number
	option_type l: list
	"""
	if option:
		if option_type == "f" or option_type == "d":
			config[group][name] = os.path.abspath(option)
		elif option_type == "l":
			if name == "sample":
				config[group][name] = getSample(option)
			elif name == "pool":
				config[group][name] = getPool(option)
			else:
				config[group][name] = option
		else:
			config[group][name] = option
	if required and name not in config[group]:
		print("[ERROR] - please define option "+ name +"=<"+option_type+"> ",file=sys.stderr)
		sys.exit()
	return config

def getConfig(args):
	"""
	Get the snakemake configuration dict from options
	"""
	
	#Get config from file
	if args.config:
		config = readjson(args.config)
	elif os.path.exists(os.path.join(current_dir(), "../../config/config.json")):
		config = readjson(os.path.join(current_dir(), "../../config/config.json"))
	else:
		config = {}

	if "tools" not in config:
		config["tools"] = {}
	if "analysis" not in config:
		config["analysis"] = {}	

	#Get config from option
	#analysis
	config=add_option(config,args.output,"output","analysis","d",True)
	config=add_option(config,args.tmp,"tmp","analysis","d")
	config=add_option(config,args.sampleList,"sample","analysis","l",True)
	config=add_option(config,args.poolList,"pool","analysis","l",True)
	config=add_option(config,args.bed,"bed","analysis","f",True)
	config=add_option(config,args.genome,"genome","analysis","f")
	config=add_option(config,args.padding,"padding","analysis","n")
	config=add_option(config,args.config,"config","analysis","f")
	#tools
	config=add_option(config,args.java,"java","tools","f",True)
	config=add_option(config,args.gatk,"gatk","tools","f",True)
	config=add_option(config,args.howard,"howard","tools","f",True)
	config=add_option(config,args.bcftools,"bcftools","tools","f",True)
	config=add_option(config,args.bgzip,"bgzip","tools","f",True)
	config=add_option(config,args.tabix,"tabix","tools","f",True)
	return config

def getSample(sampleList):
	"""
	Extract sample from the options
	"""
	sample_dict = {}
	sample_list = [ x.strip() for x in sampleList.strip().split(",") ]
	for sample in sample_list:
		sample_name = basename(sample)
		if sample_name not in sample_dict:
			sample_dict[sample_name]={}
			sample_dict[sample_name]["vcf"]=sample
	return sample_dict

def getPool(poolList):
	"""
	Extract pool from the options
	"""
	pool_dict = {}
	pool_list = [ x.strip().split(":") for x in poolList.strip().split(",") ]
	for pool in pool_list:
		pool_name = pool[0].strip()
		pool_vcf  = pool[1].strip()
		pool_vcf_name  = pool_vcf.split("/")[-1].split(".")[0]
		pool_bam  = pool[2].strip()
		if pool_name not in pool_dict:
			pool_dict[pool_name]={}
			pool_dict[pool_name]["vcf"]=pool_vcf
			pool_dict[pool_name]["bam"]=pool_bam
			pool_dict[pool_name]["name"]=pool_vcf_name
	return pool_dict

def nbCPUs():
	"""
	Get the nuber of CPUs
	"""
	command = "lscpu | grep -E '^CPU\\(s\\)' | awk '{print $NF}'"
	return systemcall(command)

def readjson(f):
	"""
	Read a json file, return dict
	"""
	with open(f) as infile:
		datastore = json.load(infile)
	return datastore

def writejson(d,f):
	"""
	Read dict, return a json file
	"""
	with open(f, 'w') as outfile:
		json.dump(d, outfile, indent=4)

def pool(config,output,jobs,dryrun,dag,extraconf):
	"""
	Launch snakefile for pool
	"""
	if not jobs:
		#jobs = str(int(nbCPUs()[0]) -1)
		jobs = str(4)
	snk_opt = ""
	if dryrun:
		snk_opt = " -np "
	if dag:
		snk_opt += " --dag | dot -Tsvg > " + output + "/dag.svg "
	current_dir = os.path.dirname(os.path.realpath(__file__))

	#extra options
	if extraconf:
		extra = " --config "+extraconf
	else:
		extra = ""

	command = "/usr/local/lib/miniconda3/bin/snakemake -p --snakefile " + current_dir + "/../snakemake/Snakefile --configfile " + config + extra+" --directory " + output + " --jobs "+ jobs + snk_opt
	print(command)
	systemcall(command)

def main_sample(args):
	"""
	Handle with sample suboptions
	"""
	config = getConfig(args)
	configfile = config["analysis"]["output"]+"/analysis.json"
	writejson(config,configfile)
	pool(configfile,config["analysis"]["output"],args.jobs,args.dryrun,args.dag, args.extraconf)

def main_run(args):
	"""
	Handle with run suboptions
	"""
	# to do
	pass

def main():
	"""
	Main function
	"""
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(help='sub-command help')

	#Sample suboptions
	parser_a = subparsers.add_parser("sample", help="add new sample to VaRank folder")
	
	group0 = parser_a.add_argument_group('Analysis options')
	group0.add_argument("-o","--output",type=str,help="The output folder",dest="output",required=True)
	group0.add_argument("-t","--tmp",type=str,help="The temporary folder",dest="tmp")
	group0.add_argument("-s","--sampleVcf", help="The sample vcf files", type=str, dest='sampleList',required=True)
	group0.add_argument("-p","--pool", help="pool vcf and bam files: poolName:poolVcf:poolBam, example: POOL_F,A.vcf,A.bam", type=str, dest='poolList',required=True)
	group0.add_argument("-b","--bed", help="analysis bed file", type=str, dest='bed')
	group0.add_argument("-g","--genome", help="analysis genome file", type=str, dest='genome')
	group0.add_argument("-n","--padding", help="The padding used for the variant calling", type=str, dest='padding')
	group0.add_argument("-c","--config", help="The json config file", type=str, dest='config')
	group0.add_argument("-x", "--extraconf", type=str, help="extraconf to pass to pool script", dest='extraconf')
	
	group1 = parser_a.add_argument_group('Ressources')
	group1.add_argument("--java",type=str,help="The path to the java executable")
	group1.add_argument("--gatk",type=str,help="The path to the GATK executable")
	group1.add_argument("--howard",type=str,help="The path to the HOWARD folder")
	group1.add_argument("--bcftools",type=str,help="The path to the bcftools executable")
	group1.add_argument("--bgzip",type=str,help="The path to the bgzip executable")
	group1.add_argument("--tabix",type=str,help="The path to the tabix folder")

	group2 = parser_a.add_argument_group('Snakemake options')
	group2.add_argument("-z","--jobs",type=int,help="number of parallel cores to use (default: nbCPUs -1)",dest="jobs")
	group2.add_argument("-d","--dag",action="store_true",help="write dag",dest="dag")
	group2.add_argument("-y","--dryrun",action="store_true",help="dry run",dest="dryrun")
	
	parser_a.set_defaults(mode=main_sample)
	
	#Run suboptions
	#to do

	args = parser.parse_args()
	args.mode(args)

if __name__ == '__main__':
	main()
