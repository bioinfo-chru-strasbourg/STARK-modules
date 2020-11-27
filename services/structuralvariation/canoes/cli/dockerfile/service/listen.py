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

def listen(app,input,output,repository,genome,day,cmd):
	list_path=[]
	for group in app:
		path=os.path.join(input, group)
		folder_to_add=systemcall("find "+path+"/*/*/STARKCopyComplete.txt -mtime -"+str(day)+" 2>/dev/null")
		if not folder_to_add == ['']:
			list_path += folder_to_add

	run_list = [ "/".join(x.split('/')[:-1]) for x in list_path ]
	for run in run_list:
		folder_output = run.replace(input, output, 1)
		folder_repository = run.replace(input, repository, 1)
		if systemcall("find "+folder_output+" -maxdepth 1 -type d -name CANOES 2>/dev/null") == [""]:
			print("[INFO] - "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+": Analysis launched: "+run)
			command = cmd+" run --run "+run+" --genome "+genome+" --output "+folder_output+"/CANOES"+" --repository "+folder_repository
			systemcall(command)
			print("[INFO] - "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+": Analysis ended: "+run)
			return True

def main():
	args = myoptions()
	print("[INFO] - "+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+": canoes-service is listening to: "+",".join([ "/data/input/"+x for x in args.group ]))
	while True:
		binary = os.path.dirname(os.path.realpath(__file__))+"/../bin/canoes"
		listen(args.group,"/data/input","/data/output","/data/repository",args.genome,args.days,binary)
		time.sleep(60.0*args.delay)

def myoptions():
	'''
	*arg parser*
	*return options*
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--nbMinDelay", type=int, default=10, help="daemon refresh delay", dest='delay')
	parser.add_argument("-t", "--nbDaysBack", type=int, default=30, help="folder older than x days", dest='days')
	parser.add_argument("-l", "--genome", type=str, help="path to genome.fa file", dest='genome')
	parser.add_argument("-g", "--group", action='append', help="STARK group", dest='group', required=True)
	return parser.parse_args()

if __name__ == '__main__':
	main()