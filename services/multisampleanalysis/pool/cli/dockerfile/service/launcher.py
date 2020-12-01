#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import os
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

def launch(run, serviceName, containersFile, montage, image, launchCommand, configFile):
	createRunningFile(run, serviceName)
	containerName = serviceName+"-NAME-"+os.path.basename(run)
	cmd = "docker run --rm --name="+containerName+" -v "+run+":"+run+" "+montage+" "+image+" "+launchCommand+" -i "+run
	# cmd = "docker run --name="+containerName+" -v "+run+":"+run+" "+montage+" "+image+" "+launchCommand+" -i "+run
	# cmd = "docker run --name="+containerName+" -v "+run+":"+run+" "+montage+" "+image+" "+launchCommand
	subprocess.call(cmd, shell = True)
	createContainerFile(containersFile, run, containerName)