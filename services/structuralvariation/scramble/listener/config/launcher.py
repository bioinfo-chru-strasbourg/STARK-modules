#! /usr/bin/env python
# -*- coding: utf-8 -*-


# launcher for SCRAMBLE


from __future__ import division
from __future__ import print_function

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


def launch(run, serviceName, containersFile, montage, image, launchCommand, configFile, repository):
	containerName = serviceName+os.path.basename(run)
	cmd = "docker run --rm --name="+containerName+" -v "+run+":"+run+" "+montage+" "+image+" "+launchCommand+" --config DATA_DIR="+os.path.basename(run)+"'"
	subprocess.call(cmd, shell = True)
	createContainerFile(containersFile, run, containerName)