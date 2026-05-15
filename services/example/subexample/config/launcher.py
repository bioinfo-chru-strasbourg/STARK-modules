import os
import subprocess

from datetime import datetime
from os.path import join as osj

def write_services_log(stark_services_dir: str, run: str, container_name: str):
	with open(osj(stark_services_dir,container_name+".log"), "w+") as f:
		f.write("RUN: "+os.path.basename(run)+"\n")
		f.write("FOLDER: "+run+"\n")
		f.write("EXEC_DATE: "+datetime.now().strftime("%d%m%Y-%H%M%S")+"\n")
		f.write("ID: "+container_name+"\n")

def create_running_file(run: str, service_name: str):
	with open(osj(run,service_name+"Running.txt"), "w+") as f:
		f.write("# ["+datetime.now().strftime("%d/%m/%Y %H:%M:%S")+"] "+os.path.basename(run)+" running with "+service_name+"\n")

def launch(run: str, service_name: str, stark_services_dir: str, mounts: str, image: str, launch_command: str, config_file: str, repository: str):
	"""
	This function's args must follow the template defined in STARK's common listener service
	"""
	create_running_file(run, service_name)
	container_name = service_name+"-NAME-"+os.path.basename(run)
	cmd = "docker run --rm --name="+container_name+" --volumes-from stark-module-example-submodule-subexample-service-cli "+image+" "+launch_command+" --fibonacci_n 10 -i "+run
	print(cmd)
	subprocess.call(cmd, shell = True)
	write_services_log(stark_services_dir, run, container_name)