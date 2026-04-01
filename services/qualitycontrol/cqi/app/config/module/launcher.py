#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import hashlib
import os
import re
import subprocess
import json
import argparse
from datetime import datetime
from os.path import join as osj

date_time = datetime.now().strftime("%Y%m%d-%H%M%S")


def readconfig(configFile, serviceName, configkey):
    """Reads a specific key configuration from a JSON file."""
    with open(configFile, "r") as f:
        json_config = json.load(f)
    try:
        outputconfig = json_config["services"][serviceName][configkey]
        return outputconfig
    except KeyError:
        raise ValueError(f"[ERROR] Missing {configkey} for service {serviceName} in the config file.")


def createContainerFile(containersFile, run, containerName):
    """Creates a log file for the container execution."""
    with open(osj(containersFile, containerName + ".log"), "w+") as file:
        file.write("RUN: " + os.path.basename(run) + "\n")
        file.write("FOLDER: " + run + "\n")
        file.write("EXEC_DATE: " + datetime.now().strftime("%d%m%Y-%H%M%S") + "\n")
        file.write("ID: " + containerName + "\n")


def getMd5(run):
    """Generates an MD5 hash for the run folder."""
    runMd5 = hashlib.md5()
    runMd5.update(run.encode('utf-8'))
    return runMd5.hexdigest()


def findAnyBed(run):
    """Finds the first .bed file within the run directory."""
    p = subprocess.Popen(f"find -L {run} -maxdepth 3 -name '*.bed'", stdout=subprocess.PIPE, shell=True)
    out = p.stdout.readlines()
    for bed in out:
        bed = bed.decode("utf-8").strip()
        r = re.match(f"{run.rstrip('/')}/(.*)/STARK/(.*).bed", bed)
        if r is None:
            continue
        elif r.group(1) == r.group(2):  # checks if (.*) == (.*)
            return bed
    return "NO_BED_FOUND"


def find_any_samplesheet(runDir, fromResDir=False):
    """Looks for SampleSheet.csv files in the run directory."""
    p = subprocess.Popen(f"find -L {runDir} -maxdepth 3 -name '*SampleSheet.csv'", stdout=subprocess.PIPE, shell=True)
    out = p.stdout.readlines()
    for ss in out:
        ss = ss.decode("utf-8").strip()
        if fromResDir:
            r = re.match(f"{runDir.rstrip('/')}/(.*)/(.*).SampleSheet.csv", ss)
        else:
            r = re.match(f"{runDir.rstrip('/')}/(.*)/STARK/(.*).SampleSheet.csv", ss)
        if r is None:
            continue
        elif r.group(1) == r.group(2):  # checks if (.*) == (.*)
            return ss
    return "NO_SAMPLESHEET_FOUND"


def createRunningFile(run, serviceName):
    """Creates a running file indicating the service is executing."""
    with open(osj(run, serviceName + "Running.txt"), "w+") as file:
        file.write(f"# [{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}] {os.path.basename(run)} running with {serviceName}\n")


def launch(run, serviceName, containersFile, montage, image, launchCommand, configFile):
    """Launches the Docker container using Docker Compose with the necessary configurations."""
    createRunningFile(run, serviceName)

    # Read the image and launch command from the configuration file
    if not launchCommand:
        launchCommand = readconfig(configFile, serviceName, "launch")
    image = readconfig(configFile, serviceName, "image")[0]  # assuming only one image in the list

    # Ensure launchCommand is valid
    if not launchCommand:
        raise ValueError(f"[ERROR] No 'launch' command specified for service '{serviceName}' in the config file.")

    # Find samplesheet and .bed files
    samplesheet = find_any_samplesheet(run)
    if samplesheet == "NO_SAMPLESHEET_FOUND":
        raise FileNotFoundError(f"[ERROR] find_any_samplesheet() couldn't find any samplesheet in run {run}.")
    
    bed = findAnyBed(run)
    md5 = getMd5(run)
    containerName = f"{serviceName}-{md5}-NAME-{os.path.basename(run)}"

    # Check if Docker Compose file exists
    compose_file = os.getenv('DOCKER_COMPOSE_FILE', 'docker-compose.yml')
    if not os.path.exists(compose_file):
        raise FileNotFoundError(f"Could not find Docker Compose file: {compose_file}")

    # Build the Docker Compose command
    cmd = f"docker compose -f {compose_file} run --rm --name={containerName} -v {run}:{run} {montage} {image} {' '.join(launchCommand)} --run={run}"

    print(f"Running command: {cmd}")
    subprocess.call(cmd, shell=True)

    # Create log file after execution
    createContainerFile(containersFile, run, containerName)


def myoptions():
    """Argument parser for command line options."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--run", type=str, default="", help="Run to launch", dest="run")
    parser.add_argument("-s", "--service", type=str, default="", help="Service name", dest="serviceName")
    parser.add_argument("-log", "--log", type=str, default="", help="Path for the log file", dest="containersFile")
    parser.add_argument("-v", "--volume", type=str, default="", help="Docker volumes to add", dest="montage")
    parser.add_argument("-i", "--image", type=str, default="", help="Docker image to use", dest="image")
    parser.add_argument("-l", "--launchcommand", type=str, default="", help="Command to launch inside the container", dest="launchCommand")
    parser.add_argument("-c", "--config", type=str, default="", help="Config file to read from", dest="configFile")
    return parser.parse_args()


if __name__ == "__main__":
    args = myoptions()

    # Launch the container with parsed arguments
    launch(
        args.run,
        args.serviceName,
        args.containersFile,
        args.montage,
        args.image,
        args.launchCommand,
        args.configFile,
    )