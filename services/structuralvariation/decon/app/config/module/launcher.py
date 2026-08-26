##########################################################################
# Launcher Version:			3.0
# Description:				Launcher to run Snakemake module
##########################################################################

# DEV version 0.1 : 10/11/2021
# INT version 0.1 : 17/03/2022
# PROD version 1 : 03/06/2022
# Authoring : Thomas LAVAUX

# PROD version 2 : 16/06/2022 changelog
# yaml files can be defined by groupe_name + project_name or groupe_name only

# PROD version 3.0 : 28/11/2023 changelog
# docker compose to run containers

# PROD version 3.0 : 17/08/2026 changelog
# hg38 ready

################## Context ##################
# type python launcher.py -h for help
# ex of command: python launcher.py -r run
####################################

# From listener.py
# launch(run, serviceName, containersFile, os.getenv('MICROSERVICE_MONTAGE'), os.getenv('MICROSERVICE_IMAGE'), os.getenv('MICROSERVICE_LAUNCH'), configFile, os.getenv('MICROSERVICE_REPOSITORY'))
# run = path of the run to analyse
# serviceName = name of the service (ex DECON)
# containersFile = part name of a log file
# configFile = config json with .conf ext containing launch command (= MICROSERVICE_LAUNCH) and image name of the docker ( = MICROSERVICE_IMAGE)
# MICROSERVICE_MONTAGE = additional volume to mount
# MICROSERVICE_IMAGE = cf. supra
# MICROSERVICE_LAUNCH = cf. supra
# MICROSERVICE_REPOSITORY = path to the repository folder

################## Import libraries ##################

import os
import re
import glob
import subprocess
import json
import argparse
import doctest
from datetime import datetime
from os.path import join as osj

date_time = datetime.now().strftime("%Y%m%d-%H%M%S")


# Function to compare versions
def version_gt(version1, version2):
    """Compare two dotted version strings numerically (e.g. '20.10.7' > '9.2' is True,
    which a naive string comparison would get wrong)."""

    def parse(version):
        parts = []
        for part in version.split("."):
            digits = "".join(ch for ch in part if ch.isdigit())
            parts.append(int(digits) if digits else 0)
        return tuple(parts)

    return parse(version1) > parse(version2)


# Function to get Docker version
def get_docker_version():
    try:
        docker_version_output = (
            subprocess.check_output(["docker", "--version"]).decode().strip()
        )
        docker_version = docker_version_output.split()[2].split(",")[0]
        return docker_version
    except (subprocess.CalledProcessError, FileNotFoundError, IndexError):
        return None


def readconfig(configFile, serviceName, configkey):
    """Function to extract a specific key configuration variable from a json configuration file"""
    with open(configFile, "r") as f:
        json_config = json.load(f)
    outputconfig = json_config["services"][serviceName][configkey]
    return outputconfig


# Matches an ASSEMBLY=value key, however it's prefixed (ex "#[INFO]        ASSEMBLY=hg19" in a STARK
# run report, or a plain "ASSEMBLY=hg19" line). Anchored on "=" right after the word so it can't match
# prose that merely mentions "genome assembly" or "reference assembly" without an "=" following it.
ASSEMBLY_PATTERN = re.compile(r"(?<![A-Za-z0-9_])ASSEMBLY\s*=\s*([^\s;#]+)", re.IGNORECASE)


def find_assembly_from_config(run):
    """
    Search the run folder for STARK '*.config' analysis-report files (the per-run ini-like
    dumps STARK writes, ex KLA2602897_20260723-084645.config) and return the ASSEMBLY value
    (ex hg19, hg38) from the most recently modified file that has one.
    Returns None if run doesn't exist, no config file is found, or none contain an ASSEMBLY key.
    """
    if not run or not os.path.isdir(run):
        return None

    # Most runs have the report config directly at the top level; fall back to a recursive
    # search only if that turns up nothing, to keep the common case fast.
    config_files = glob.glob(os.path.join(run, "*.config"))
    if not config_files:
        config_files = glob.glob(os.path.join(run, "**", "*.config"), recursive=True)

    for config_file in sorted(config_files, key=os.path.getmtime, reverse=True):
        try:
            with open(config_file, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    match = ASSEMBLY_PATTERN.search(line)
                    if match:
                        return match.group(1).strip()
        except OSError:
            continue
    return None


def find_yaml_config(yaml_path, group_name, project_name, assembly=None):
    """
    Return the most specific existing yaml config file, trying in order:
      {group}_{project}_{assembly}.yaml
      {group}_{project}.yaml
      {group}_{assembly}.yaml
      {group}.yaml
    so an assembly-specific override is preferred when present, but everything falls back
    to the pre-existing group/project convention when there's no assembly-specific file
    (or no assembly could be determined at all).
    Returns None if none of these exist.
    """
    candidates = []
    if assembly:
        candidates.append(f"{yaml_path}/{group_name}_{project_name}_{assembly}.yaml")
    candidates.append(f"{yaml_path}/{group_name}_{project_name}.yaml")
    if assembly:
        candidates.append(f"{yaml_path}/{group_name}_{assembly}.yaml")
    candidates.append(f"{yaml_path}/{group_name}.yaml")

    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    return None


def launch(
    run,
    serviceName,
    containersFile=None,
    montage=None,
    image=None,
    launchCommand=None,
    configFile=None,
    microserviceRepo=None,
    assembly=None,
):
    """Function to start a docker container with a specific command"""
    if configFile:
        launchCommand = readconfig(configFile, serviceName, "launch")
        image = readconfig(configFile, serviceName, "image")

    # group_name/project_name previously stayed undefined whenever run was falsy, which
    # crashed the "if group_name and project_name" check below with a NameError.
    group_name = None
    project_name = None
    if run:
        containerName = f"{serviceName}_{date_time}_{os.path.basename(run)}"
        run_parts = run.split("/")
        if len(run_parts) > 4:
            group_name, project_name = run_parts[3], run_parts[4]
    else:
        containerName = f"{serviceName}_{date_time}"

    if group_name and project_name:
        yaml_path = (
            f"{os.getenv('DOCKER_STARK_MODULE_SUBMODULE_INNER_FOLDER_CONFIG')}/cli"
        )
        # ASSEMBLY (hg19/hg38/...) is auto-detected from the run's STARK '*.config' report
        # file unless explicitly overridden; an assembly-specific yaml is preferred when one
        # exists, falling back to the pre-existing group/project convention otherwise.
        detected_assembly = assembly or find_assembly_from_config(run)
        if detected_assembly:
            print(f"[INFO] Using ASSEMBLY={detected_assembly} for yaml selection")
        else:
            print("[INFO] No ASSEMBLY value found for this run, falling back to group/project yaml only")
        yaml_config_file = find_yaml_config(yaml_path, group_name, project_name, detected_assembly)
    else:
        yaml_config_file = None

    COMPOSE_PATH = (
        f"{os.getenv('DOCKER_STARK_MODULE_SUBMODULE_INNER_FOLDER_CONFIG')}/listener/"
    )

    docker_version = get_docker_version()
    if not docker_version:
        # DOCKER_COMMAND was previously left unset on this path, so building cmd below would
        # crash with a NameError instead of failing on this clear, actionable message.
        print("[ERROR] Docker is not installed or not reachable, aborting launch.")
        return

    # Docker Compose V2 ("docker compose", a CLI plugin) ships with Docker Engine 20.10+;
    # older engines need the standalone V1 binary ("docker-compose"). The two branches were
    # previously swapped - verify this matches what's actually installed on your host.
    if version_gt(docker_version, "20"):
        DOCKER_COMMAND = "docker compose"
    else:
        DOCKER_COMMAND = "docker-compose"
    print(f"Using {DOCKER_COMMAND} for Docker commands.")

    if yaml_config_file and os.path.exists(yaml_config_file):
        cmd = f"{DOCKER_COMMAND} -f {COMPOSE_PATH}/STARK.docker-compose.yml run --rm --name={containerName} {image} '{launchCommand} --config run={run} --configfile {yaml_config_file}'"
        print(cmd)
    else:
        cmd = f"{DOCKER_COMMAND} -f {COMPOSE_PATH}/STARK.docker-compose.yml run --rm --name={containerName} {image} '{launchCommand} --config run={run}'"
        print(cmd)
    subprocess.call(cmd, shell=True)


def myoptions():
    """
    *arg parser*
    *return options*
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r", "--run", type=str, default="", help="Run to launch", dest="run"
    )
    parser.add_argument(
        "-s", "--service", type=str, default="", help="Service name", dest="serviceName"
    )
    parser.add_argument(
        "-log",
        "--log",
        type=str,
        default="",
        help="Path for the log file",
        dest="containersFile",
    )
    parser.add_argument(
        "-v",
        "--volume",
        type=str,
        default="",
        help="Docker volumes to add",
        dest="montage",
    )
    parser.add_argument(
        "-i", "--image", type=str, default="", help="Docker image to use", dest="image"
    )
    parser.add_argument(
        "-l",
        "--launchcommand",
        type=str,
        default="",
        help="Command to launch inside the container",
        dest="launchCommand",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default="",
        help="Config file to read from",
        dest="configFile",
    )
    parser.add_argument(
        "-repo",
        "--repo",
        type=str,
        default="",
        help="Microservice repository name",
        dest="microserviceRepo",
    )
    parser.add_argument(
        "-a",
        "--assembly",
        type=str,
        default="",
        help="Override the auto-detected ASSEMBLY value (ex hg19, hg38) used to pick a yaml config file",
        dest="assembly",
    )
    return parser.parse_args()


if __name__ == "__main__":
    doctest.testmod()
    args = myoptions()
    launch(
        args.run,
        args.serviceName,
        args.containersFile,
        args.montage,
        args.image,
        args.launchCommand,
        args.configFile,
        args.microserviceRepo,
        args.assembly or None,
    )

