#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@Goal: 
- Listener for STARK18 modules
- Launch analysis if triggers are passed

@Author: Victor Grentzinger (2020)
"""

import argparse
import doctest
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
import time

from os.path import join as osj
from datetime import datetime
import pandas as pd

sys.path.insert(1, "/app/bin/config")
from launcher import launch

def create_container_file(containers_file, run, container_name):
    with open(osj(containers_file, f"{container_name}.log"), "w+") as file:
        file.write(f"RUN: {os.path.basename(run)}\n")
        file.write(f"FOLDER: {run}\n")
        file.write(f"EXEC_DATE: {datetime.now().strftime('%d%m%Y-%H%M%S')}\n")
        file.write(f"ID: {container_name}\n")

def get_md5(run):
    return hashlib.md5(run.encode()).hexdigest()

def create_running_file(run, service_name):
    with open(osj(run, f"{service_name}Running.txt"), "w+") as file:
        file.write(f"# [{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}] {os.path.basename(run)} running with {service_name}\n")

def check_projects(service_project, sample_projects_list):
    return any(service_project in project.split("-")[0] for project in sample_projects_list)

def check_group(service_group, sample_groups_list):
    return any(service_group in group.split("-")[0] for group in sample_groups_list)

def assert_file_exists_and_is_readable(file_path):
    assert os.path.isfile(file_path) and os.access(file_path, os.R_OK), f"[ERROR] File {file_path} doesn't exist or isn't readable"

def get_sample_project_from_samplesheet(samplesheet_path):
    if samplesheet_path == "NO_SAMPLESHEET_FOUND":
        return []
    assert_file_exists_and_is_readable(samplesheet_path)

    df = pd.read_csv(samplesheet_path)
    if 'Sample_Project' in df.columns:
        sample_projects = df['Sample_Project'].dropna().tolist()
    else:
        sample_projects = df[df.columns[1]].tolist() if 'Investigator Name' in df.columns else []
    
    return sample_projects

def find_any_samplesheet(run_dir, from_res_dir=False):
    pattern = "*SampleSheet.csv"
    for ss in glob.glob(osj(run_dir, '**', pattern), recursive=True):
        ss = ss.strip()
        if from_res_dir:
            r = re.match(run_dir.rstrip("/") + "/(.*)/(.*).SampleSheet.csv", ss)
        else:
            r = re.match(run_dir.rstrip("/") + "/(.*)/STARK/(.*).SampleSheet.csv", ss)
        if r and r.group(1) == r.group(2):
            return ss
    return "NO_SAMPLESHEET_FOUND"

def check_tags(tag, sample_tags_list, analysis_tags_list):
    """
    >>> check_tags("!POOL", ['SEX#F!PLUGAPP#POOL!', 'SEX#F!', 'SEX#M!PLUGAPP#POOL!', 'SEX#M!', 'SEX#F!'], ['APP#DIAG.DI', 'APP#DIAG.DI', 'APP#DIAG.DI', 'APP#DIAG.DI', 'APP#DIAG.DI'])
    False
    >>> check_tags("POOL", ['SEX#F!PLUGAPP#POOL!', 'SEX#F!', 'SEX#M!PLUGAPP#POOL!', 'SEX#M!', 'SEX#F!'], ['APP#DIAG.DI', 'APP#DIAG.DI', 'APP#DIAG.DI', 'APP#DIAG.DI', 'APP#DIAG.DI'])
    True
    """
    tags_to_check = set(re.split("[!#]", tag.strip('!')) for tag_list in (sample_tags_list, analysis_tags_list) for tag in tag_list)
    if tag.startswith("!"):
        return tag[1:] not in tags_to_check
    else:
        return tag in tags_to_check

def get_analysis_tags_from_sample_list(sample_list, run):
    analysis_tags_list = []
    for s in sample_list:
        with open(osj(run, s, "STARK", f"{s}.analysis.tag"), "r") as tags_file:
            analysis_tags_list.extend(tag.strip() for tag in tags_file)
    return analysis_tags_list

def get_sample_tags_from_sample_list(sample_list, run):
    sample_tags_list = []
    for s in sample_list:
        with open(osj(run, s, "STARK", f"{s}.tag"), "r") as tags_file:
            sample_tags_list.extend(tag.strip() for tag in tags_file)
    return sample_tags_list

def get_sample_list_from_run_path(run):
    return [os.path.basename(os.path.normpath(patient)) for patient in glob.glob(osj(run, "*/")) 
            if os.path.exists(osj(patient, "STARKCopyComplete.txt")) or os.path.exists(osj(patient, os.path.basename(os.path.normpath(patient)) + ".SampleSheet.csv"))]

def complete(run, service_name):
    return not bool(glob.glob(osj(run, f"{service_name}Complete.txt")))

def running(run, service_name):
    return not bool(glob.glob(osj(run, f"{service_name}Running.txt")))

def failed(run, service_name):
    return not bool(glob.glob(osj(run, f"{service_name}Failed.txt")))

def check_and_condition(value, service_name, run):
    return all(check_triggers({and_key: value[and_key]}, service_name, run) for and_key in value.keys())

def check_or_condition(value, service_name, run):
    return any(check_triggers({or_key: value[or_key]}, service_name, run) for or_key in value.keys())

def check_file_condition(value, run):
    return all(
        [
            not running(run, file[1:-11]) if file.startswith("!") and "Running.txt" in file else
            not complete(run, file[1:-12]) if file.startswith("!") and "Complete.txt" in file else
            not complete(run, "STARK") if file.startswith("!") and "STARKCopyComplete.txt" in file else
            running(run, file[:-11]) if "Running.txt" in file else
            complete(run, file[:-12]) if "Complete.txt" in file else
            failed(run, file[:-12]) if "Failed.txt" in file else
            complete(run, "STARK") if "STARKCopyComplete.txt" in file else
            False for file in value
        ]
    )

def check_tags_condition(value, run):
    sample_list = get_sample_list_from_run_path(run)
    sample_tags_list = get_sample_tags_from_sample_list(sample_list, run)
    analysis_tags_list = get_analysis_tags_from_sample_list(sample_list, run)
    return all(check_tags(tag, sample_tags_list, analysis_tags_list) for tag in value)

def check_group_or_project_condition(key, value, run):
    run_group = os.path.abspath(run).split("/")[-3] if key == "group" else os.path.abspath(run).split("/")[-2]
    return run_group in value

def check_triggers(jconfig, service_name, run):
    if not jconfig:
        return True

    for key, value in jconfig.items():
        if key.startswith("AND"):
            if check_and_condition(value, service_name, run):
                return True
        elif key.startswith("OR"):
            if check_or_condition(value, service_name, run):
                return True
        elif key == "file":
            if check_file_condition(value, run):
                return True
        elif key == "tags":
            if check_tags_condition(value, run):
                return True
        elif key in ("group", "project"):
            if check_group_or_project_condition(key, value, run):
                return True
    return False

def get_data_from_json(json_file):
    with open(json_file, 'r') as json_analysis_file:
        return json.load(json_analysis_file)

def check_triggers_from_config(json_file, service_name, run):
    jconfig = get_data_from_json(json_file)['services'][service_name]['triggers']
    return check_triggers(jconfig, service_name, run)

def verify_triggers(run, service_name, json_file):
    return check_triggers_from_config(json_file, service_name, run)

def get_run_name(stark_complete):
    assert stark_complete.endswith('/STARKCopyComplete.txt'), "[ERROR]"
    return stark_complete[:-22]

def get_stark_copy_complete(group_input, days):
    cmd = f"find {group_input}/*/*/STARKCopyComplete.txt -maxdepth 5 -mtime -{days} 2>/dev/null"
    out = subprocess.check_output(cmd, shell=True).decode("utf-8").splitlines()
    return out

def main(group_input_list, service_name, json_file, days, delay, containers_file, config_file):
    while True:
        for group_input in group_input_list:
            stark_complete_list = get_stark_copy_complete(group_input, days)
            for stark_complete in stark_complete_list:
                run_dir = os.path.dirname(stark_complete)
                if any(os.path.basename(path.strip()) == "STARK" for path in glob.glob(osj(run_dir, "*"))):
                    run = get_run_name(stark_complete)
                    if verify_triggers(run, service_name, json_file):
                        print(f"[INFO] Launching {service_name} analysis for run {run}")
                        launch(run, service_name, containers_file, os.getenv('MICROSERVICE_MONTAGE'), os.getenv('MICROSERVICE_IMAGE'), os.getenv('MICROSERVICE_LAUNCH'), config_file, os.getenv('MICROSERVICE_REPOSITORY'))
                else:
                    print(f"[INFO] Run {stark_complete} ignored as no STARK directory was found")
        time.sleep(60.0 * delay)

def myoptions():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default=" ", help="list of runs to listen to : <PATH_RUN_1>,<PATH_RUN_2>...", dest='listInput')
    parser.add_argument("-n", "--servicename", type=str, default="TEST", help="name of the service to use", dest='serviceName')
    parser.add_argument("-j", "--json", type=str, default="config/listener.json", help="path to the triggers file listener.json", dest='jsonFile')
    parser.add_argument("-c", "--config", type=str, default="", help="path to the config file listener.conf", dest='configFile')
    parser.add_argument("-f", "--containersfile", type=str, default="", help="path to the container's file folder", dest='containersFile')
    parser.add_argument("-t", "--nbDaysBack", type=int, default=30, help="folder older than x days", dest='days')
    parser.add_argument("-r", "--repository", type=int, default=5, help="path to repository", dest='repository')
    return parser.parse_args()

if __name__ == "__main__":
    doctest.testmod()
    args = myoptions()
    group_input_list = args.listInput.split(",")
    main(group_input_list, args.serviceName, args.jsonFile, args.days, args.delay, args.containersFile, args.configFile, args.repository)
