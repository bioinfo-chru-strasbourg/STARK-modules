# -*- coding: utf-8 -*-
import glob
from os.path import join as osj
import os
import subprocess
import logging as log
import shutil
import json

def define_howard_image():
    module_config = osj(os.environ["DOCKER_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        howard_image = data["howard_image"]
    log.info(f"Using {howard_image}")
    return howard_image

def convert_vcf_parquet(run_informations):
    print(run_informations["parquet_db_project_folder"])
    howard_image = define_howard_image()

def calculate_dejavu(run_informations):
    print(run_informations["parquet_db_project_folder"])
    howard_image = define_howard_image()