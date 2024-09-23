# -*- coding: utf-8 -*-
from os.path import join as osj
import os
import subprocess
import logging as log
import json


def launch(container_name, launch_arguments):
    module_config = osj(os.environ["DOCKER_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    howard_config = osj(os.environ["DOCKER_MODULE_CONFIG"], "howard_config.json")
    if not os.path.isfile(module_config):
        log.error(f"{module_config} do not exist, primordial file, check its existence")
        raise ValueError(module_config)
    elif not os.path.isfile(howard_config):
        log.error(f"{howard_config} do not exist, primordial file, check its existence")
        raise ValueError(howard_config)
      
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        howard_image = data["howard_image"]
    log.info(f"Using {howard_image}")
    
    command_list = [
        "docker", 
        "run", 
        "--rm", 
        "--name", 
        container_name, 
        "--env", 
        f"http_proxy={os.environ["http_proxy"]}", 
        "--env", 
        f"ftp_proxy={os.environ["ftp_proxy"]}", 
        "--env", 
        f"https_proxy={os.environ["https_proxy"]}",
        "-v",
        f"{os.environ["HOST_TMP"]}:{os.environ["DOCKER_TMP"]}",
        "-v",
        f"{os.environ["HOST_DATABASES"]}:/databases/",
        "-v",
        f"{os.environ["HOST_SERVICES"]}:{os.environ["DOCKER_SERVICES"]}",
        "-v",
        f"{os.environ["HOST_CONFIG"]}:{os.environ["DOCKER_CONFIG"]}",
        "-v",
        f"{howard_config}:/tools/devel/config/config.json",
        howard_image
        ]
    
    command_list = command_list + launch_arguments
    log.debug(" ".join(command_list))
    subprocess.call(command_list, universal_newlines=True)