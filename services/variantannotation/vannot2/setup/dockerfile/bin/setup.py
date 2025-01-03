# !/usr/bin/python
# !/bin/sh
# !/usr/tcl/bin
# -*- coding: utf-8 -*-
"""
###############################
##							 ##
##		VARANK ANALYSIS		 ##
## need alamut-batch license ##
##	Author : Mateusz RAUCH	 ##
##							 ##
###############################
"""
import os
import shutil
import subprocess
import json
from datetime import datetime
from os.path import join as osj


def main():
    if os.path.isfile(osj(config_varank, "SetupComplete.txt")):
        os.remove(osj(config_varank, "SetupComplete.txt"))
    global_config()
    alamut_licence()
    varank_configuration()
    varank_listener_configuration()
    varank_services()
    alamutdb_checker()


def global_config():
    main_folders = ["alamut-batch-license", "configfiles", "listener"]
    for i in main_folders:
        if_not_folder_create(config_varank, i)
    if_not_file_copy("", "varank_config.json")


def alamut_licence():
    if_not_file_copy("alamut-batch-license", "alamut-batch.ini")


def varank_configuration():
    if_not_folder_create(osj(config_varank, "configfiles"), "extanns")

    configuration_files = [
        "configfile.default",
        "non_redundant_config.txt",
        "extann_config_file.tsv",
        "deleted_samples_config.tsv",
    ]
    for i in configuration_files:
        if_not_file_copy("configfiles", i)


def varank_listener_configuration():
    configuration_files = ["launcher.py", "varank.json"]
    for i in configuration_files:
        if_not_file_copy("listener", i)


def varank_services():
    main_folders = ["Archives", "listener", "logs"]
    for i in main_folders:
        if_not_folder_create(services_varank, i)

    stark_module = osj(services_varank[:-6], "STARK.module")
    stark_submodule = osj(services_varank, "STARK.module")

    if not os.path.isfile(stark_module):
        json = {
            "code": "v",
            "name": "Variant Annotation",
            "fullname": "Variant Annotation Services",
            "release": os.environ["DOCKER_STARK_MODULE_SUBMODULE_RELEASE"],
            "available": "true",
        }

        with open(stark_module, "w") as write_file:
            json.dump(json, write_file, indent=4)

    if not os.path.isfile(stark_submodule):
        json = {
            "submodules": {
                "varank": {
                    "description": "Only VaRank annotation tool",
                    "website": "https://www.lbgi.fr/VaRank/",
                    "services": {
                        "cli": {
                            "code": "cli",
                            "name": "VaRank cli",
                            "type": "cli",
                            "description": "A VaRank cli",
                            "release": os.environ[
                                "DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_RELEASE"
                            ],
                            "available": "true",
                        },
                        "listener": {
                            "code": "listener",
                            "name": "VaRank listener",
                            "type": "listener",
                            "description": "A VaRank listener service is started as a daemon, listening for new STARK analyzed run (new folder in output/repository) and well configured (STARKComplete.txt without VaRankRunning.txt or VaRankComplete.txt), and check triggers to launch or not analysis",
                            "release": os.environ[
                                "DOCKER_STARK_MODULE_SUBMODULE_SERVICE_LISTENER_RELEASE"
                            ],
                            "available": "true",
                        },
                        "setup": {
                            "code": "setup",
                            "name": "VaRank setup",
                            "type": "setup",
                            "description": "A VaRank setup, prepare all the folders and configurations files for VaRank to work. It will also download AlamutDB if not available",
                            "release": os.environ[
                                "DOCKER_STARK_MODULE_SUBMODULE_SERVICE_SETUP_RELEASE"
                            ],
                            "available": "true",
                        },
                    },
                }
            }
        }

        with open(stark_submodule, "w") as write_file:
            json.dump(json, write_file, indent=4)


def if_not_file_copy(inner_folder, inner_file):
    if not os.path.isfile(osj(config_varank, inner_folder, inner_file)):
        shutil.copy(
            osj("/setup", inner_file),
            osj(config_varank, inner_folder),
        )


def if_not_folder_create(path, new_folder):
    if not os.path.isdir(osj(path, new_folder)):
        os.mkdir(osj(path, new_folder))


if __name__ == "__main__":
    config_varank = os.environ[
        "DOCKER_STARK_MODULE_SUBMODULE_SERVICE_SETUP_INNER_FOLDER_CONFIG"
    ]
    services_varank = os.environ[
        "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES"
    ]
    main()
