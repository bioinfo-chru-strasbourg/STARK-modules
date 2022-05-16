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
from os.path import join as osj


def main():
    global_config()
    alamut_licence()
    varank_configuration()
    varank_listener_configuration()


def global_config():
    main_folders = ["alamut-batch-license", "configfiles", "listener"]
    for i in main_folders:
        if not os.path.isdir(osj(config_varank, i)):
            os.mkdir(osj(config_varank, i))


def alamut_licence():
    if_not_file_copy("alamut-batch-license", "alamut-batch.ini")


def varank_configuration():
    if not os.path.isdir(osj(config_varank, "configfiles", "extanns")):
        os.mkdir(osj(config_varank, "configfiles", "extanns"))

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


def if_not_file_copy(inner_folder, inner_file):
    if not os.path.isfile(osj(inner_folder, inner_file)):
        shutil.copy(
            osj("/setup", inner_file),
            osj(config_varank, inner_folder),
        )


if __name__ == "__main__":
    config_varank = os.environ[
        "DOCKER_STARK_MODULE_SUBMODULE_SERVICE_SETUP_INNER_FOLDER_CONFIG"
    ]
    main()
