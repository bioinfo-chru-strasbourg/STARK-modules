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
    # alamut_licence()
    # varank_configuration()
    # varank_listener_configuration()


def global_config():
    if not os.path.isdir(osj(config_varank, "alamut-batch-license")):
        os.mkdir(osj(config_varank, "alamut-batch-license"))

    if not os.path.isdir(osj(config_varank, "configfiles")):
        os.mkdir(osj(config_varank, "configfiles"))

    if not os.path.isdir(osj(config_varank, "listener")):
        os.mkdir(osj(config_varank, "listener"))


def alamut_licence():
    if not os.path.isfile(osj(config_varank, "alamut-batch.ini")):
        shutil.copy(
            osj("/setup", "alamut-batch.ini"),
            osj(config_varank, "alamut-batch-license"),
        )


# def varank_configuration():

# def varank_listener_configuration():


if __name__ == "__main__":
    config_varank = os.environ[
        "DOCKER_STARK_MODULE_SUBMODULE_SERVICE_SETUP_INNER_FOLDER_CONFIG"
    ]
    main()
