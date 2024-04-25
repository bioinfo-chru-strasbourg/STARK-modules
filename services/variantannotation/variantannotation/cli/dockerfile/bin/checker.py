# -*- coding: utf-8 -*-
import os
import json
import logging as log
import glob
from os.path import join as osj
import commons
import re


def absolute_folder_path(path):
    if (
        os.path.isabs(path)
        and os.path.isdir(path)
        and not path.startswith(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARIANTANNOTATION_FOLDER_REPOSITORY"
            ]
        )
    ):
        return path
    elif (
        not os.path.isabs(path)
        and os.path.isdir(path)
        and not os.path.abspath(path).starstwith(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARIANTANNOTATION_FOLDER_REPOSITORY"
            ]
        )
    ):
        return os.path.abspath(path)
    else:
        raise ValueError(path)


def absolute_run_path(path):
    if (
        os.path.isabs(path)
        and os.path.isdir(path)
        and path.startswith(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARIANTANNOTATION_FOLDER_REPOSITORY"
            ]
        )
    ):
        return path
    elif (
        not os.path.isabs(path)
        and os.path.isdir(path)
        and os.path.abspath(path).startswith(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARIANTANNOTATION_FOLDER_REPOSITORY"
            ]
        )
    ):
        return os.path.abspath(path)
    else:
        raise ValueError(path)

def depository_checker(run_informations):
    run_repository = run_informations["run_repository"]
    run_depository = run_informations["run_depository"]

    if not os.path.isdir(run_depository):
        log.warning(
            f"Specified depository folder doesn't exists, maybe it was sent to the Archives ? Creating a new folder with all subfolders tree."
        )
        run_repository_sample_folders = glob.glob(osj(run_repository, "*", ""))
        os.makedirs(run_depository, 0o775)

        for sample_folder in run_repository_sample_folders:
            sample_folder = osj(run_depository, os.path.basename(sample_folder[:-1]))
            if not os.path.isdir(sample_folder):
                os.mkdir(sample_folder, 0o755)


def pattern_checker(run_informations):
    run_repository = run_informations["run_repository"]
    pattern = run_informations["run_pattern"]

    for element in pattern:
        vcf_files = glob.glob(osj(run_repository, element))
        if len(vcf_files) == 0 and element != commons.default_pattern:
            log.error(
                f"There is no vcf files with the specified pattern {element}, please check your command-line"
            )
            raise ValueError(element)

        elif len(vcf_files) == 0 and element == commons.default_pattern:
            log.error(
                f"There is no vcf files with the default STARK analysis pattern {element}, please check the analysis integrity"
            )
            raise ValueError(element)


# def configfile(run_informations):

def logfile(run_informations):
    variantannotation_error = True
    logfile = osj(run_informations["run_processing_folder"], "variantannotation.log")

    with open(logfile, "r") as read_file:
        for line in read_file.readlines():
            if re.match(r"^\.\.\.variantannotation is done with the analysis", line):
                variantannotation_error = False
                log.info(
                    f"variantannotation is done with the analysis, non redundant will be launched {logfile}"
                )

    if variantannotation_error is True:
        log.error(f"Unexpected error just happened, please check {logfile}")


if __name__ == "__main__":
    pass