# -*- coding: utf-8 -*-
from multiprocessing.sharedctypes import Value
import os
import json
import logging as log
from os.path import join as osj


def absolute_folder_path(path):
    if os.path.isabs(path) and os.path.isdir(path):
        return path
    elif not os.path.isabs(path) and os.path.isdir(path):
        return os.path.abspath(path)
    else:
        raise ValueError(path)


def absolute_run_path(path):
    if (
        os.path.isabs(path)
        and os.path.isdir(path)
        and path.startswith(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_REPOSITORY"
            ]
        )
    ):
        return path
    elif (
        not os.path.isabs(path)
        and os.path.isdir(path)
        and os.path.abspath(path).startswith(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_REPOSITORY"
            ]
        )
    ):
        return os.path.abspath(path)
    else:
        raise ValueError(path)


def varank_config_json_checker():
    varank_json = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG"],
        "variantannotation",
        "varank",
        "varank_config.json",
    )
    with open(varank_json, "r") as read_file:
        data = json.load(read_file)
        for i in data["proxy"]:
            if data["proxy"][i] == "...":
                log.warning(
                    f"{i} is not defined, if you are using a proxy, please modify values in {varank_json}"
                )

        for i in data["databases"]:
            value = data["databases"][i]
            if value == "...":
                log.error(f"{i} path is not defined, please modify in {varank_json}")
                raise ValueError(value)
            elif not os.path.isfile(value):
                log.error(f"{i} file do not exists, please modify in {varank_json}")
                raise ValueError(value)

        threads = data.get("threads")
        if threads == "all":
            log.info(f"{threads} threads will be used according to {varank_json}")
        elif not threads.isnumeric:
            log.warning(f"Not a number, all available threads will be used")
        elif threads.isnumeric:
            log.info(f"{threads} threads will be used according to {varank_json}")


def alamut_license_checker():
    alamut_license = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG"],
        "variantannotation",
        "varank",
        "alamut-batch-license",
        "alamut-batch.ini",
    )

    license_parameters = {
        "File": "database path",
        "Name": "database file name",
        "Institution": "institution",
        "LicenceKey": "license key",
        "User": "username",
    }

    with open(alamut_license, "r") as read_file:
        for line in read_file:
            line = line.rstrip()

            for key, value in license_parameters.items():
                if line == f"{key}=...":
                    log.error(
                        f"AlamutDB {value} is not set in the license, please check the file in the alamut config file {alamut_license}"
                    )
                    raise ValueError(line)

            if line.startswith("File="):
                line = line[5:]
                if not os.path.isfile(line):
                    log.error(
                        f"AlamutDB database file {line} was not found, please check the file in the alamut config file {alamut_license}"
                    )
                    raise ValueError(line)
                else:
                    log.info(f"AlamutDB database file {line} exists")


if __name__ == "__main__":
    pass
