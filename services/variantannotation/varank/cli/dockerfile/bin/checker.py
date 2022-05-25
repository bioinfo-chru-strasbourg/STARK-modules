# -*- coding: utf-8 -*-
import os
import json
import logging as log
from os.path import join as osj


def absolute_folder_path(path):
    if (
        os.path.isabs(path)
        and os.path.isdir(path)
        and not path.startswith(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_REPOSITORY"
            ]
        )
    ):
        return path
    elif (
        not os.path.isabs(path)
        and os.path.isdir(path)
        and not os.path.abspath(path).starstwith(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_REPOSITORY"
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

        threads = data.get("threads")
        if threads == "all":
            log.info(
                f"{threads} available threads will be used according to {varank_json}"
            )
        elif not threads.isnumeric:
            log.warning(f"Not a number, all available threads will be used")
        elif threads.isnumeric:
            log.info(f"{threads} threads will be used according to {varank_json}")

        for i in data["databases"]:
            if i == "alamutDB":
                if not os.path.isfile(data["databases"]["alamutDB"]):
                    log.error(f"{i} file do not exists, please modify in {varank_json}")
                    raise ValueError(data["databases"]["alamutDB"])
                else:
                    alamutdb_path = data["databases"]["alamutDB"]
                    log.info(f"Alamut database file {alamutdb_path} exists")

            elif i == "OMIM":
                api_token = data["databases"]["OMIM"]["api_token"]
                OMIM_1 = data["databases"]["OMIM"]["OMIM_1"]
                OMIM_2 = data["databases"]["OMIM"]["OMIM_2"]
                # if api_token == "...":
                #     log.error(
                #         f"{i} license was not parametered, please modify in {varank_json}"
                #     )
                #     raise ValueError(data["databases"]["alamutDB"]["api_token"])
                if not os.path.isfile(OMIM_1):
                    log.error(
                        f"{OMIM_1} file do not exists, please modify the path in {varank_json} or download the latest version with the omim command"
                    )
                    raise ValueError(OMIM_1)
                else:
                    log.info(f"OMIM database file {OMIM_1} exists")

                if not os.path.isfile(OMIM_2):
                    log.error(
                        f"{OMIM_2} file do not exists, please modify the path in {varank_json} or download the latest version with the omim command"
                    )
                    raise ValueError(OMIM_2)
                else:
                    log.info(f"OMIM database file {OMIM_2} exists")

        alamut_license_checker(alamutdb_path)


def alamut_license_checker(alamutdb_path):
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

            if line.startswith("File="):
                line = line[5:]
                if not os.path.isfile(line):
                    alamut_license_tmp = alamut_license + ".tmp"
                    with open(alamut_license, "r") as read_file:
                        with open(alamut_license_tmp, "w") as write_file:
                            for line2 in read_file:
                                if line2.startswith("File="):
                                    write_file.write(f"File={alamutdb_path}\n")
                                else:
                                    write_file.write(line2)
                    os.remove(alamut_license)
                    os.rename(alamut_license_tmp, alamut_license)

                    log.warning(
                        f"AlamutDB database file {line} was not found, it was replace with the data in the varank_config.json"
                    )

    with open(alamut_license, "r") as read_file:
        for line in read_file:
            line = line.rstrip()

            for key, value in license_parameters.items():
                if line == f"{key}=...":
                    log.error(
                        f"AlamutDB {value} is not set in the license, please check the file in the alamut config file {alamut_license}"
                    )
                    raise ValueError(line)


if __name__ == "__main__":
    pass
