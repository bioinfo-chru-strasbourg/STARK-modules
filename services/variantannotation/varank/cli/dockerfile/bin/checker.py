# -*- coding: utf-8 -*-
import os
import json
import logging as log
import requests
import subprocess
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


def alamutdb_checker():
    alamutdb_folder = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DATABASES"],
        "alamutDB",
    )

    if not os.path.isdir(alamutdb_folder):
        os.mkdir(alamutdb_folder)
        log.warning(
            "AlamutDB do not exist, the newer version will be downloaded, the process can take few hours"
        )
        alamut_source = osj(alamutdb_folder, "source.txt")
        alamut_req = requests.get(
            "https://downloads.interactive-biosoftware.com/", "html.parser"
        )
        with open(alamut_source, "w") as write_file:
            write_file.write(alamut_req.text)

        with open(alamut_source, "r") as read_file:
            for line in read_file:
                line.rstrip()
                if (
                    "http://downloads.interactive-biosoftware.com/sources/files/Alamut-batch-standalone/"
                    in line
                ):
                    alamut_download_url = line.split("'")[1]
                    version = alamut_download_url.split("/")[6].split("-")[2]
                    year = version.split(".")[0]
                    month = version.split(".")[1]
                    day = version.split(".")[2]
                    version = year + month + day

        subprocess.run(
            [
                "wget",
                alamut_download_url,
                osj(alamutdb_folder, version),
            ]
        )

    if not os.path.isfile(osj(alamutdb_folder, "STARK.database")):
        with open(osj(alamutdb_folder, "STARK.database"), "w") as write_file:
            write_file.write("coucou")


# def omim_checker():


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
            if i == "OMIM":
                print(data["databases"]["OMIM"])
            elif i == "alamutDB":
                if not os.path.isfile(data["databases"]["alamutDB"]):
                    log.error(f"{i} file do not exists, please modify in {varank_json}")
                    raise ValueError(data["databases"]["alamutDB"])
                else:
                    alamutdb_path = data["databases"]["alamutDB"]

        # for i in data["databases"]:
        #     print(data["databases"][i])

        # for i in data["databases"]:
        #     value = data["databases"][i]
        #     if value == "...":
        #         log.error(f"{i} path is not defined, please modify in {varank_json}")
        #         raise ValueError(value)
        #     elif not os.path.isfile(value):
        #         log.error(f"{i} file do not exists, please modify in {varank_json}")
        #         raise ValueError(value)

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
                else:
                    log.info(f"AlamutDB database file {line} exists")

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
