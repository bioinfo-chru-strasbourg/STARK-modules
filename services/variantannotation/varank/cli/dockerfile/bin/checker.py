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
                if api_token == "...":
                    log.info(
                        f"{i} license was not parametered, please modify asap in {varank_json}"
                    )
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
                        f"AlamutDB database file {line} was not found in the license, it was replace with the data in the varank_config.json"
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
    pattern = run_informations["vcf_pattern"]

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


def configfile(run_informations):
    configfile = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG"],
        "variantannotation",
        "varank",
        "configfiles",
        f"configfile.{run_informations['run_platform_application']}",
    )

    extann_config_file = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG"],
        "variantannotation",
        "varank",
        "configfiles",
        "extann_config_file.tsv",
    )
    extann_directory = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG"],
        "variantannotation",
        "varank",
        "configfiles",
        "extanns",
    )
    vcf_fields = False

    family_list = []
    tmp_configfile = osj(
        os.path.dirname(configfile), f"tmp.{os.path.basename(configfile)}"
    )

    if not os.path.isfile(configfile):
        log.error(f"There is no configfile linked to this application, please create one and launch the analysis again")
        raise ValueError(configfile)

    with open(configfile, "r") as read_file:
        for line in read_file.readlines():
            line = line.strip()
            if re.match(r"^-vcfFields:", line):
                vcf_fields = True

    with open(tmp_configfile, "w") as write_file, open(
        configfile, "r"
    ) as read_file, open(extann_config_file, "r") as read_file2:
        for line in read_file.readlines():
            writted = False
            line = line.strip()

    with open(tmp_configfile, "w") as write_file, open(
        configfile, "r"
    ) as read_file, open(extann_config_file, "r") as read_file2:
        for line in read_file.readlines():
            writted = False
            line = line.strip()

            if re.match(r"-vcfInfo:|#-vcfInfo:", line):
                if re.match(r"#-vcfInfo:", line):
                    line = re.sub("#", "", line)
                if re.search("no", line):
                    line = re.sub("no", "yes", line)
                if vcf_fields is False:
                    write_file.write(
                        line
                        + "\n"
                        + '-vcfFields: "FindByPipelines GenotypeConcordance POOL_F_Depth POOL_M_Depth POOL_F_base_counts POOL_M_base_counts BARCODE trio_variant_type"'
                        + "\n"
                    )
                    writted = True
                elif vcf_fields is True:
                    write_file.write(line + "\n")
                    writted = True
                vcf_fields = True

            if re.match(r"-metrics:|#-metrics:", line):
                if re.match(r"#-metrics:", line):
                    line = re.sub("#", "", line)
                if re.search("us", line):
                    line = re.sub("us", "fr", line)
                write_file.write(line + "\n")
                writted = True

            if re.match(r"#-uniprot:", line):
                line = re.sub("#", "", line)
                write_file.write(line + "\n")
                writted = True

            if re.match(r"#-refseq:", line):
                line = re.sub("#", "", line)
                write_file.write(line + "\n")
                writted = True

            if re.match(r"fam[0-9]+:", line):
                fam = line.split(" ")
                fam = fam[1:]
                for id in fam:
                    family_list.append(id)

            if re.match(r"-extann:|#-extann:", line):
                extann_files = line.split('"')
                for extann_file in reversed(extann_files):
                    if not extann_file.startswith("/"):
                        extann_files.remove(extann_file)

                if os.path.basename(configfile) != "configfile.default":
                    if len(extann_files) != 0:
                        for extann_file in extann_files:
                            extann_file_list = extann_file.split(" ")
                    else:
                        extann_file_list = []

                    for extann_file in reversed(extann_file_list):
                        if re.search(r"^.*$", extann_file):
                            extann_file_list.remove(extann_file)

                    next(read_file2, 1)
                    for line2 in read_file2.readlines():
                        line2 = line2.strip()
                        line2 = line2.split("\t")
                        extann_platform_application = line2[0] + "." + line2[1]
                        extann_name = osj(extann_directory, line2[2])
                        extann_used = line2[3]
                        configfile_platform_application = run_informations[
                            "run_platform_application"
                        ]

                        if (
                            extann_platform_application
                            == configfile_platform_application
                            and extann_used == "yes"
                        ):
                            extann_file_list.append(extann_name)

                    if len(extann_file_list) == 0:
                        log.info(
                            f"Not using any additional extann file for {configfile}"
                        )

                    for extann_file in reversed(extann_file_list):
                        if os.path.isfile(extann_file) and len(extann_file_list) != 0:
                            log.info(
                                f"Using the following for your analysis : {extann_file}"
                            )
                            continue
                        elif (
                            not os.path.isfile(extann_file)
                            and len(extann_file_list) != 0
                        ):
                            extann_file_list.remove(extann_file)
                            log.error(
                                f"Not found any additional extann(s) file(s) at {extann_file} directory, please check  your configfile"
                            )

                    if len(extann_file_list) == 0:
                        line = '#-extann:\t\t""'
                    else:
                        line = '-extann:\t\t"' + " ".join(extann_file_list) + '"'
                        write_file.write(line + "\n")
                        writted = True

            proxy_user = os.environ["http_proxy"].split(":")[1][2:]
            proxy_passwd = os.environ["http_proxy"].split(":")[2].split("@")[0]
            proxy_server = os.environ["http_proxy"].split(":")[2].split("@")[1]
            proxy_port = os.environ["http_proxy"].split(":")[3]

            if re.match(r"#-proxyUser:", line):
                if re.search(r'"[a-zA-Z0-9]+"', line):
                    line = re.sub('"[a-zA-Z0-9]+"', '"' + proxy_user + '"', line)
                else:
                    line = re.sub('""', '"' + proxy_user + '"', line)
                line = re.sub("#", "", line)
                write_file.write(line + "\n")
                writted = True

            if re.match(r"#-proxyPasswd:", line):
                if re.search(r'"[a-zA-Z0-9]+"', line):
                    line = re.sub(
                        '"[a-zA-Z0-9]+"',
                        '"' + proxy_passwd + '"',
                        line,
                    )
                else:
                    line = re.sub('""', '"' + proxy_passwd + '"', line)
                line = re.sub("#", "", line)
                write_file.write(line + "\n")
                writted = True

            if re.match(r"#-proxyServer:", line):
                if re.search(r'"[a-zA-Z0-9]+"', line):
                    line = re.sub(
                        '"[a-zA-Z0-9]+"',
                        '"' + proxy_server + '"',
                        line,
                    )
                else:
                    line = re.sub('""', '"' + proxy_server + '"', line)
                line = re.sub("#", "", line)
                write_file.write(line + "\n")
                writted = True

            if re.match(r"#-proxyPort:", line):
                if re.search(r'"[a-zA-Z0-9]+"', line):
                    line = re.sub('"[a-zA-Z0-9]+"', '"' + proxy_port + '"', line)
                else:
                    line = re.sub('""', '"' + proxy_port + '"', line)
                line = re.sub("#", "", line)
                write_file.write(line + "\n")
                writted = True

            if re.match(r"-proxyUser:", line):
                if re.search(r'"[a-zA-Z0-9]+"', line):
                    line = re.sub('"[a-zA-Z0-9]+"', '"' + proxy_user + '"', line)
                else:
                    line = re.sub('""', '"' + proxy_user + '"', line)
                write_file.write(line + "\n")
                writted = True

            if re.match(r"-proxyPasswd:", line):
                if re.search(r'"[a-zA-Z0-9]+"', line):
                    line = re.sub(
                        '"[a-zA-Z0-9]+"',
                        '"' + proxy_passwd + '"',
                        line,
                    )
                else:
                    line = re.sub('""', '"' + proxy_passwd + '"', line)
                write_file.write(line + "\n")
                writted = True

            if re.match(r"-proxyServer:", line):
                if re.search(r'"[a-zA-Z0-9]+"', line):
                    line = re.sub(
                        '"[a-zA-Z0-9]+"',
                        '"' + proxy_server + '"',
                        line,
                    )
                else:
                    line = re.sub('""', '"' + proxy_server + '"', line)
                write_file.write(line + "\n")
                writted = True

            if re.match(r"-proxyPort:", line):
                if re.search(r'"[a-zA-Z0-9]+"', line):
                    line = re.sub('"[a-zA-Z0-9]+"', '"' + proxy_port + '"', line)
                else:
                    line = re.sub('""', '"' + proxy_port + '"', line)
                write_file.write(line + "\n")
                writted = True

            if writted is False:
                write_file.write(line + "\n")

    os.remove(configfile)
    os.rename(tmp_configfile, configfile)
    os.chmod(configfile, 0o775)


def logfile(run_informations):
    varank_error = True
    logfile = osj(run_informations["archives_project_folder"], "VaRank.log")

    with open(logfile, "r") as read_file:
        for line in read_file.readlines():
            if re.match(r"^\.\.\.VaRank is done with the analysis", line):
                varank_error = False
                log.info(
                    f"VaRank is done with the analysis, non redundant will be launched {logfile}"
                )

    if varank_error is True:
        log.error(f"Unexpected error just happened, please check {logfile}")


if __name__ == "__main__":
    pass
