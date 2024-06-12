# -*- coding: utf-8 -*-
import glob
from os.path import join as osj
import os
import subprocess
import logging as log
import shutil
import json
import multiprocessing


def folder_initialisation(run_informations):
    input_files = glob.glob(
        osj(run_informations["archives_project_folder"], "*.final.vcf*")
    ) + glob.glob(osj(run_informations["archives_project_folder"], "*", "*.final.vcf*"))

    if len(input_files) != 0:
        if not os.path.isdir(
            osj(run_informations["archives_project_folder"], "VCF", "DEFAULT")
        ):
            os.makedirs(
                osj(run_informations["archives_project_folder"], "VCF", "DEFAULT", ""),
                0o775,
            )
        archives_default_vcf_folder = osj(
            run_informations["archives_project_folder"], "VCF", "DEFAULT", ""
        )

        for input_file in input_files:
            subprocess.run(
                [
                    "rsync",
                    "-rp",
                    input_file,
                    archives_default_vcf_folder,
                ]
            )
            os.remove(input_file)

    vcf_file_list = glob.glob(
        osj(run_informations["archives_project_folder"], "VCF", "*", "*.final.vcf*")
    )

    if os.path.isdir(run_informations["tmp_analysis_folder"]):
        log.error(f"{run_informations["tmp_analysis_folder"]} already existing, be careful that the analysis is not already running")
        raise ValueError(run_informations["tmp_analysis_folder"])
    else:
        os.mkdir(run_informations["tmp_analysis_folder"])

    for processed_vcf_file in vcf_file_list:
        subprocess.run(
            [
                "rsync",
                "-rp",
                processed_vcf_file,
                run_informations["tmp_analysis_folder"],
            ]
        )

    vcf_file_to_analyse = glob.glob(
        osj(run_informations["tmp_analysis_folder"], "*.final.vcf*")
    )
    for vcf_file in vcf_file_to_analyse:
        howard_launcher(run_informations, vcf_file)
        os.remove(vcf_file)

def run_initialisation(run_informations):
    vcf_file_list = glob.glob(
        osj(run_informations["archives_run_folder"], "*.final.vcf*")
    )
    if os.path.isdir(run_informations["tmp_analysis_folder"]):
        shutil.rmtree(run_informations["tmp_analysis_folder"])
        os.mkdir(run_informations["tmp_analysis_folder"])
    else:
        os.mkdir(run_informations["tmp_analysis_folder"])

    for processed_vcf_file in vcf_file_list:
        subprocess.run(
            [
                "rsync",
                "-rp",
                processed_vcf_file,
                osj(run_informations["tmp_analysis_folder"], ""),
            ]
        )

    vcf_file_to_analyse = glob.glob(
        osj(run_informations["tmp_analysis_folder"], "*.final.vcf*")
    ) 
    for vcf_file in vcf_file_to_analyse:
        howard_launcher(run_informations, vcf_file)
        os.remove(vcf_file)


def howard_launcher(run_informations, vcf_file):
    log.info(f"Launching for {vcf_file}")
    logfile = osj(
        run_informations["tmp_analysis_folder"], f"VA_{os.path.basename(vcf_file.split(".")[0])}.log"
    )

    if run_informations["output_format"] != None:
        output_file = osj(
            run_informations["tmp_analysis_folder"], f"VA_{os.path.basename(vcf_file).split(".")[0]}.{run_informations["output_format"]}"
        )
    else:
        output_file = osj(
            run_informations["tmp_analysis_folder"], f"VA_{os.path.basename(vcf_file)}"
        )

    if run_informations["parameters_file"] == None:
        param_file = osj(os.environ["DOCKER_MODULE_CONFIG"], "configfiles", f"param.{run_informations["run_platform_application"]}.json")
    else:
        param_file = run_informations["parameters_file"]
    
    module_config = osj(os.environ["DOCKER_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        howard_image = data["howard_image"]
    
    container_name = f"VA_{run_informations['run_name']}_{os.path.basename(vcf_file).split('.')[0]}"

    log.info("Generating results files")
    
    with open(logfile, "w") as f:
        subprocess.call(
            [
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
                f"{os.environ["HOST_CONFIG"]}:{os.environ["DOCKER_CONFIG"]}",
                howard_image,
                "annotation",
                "--input",
                vcf_file,
                "--output",
                output_file,
                "--param",
                param_file,
                "--assembly",
                run_informations["assembly"],
            ],
            stdout=f,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
        )
    os.chmod(logfile, 0o777)


def cleaner(run_informations):
    log.info("Cleaning the analysis folder")
    results_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*"))

    if os.path.isdir(run_informations["archives_results_folder"]):
        shutil.rmtree(run_informations["archives_results_folder"])
        os.mkdir(run_informations["archives_results_folder"])
    else:
        os.mkdir(run_informations["archives_results_folder"])

    for results_file in results_files:
        os.chmod(results_file, 0o777)
        shutil.move(results_file, run_informations["archives_results_folder"])
    log.info("Moved generated tsv files")

    with open(
        osj(run_informations["archives_results_folder"], "VACopyComplete.txt"),
        mode="a",
    ):
        pass

    shutil.rmtree(run_informations["tmp_analysis_folder"])
    log.info("Deleted temporary analysis folder")


if __name__ == "__main__":
    pass
