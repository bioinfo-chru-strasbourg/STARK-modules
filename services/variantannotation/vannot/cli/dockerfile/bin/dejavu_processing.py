# -*- coding: utf-8 -*-
import glob
from os.path import join as osj
import os
import subprocess
import logging as log
import shutil
import json

def define_howard_image():
    module_config = osj(os.environ["DOCKER_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        howard_image = data["howard_image"]
    log.info(f"Using {howard_image}")
    return howard_image

def convert_vcf_parquet(run_informations):
    howard_image = define_howard_image()

    if not os.path.isdir(run_informations["parquet_db_run_folder"]):
        os.makedirs(run_informations["parquet_db_run_folder"], 0o775)
    else:
        shutil.rmtree(run_informations["parquet_db_run_folder"])
        os.makedirs(run_informations["parquet_db_run_folder"], 0o775)

    if not os.path.isdir(run_informations["tmp_analysis_folder"]):
        os.makedirs(run_informations["tmp_analysis_folder"], 0o775)
    else:
        shutil.rmtree(run_informations["tmp_analysis_folder"])
        os.makedirs(run_informations["tmp_analysis_folder"], 0o775)


    for vcf_file in glob.glob(osj(run_informations["archives_run_folder"], "*.final.vcf.gz")):
        container_name = f"VANNOT_dejavu_{run_informations["run_name"]}_{os.path.basename(vcf_file).split(".")[0]}"
        vcf_minimalize = osj(run_informations["tmp_analysis_folder"], f"{os.path.basename(vcf_file).split(".")[0]}.mini.parquet")
        output_parquet = osj(run_informations["parquet_db_howard_run_folder"], f"{os.path.basename(vcf_file).split(".")[0]}.parquet")

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
                f"{os.environ["HOST_SERVICES"]}:{os.environ["DOCKER_SERVICES"]}",
                howard_image,
                "minimalize",
                "--input",
                vcf_file,
                "--output",
                vcf_minimalize,
                "--minimalize_info",
                "--minimalize_samples",
            ],
            universal_newlines=True,
        )
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
                f"{os.environ["HOST_SERVICES"]}:{os.environ["DOCKER_SERVICES"]}",
                howard_image,
                "process",
                "--input",
                vcf_minimalize,
                "--output",
                output_parquet,
                "--calculations=BARCODE",
                "--explode_infos",
                '--explode_infos_fields=barcode',
                '--query=SELECT \"#CHROM\", POS, ID, REF, ALT, QUAL, FILTER, INFO, barcode FROM variants',  
            ],
            universal_newlines=True,
        )

    shutil.rmtree(run_informations["tmp_analysis_folder"])

def calculate_dejavu(run_informations):
    howard_image = define_howard_image()
    container_name = f"VANNOT_dejavu_{run_informations["run_name"]}"
    parquet_db_project = run_informations["parquet_db_howard_project_folder"]

    inner_dejavu_output_parquet = osj(run_informations["parquet_db_howard_project_folder"], f"dejavu.{run_informations["run_application"]}.parquet")
    inner_dejavu_output_vcf = osj(run_informations["parquet_db_howard_project_folder"], f"dejavu.{run_informations["run_application"]}.vcf.gz")

    dejavu_output_parquet = osj(run_informations["parquet_db_project_folder"], f"dejavu.{run_informations["run_application"]}.parquet")
    dejavu_output_vcf = osj(run_informations["parquet_db_project_folder"], f"dejavu.{run_informations["run_application"]}.vcf.gz")
    dejavu_previous_output_parquet = osj(run_informations["parquet_db_project_folder"], f"previous.dejavu.{run_informations["run_application"]}.parquet")
    dejavu_previous_output_vcf = osj(run_informations["parquet_db_project_folder"], f"previous.dejavu.{run_informations["run_application"]}.vcf.gz")

    # if os.path.isfile(dejavu_output_parquet):
    #     os.rename(dejavu_output_parquet, dejavu_previous_output_parquet)
    # if os.path.isfile(dejavu_output_vcf):
    #     os.rename(dejavu_output_vcf, dejavu_previous_output_vcf)

    sample_count = len(glob.glob(osj(run_informations["parquet_db_project_folder"], "*", "*.parquet")))
    
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
            f"{os.environ["HOST_SERVICES"]}:{os.environ["DOCKER_SERVICES"]}",
            howard_image,
            "query",
            "--input",
            parquet_db_project,
            "--query",
            f"SELECT \"#CHROM\", POS, REF, ALT, count(barcode)/({sample_count}*2) AS FREQ FROM variants GROUP BY \"#CHROM\", POS, REF, ALT",
            "--output",
            inner_dejavu_output_parquet,
        ],
        universal_newlines=True,
    )
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
            f"{os.environ["HOST_SERVICES"]}:{os.environ["DOCKER_SERVICES"]}",
            howard_image,
            "convert",
            "--input",
            inner_dejavu_output_parquet,
            "--output",
            inner_dejavu_output_vcf,
        ],
        universal_newlines=True,
    )