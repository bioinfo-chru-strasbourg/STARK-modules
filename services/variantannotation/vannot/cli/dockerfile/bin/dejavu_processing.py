# -*- coding: utf-8 -*-
import glob
from os.path import join as osj
import os
import logging as log
import shutil
import subprocess

import howard_launcher


def convert_vcf_parquet(run_informations):
    log.info(f"Generating new dejavu database for {run_informations["run_platform_application"]}")
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

    vcf_files = glob.glob(osj(run_informations["archives_run_folder"], "*.vcf.gz"))
    for vcf_file in vcf_files:
        subprocess.call(["rsync", "-rvt", vcf_file, run_informations["tmp_analysis_folder"]], universal_newlines=True)
    
    for vcf_file in glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz")):
        container_name = f"VANNOT_dejavu_{run_informations["run_name"]}_{os.path.basename(vcf_file).split(".")[0]}"
        vcf_minimalize = osj(run_informations["tmp_analysis_folder"], f"{os.path.basename(vcf_file).split(".")[0]}.mini.parquet")
        output_parquet = osj(run_informations["tmp_analysis_folder"], f"{os.path.basename(vcf_file).split(".")[0]}.parquet")
        
        launch_minimalize_arguments = ["minimalize", "--input", vcf_file, "--output", vcf_minimalize, "--minimalize_info", "--minimalize_samples"]
        launch_parquet_arguments = ["process", "--input", vcf_minimalize, "--output", output_parquet, "--calculations=BARCODE", "--explode_infos", '--explode_infos_fields=barcode', '--query=SELECT \"#CHROM\", POS, ID, REF, ALT, QUAL, FILTER, INFO, barcode FROM variants']
        log.info("Minimalizing vcfs")
        howard_launcher.launch(container_name, launch_minimalize_arguments)
        log.info("Converting to parquet format")
        howard_launcher.launch(container_name, launch_parquet_arguments)
        os.remove(vcf_minimalize)
        os.remove(vcf_minimalize + ".hdr")

    parquet_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.parquet*"))
    
    for parquet_file in parquet_files:
        if "hdr" not in os.path.basename(parquet_file):
            sample_dejavu_db_folder = osj(run_informations["parquet_db_run_folder"], f"SAMPLE={".".join(os.path.basename(parquet_file).split(".")[:-1])}")
        else:
            sample_dejavu_db_folder = osj(run_informations["parquet_db_run_folder"], f"SAMPLE={".".join(os.path.basename(parquet_file).split(".")[:-2])}")

        if not os.path.isdir(sample_dejavu_db_folder):
            os.mkdir(sample_dejavu_db_folder)
        subprocess.call(["rsync", "-rvt", parquet_file, sample_dejavu_db_folder], universal_newlines=True)

    shutil.rmtree(run_informations["tmp_analysis_folder"])

def calculate_dejavu(run_informations):
    container_name = f"VANNOT_dejavu_{run_informations["run_name"]}"
    parquet_db_project = run_informations["parquet_db_howard_folder"]
    project = run_informations["run_application"]

    inner_dejavu_root_folder = os.path.dirname(parquet_db_project)
    inner_dejavu_output_parquet = osj(inner_dejavu_root_folder, f"dejavu.{run_informations["run_application"]}.parquet")

    output_root_folder = os.path.dirname(run_informations["parquet_db_folder"])
    dejavu_output_parquet_hdr = osj(output_root_folder, f"dejavu.{run_informations["run_application"]}.parquet.hdr")
    dejavu_output_parquet = osj(output_root_folder, f"dejavu.{run_informations["run_application"]}.parquet")
    dejavu_previous_output_parquet_hdr = osj(output_root_folder, f"previous.dejavu.{run_informations["run_application"]}.parquet.hdr")
    dejavu_previous_output_parquet = osj(output_root_folder, f"previous.dejavu.{run_informations["run_application"]}.parquet")

    if os.path.isfile(dejavu_previous_output_parquet):
        os.remove(dejavu_previous_output_parquet)
    if os.path.isfile(dejavu_previous_output_parquet_hdr):
        os.remove(dejavu_previous_output_parquet_hdr)
    if os.path.isfile(dejavu_output_parquet):
        os.rename(dejavu_output_parquet, dejavu_previous_output_parquet)
    if os.path.isfile(dejavu_output_parquet_hdr):
        os.rename(dejavu_output_parquet_hdr, dejavu_previous_output_parquet_hdr)

    sample_count = len(glob.glob(osj(run_informations["parquet_db_project_folder"], "*", "*", "*.parquet")))
    
    log.info("Calculating new frequencies")
    launch_query_arguments = ["query", "--input", parquet_db_project, "--query", f"SELECT '#CHROM', POS, REF, ALT, count(barcode)/({sample_count}*2) AS FREQ FROM variants WHERE PROJECT='{project}' GROUP BY '#CHROM', POS, REF, ALT", "--output", inner_dejavu_output_parquet]
    howard_launcher.launch(container_name, launch_query_arguments)