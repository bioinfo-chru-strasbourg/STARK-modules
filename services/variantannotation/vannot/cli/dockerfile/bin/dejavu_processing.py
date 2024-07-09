# -*- coding: utf-8 -*-
import glob
from os.path import join as osj
import os
import logging as log
import shutil

import howard_launcher


def convert_vcf_parquet(run_informations):
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
        
        launch_minimalize_arguments = ["minimalize", "--input", vcf_file, "--output", vcf_minimalize, "--minimalize_info", "--minimalize_samples"]
        launch_parquet_arguments = ["process", "--input", vcf_minimalize, "--output", output_parquet, "--calculations=BARCODE", "--explode_infos", '--explode_infos_fields=barcode', '--query=SELECT \"#CHROM\", POS, ID, REF, ALT, QUAL, FILTER, INFO, barcode FROM variants']
        howard_launcher.launch(container_name, launch_minimalize_arguments)
        howard_launcher.launch(container_name, launch_parquet_arguments)

    shutil.rmtree(run_informations["tmp_analysis_folder"])

def calculate_dejavu(run_informations):
    container_name = f"VANNOT_dejavu_{run_informations["run_name"]}"
    parquet_db_project = run_informations["parquet_db_howard_project_folder"]

    inner_dejavu_output_parquet = osj(run_informations["parquet_db_howard_folder"], f"dejavu.{run_informations["run_application"]}.parquet")

    dejavu_output_parquet_hdr = osj(run_informations["parquet_db_folder"], f"dejavu.{run_informations["run_application"]}.parquet.hdr")
    dejavu_output_parquet = osj(run_informations["parquet_db_folder"], f"dejavu.{run_informations["run_application"]}.parquet")
    dejavu_previous_output_parquet_hdr = osj(run_informations["parquet_db_folder"], f"previous.dejavu.{run_informations["run_application"]}.parquet.hdr")
    dejavu_previous_output_parquet = osj(run_informations["parquet_db_folder"], f"previous.dejavu.{run_informations["run_application"]}.parquet")

    if os.path.isfile(dejavu_previous_output_parquet):
        os.remove(dejavu_previous_output_parquet)
    if os.path.isfile(dejavu_previous_output_parquet_hdr):
        os.remove(dejavu_previous_output_parquet_hdr)
    if os.path.isfile(dejavu_output_parquet):
        os.rename(dejavu_output_parquet, dejavu_previous_output_parquet)
    if os.path.isfile(dejavu_output_parquet_hdr):
        os.rename(dejavu_output_parquet_hdr, dejavu_previous_output_parquet_hdr)

    sample_count = len(glob.glob(osj(run_informations["parquet_db_project_folder"], "*", "*.parquet")))
    
    launch_query_arguments = ["query", "--input", parquet_db_project, "--query", f"SELECT \"#CHROM\", POS, REF, ALT, count(barcode)/({sample_count}*2) AS FREQ FROM variants GROUP BY \"#CHROM\", POS, REF, ALT", "--output", inner_dejavu_output_parquet]
    howard_launcher.launch(container_name, launch_query_arguments)