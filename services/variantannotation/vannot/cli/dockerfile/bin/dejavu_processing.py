# -*- coding: utf-8 -*-
import glob
from os.path import join as osj
import os
import logging as log
import shutil
import subprocess
import commons

import howard_launcher


def convert_vcf_parquet(run_informations):
    threads = commons.get_threads("threads_dejavu")
    memory = commons.get_memory("memory_dejavu")

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
        
        launch_minimalize_arguments = ["minimalize", "--input", vcf_file, "--output", vcf_minimalize, "--minimalize_info", "--minimalize_samples", "--threads", threads, "--memory", memory]
        launch_parquet_arguments = ["process", "--input", vcf_minimalize, "--output", output_parquet, "--calculations=BARCODE", "--explode_infos", '--explode_infos_fields=barcode', '--query=SELECT \"#CHROM\", POS, ID, REF, ALT, QUAL, FILTER, INFO, barcode FROM variants', "--threads", threads, "--memory", memory]
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
    threads = commons.get_threads("threads_dejavu")
    memory = commons.get_memory("memory_dejavu")

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
    allelecount = "ALLELECOUNT"
    hetcount = "HETCOUNT"
    homcount = "HOMCOUNT"
    allelefreq = "ALLELEFREQ"
    query = f"SELECT \"#CHROM\", POS, ANY_VALUE(ID) AS ID, REF, ALT, ANY_VALUE(QUAL) AS QUAL, ANY_VALUE(FILTER) AS FILTER, ANY_VALUE(INFO) AS INFO, sum(CAST(barcode AS INT)) AS {allelecount}, count(barcode) FILTER(barcode=1) AS {hetcount}, count(barcode) FILTER(barcode=2) AS {homcount}, sum(CAST(barcode AS INT))/({sample_count}*2) AS {allelefreq} FROM variants WHERE PROJECT='{project}' GROUP BY \"#CHROM\", POS, REF, ALT"
    
    launch_query_arguments = ["query", "--input", parquet_db_project, "--query", query , "--output", inner_dejavu_output_parquet, "--threads", threads, "--memory", memory]

    howard_launcher.launch(container_name, launch_query_arguments)
    os.remove(dejavu_output_parquet_hdr)
    with open(dejavu_output_parquet_hdr, "w") as writefile:
        writefile.write('##fileformat=VCFv4.2\n')
        writefile.write('##INFO=<ID=ALLELECOUNT,Number=.,Type=Float,Description=\"allele count annotation\">\n')
        writefile.write('##INFO=<ID=HETCOUNT,Number=.,Type=Float,Description=\"heterozygot count annotation\">\n')
        writefile.write('##INFO=<ID=HOMCOUNT,Number=.,Type=Float,Description=\"homozygot count annotation\">\n')
        writefile.write('##INFO=<ID=ALLELEFREQ,Number=.,Type=Float,Description=\"allele frequency annotation\">\n')
        writefile.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t{allelecount}\t{hetcount}\t{homcount}\t{allelefreq}\n')
