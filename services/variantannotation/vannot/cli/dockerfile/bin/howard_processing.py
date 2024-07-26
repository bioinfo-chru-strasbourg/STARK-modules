# -*- coding: utf-8 -*-
import glob
from os.path import join as osj
import os
import subprocess
import logging as log
import shutil
import json
import re

import howard_launcher


def folder_initialisation(run_informations):
    input_files = glob.glob(
        osj(run_informations["archives_project_folder"], "*vcf*")
    ) + glob.glob(osj(run_informations["archives_project_folder"], "*", "*vcf*"))

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
        osj(run_informations["archives_project_folder"], "VCF", "*", "*vcf*")
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
        osj(run_informations["tmp_analysis_folder"], "*vcf*")
    )
    for vcf_file in vcf_file_to_analyse:
        howard_proc(run_informations, vcf_file)
        os.remove(vcf_file)

def run_initialisation(run_informations):
    vcf_file_list = glob.glob(
        osj(run_informations["archives_run_folder"], "*vcf*")
    )
    if os.path.isdir(run_informations["tmp_analysis_folder"]):
        log.info("Cleaning temporary analysis folder")
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
        osj(run_informations["tmp_analysis_folder"], "*vcf*")
    )
    
    module_config = osj(os.environ["DOCKER_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        annotations_to_keep = data["keep_vcf_info"][run_informations["run_platform_application"]]

    for vcf_file in vcf_file_to_analyse:
        cleaned_vcf = cleaning_annotations(vcf_file, annotations_to_keep)
        howard_proc(run_informations, cleaned_vcf)
        os.remove(cleaned_vcf)

def cleaning_annotations(vcf_file, annotations_to_keep):
    actual_info_fields = subprocess.run(["zgrep", "##INFO", vcf_file], capture_output=True, text=True)
    actual_info_fields = actual_info_fields.stdout.strip().split("##")
    if not vcf_file.endswith(".vcf.gz"):
        subprocess.call(["/tools/bin/bgzip", vcf_file])
        vcf_file = vcf_file + ".gz"

    cleaned_vcf = osj(os.path.dirname(vcf_file), "cleaned_" + os.path.basename(vcf_file)[:-3])
    info_to_keep = []
    for i in annotations_to_keep:
        for j in actual_info_fields:
            if re.search(i, j):
                info_to_keep.append("INFO/" + j.split(",")[0].split("=")[-1])

    if len(info_to_keep) == 0:
        log.info("No annotations to keep were found, deleting all annotations")
        info_to_keep = "INFO"
    else:
        log.info(f"Keeping following annotations: {' '.join(info_to_keep)} for sample {os.path.basename(vcf_file)}")
        info_to_keep = "^" + ",".join(info_to_keep)

    cmd = ["/tools/bin/bcftools", "annotate", "-x"]
    cmd.append(info_to_keep)
    cmd.append(vcf_file)

    with open(cleaned_vcf, "w") as output:
        subprocess.call(cmd, stdout=output, universal_newlines=True)
    os.remove(vcf_file)
        
    subprocess.call(["/tools/bin/bgzip", cleaned_vcf], universal_newlines=True)
    renamed_clean = osj(os.path.dirname(cleaned_vcf), os.path.basename(cleaned_vcf).replace("cleaned_", "") + ".gz")
    cleaned_vcf = cleaned_vcf + ".gz"
    os.rename(cleaned_vcf, renamed_clean)

    return renamed_clean

def howard_proc(run_informations, vcf_file):
    log.info(f"Launching for {vcf_file}")
    # logfile = osj(
    #     run_informations["tmp_analysis_folder"], f"VANNOT_annotate_{os.path.basename(vcf_file.split(".")[0])}.log"
    # )

    if run_informations["output_format"] != None:
        output_file = osj(
            run_informations["tmp_analysis_folder"], f"VANNOT_{os.path.basename(vcf_file).split(".")[0]}.{run_informations["output_format"]}"
        )
    else:
        output_file = osj(
            run_informations["tmp_analysis_folder"], f"VANNOT_{os.path.basename(vcf_file)}"
        )

    if run_informations["parameters_file"] == None:
        for option in [run_informations["run_platform_application"], run_informations["run_platform"], "default"]:
            configfile = osj(os.environ["DOCKER_MODULE_CONFIG"], "configfiles", f"param.{option}.json")
            if os.path.isfile(configfile):
                break
            else:
                continue

    log.info(f"Using {configfile} as parameter for HOWARD analysis")
    if not os.path.isfile(configfile):
        log.error("param.default.json not found, please check your config directory")

    container_name = f"VANNOT_annotate_{run_informations['run_name']}_{os.path.basename(vcf_file).split('.')[0]}"
    launch_annotate_arguments = ["annotation", "--input", vcf_file, "--output", output_file, "--param", configfile, "--assembly", run_informations["assembly"]]

    log.info("Annotating input files with HOWARD")
    
    howard_launcher.launch(container_name, launch_annotate_arguments)
    convert_to_final_tsv(run_informations, output_file)

def merge_vcf_files(run_informations):
    log.info("Merging all vcfs into one")
    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))
    logfile = osj(run_informations["tmp_analysis_folder"], f"VANNOT_merged_{run_informations["run_name"]}.log")
    output_file = osj(run_informations["tmp_analysis_folder"], f"merged_VANNOT_{run_informations["run_name"]}.vcf.gz")
    
    for vcf_file in vcf_files:
        with open(logfile, "a") as f:
            subprocess.call(["/tools/bin/tabix", vcf_file], stdout=f, stderr=subprocess.STDOUT, universal_newlines=True)

    cmd = ["/tools/bin/bcftools", "merge", "-o", output_file, "-O", "z"]
    for vcf_file in vcf_files:
        cmd.append(vcf_file)
    with open(logfile, "a") as f:
        subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT, universal_newlines=True)
    

def convert_to_final_tsv(run_informations, input_file):
    log.info("Converting output file into readable tsv")
    # logfile = osj(run_informations["tmp_analysis_folder"], f"VANNOT_convert_{os.path.basename(input_file.split(".")[0])[7:]}.log")
    container_name = f"VANNOT_convert_{run_informations['run_name']}_{os.path.basename(input_file).split('.')[0]}"
    output_file = osj(run_informations["tmp_analysis_folder"], f"{os.path.basename(input_file).split(".")[0]}.tsv")
    
    launch_convert_arguments = ["convert", "--input", input_file, "--output", output_file, "--explode_infos"]
    howard_launcher.launch(container_name, launch_convert_arguments)

def cleaner(run_informations):
    log.info("Moving results from temporary folder")
    results_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*tsv")) + glob.glob(osj(run_informations["tmp_analysis_folder"], "*log")) + glob.glob(osj(run_informations["tmp_analysis_folder"], "*vcf.gz"))

    if os.path.isdir(run_informations["archives_results_folder"]):
        log.info("Removing old results folder in archives")
        shutil.rmtree(run_informations["archives_results_folder"])
        os.mkdir(run_informations["archives_results_folder"])
    else:
        os.mkdir(run_informations["archives_results_folder"])

    for results_file in results_files:
        os.chmod(results_file, 0o777)
        log.info(f"Moving {results_file} to {run_informations["archives_results_folder"]}")
        shutil.move(results_file, run_informations["archives_results_folder"])

    with open(
        osj(run_informations["archives_results_folder"], "VANNOTCopyComplete.txt"),
        mode="a",
    ):
        pass

    shutil.rmtree(run_informations["tmp_analysis_folder"])
    log.info("Deleted temporary analysis folder")


if __name__ == "__main__":
    pass
