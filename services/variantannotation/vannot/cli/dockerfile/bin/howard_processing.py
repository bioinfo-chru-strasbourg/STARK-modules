# -*- coding: utf-8 -*-
import glob
from os.path import join as osj
import os
import subprocess
import logging as log
import shutil
import json
import re
from tools import bcftools, bgzip, tabix, intersectBed

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

    log.info("Copying vcf files to temporary folder")
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

    for vcf_file in vcf_file_to_analyse:
        cleaned_vcf = cleaning_annotations(vcf_file, run_informations)
        howard_proc(run_informations, cleaned_vcf)
        os.remove(cleaned_vcf)

def cleaning_annotations(vcf_file, run_informations):
    module_config = osj(os.environ["DOCKER_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        annotations_to_keep = data["keep_vcf_info"][run_informations["run_platform_application"]]
        
    log.info("Cleaning INFO column in the provided vcfs")
    log.info(f"Kept informations : {", ".join(annotations_to_keep)}")
    
    actual_info_fields = subprocess.run(["zgrep", "##INFO", vcf_file], capture_output=True, text=True)
    actual_info_fields = actual_info_fields.stdout.strip().split("##")
    if not vcf_file.endswith(".vcf.gz"):
        subprocess.call([bgzip, vcf_file])
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

    cmd = [bcftools, "annotate", "-x"]
    cmd.append(info_to_keep)
    cmd.append(vcf_file)

    with open(cleaned_vcf, "w") as output:
        subprocess.call(cmd, stdout=output, universal_newlines=True)
    os.remove(vcf_file)
        
    subprocess.call([bgzip, cleaned_vcf], universal_newlines=True)
    renamed_clean = osj(os.path.dirname(cleaned_vcf), os.path.basename(cleaned_vcf).replace("cleaned_", "") + ".gz")
    cleaned_vcf = cleaned_vcf + ".gz"
    os.rename(cleaned_vcf, renamed_clean)

    return renamed_clean

def howard_proc(run_informations, vcf_file):
    log.info(f"Launching HOWARD analysis for {vcf_file}")

    if run_informations["output_format"] != None:
        output_file = osj(
            run_informations["tmp_analysis_folder"], f"VANNOT_{os.path.basename(vcf_file).split(".")[0]}.design.{run_informations["output_format"]}"
        )
    else:
        output_file = osj(
            run_informations["tmp_analysis_folder"], f"VANNOT_{os.path.basename(vcf_file).split(".")[0]}.design.{".".join(os.path.basename(vcf_file).split(".")[1:])}"
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
        raise ValueError(configfile)

    container_name = f"VANNOT_annotate_{run_informations['run_name']}_{os.path.basename(vcf_file).split('.')[0]}"
    launch_annotate_arguments = ["annotation", "--input", vcf_file, "--output", output_file, "--param", configfile, "--assembly", run_informations["assembly"]]

    log.info("Annotating input files with HOWARD")
    
    howard_launcher.launch(container_name, launch_annotate_arguments)
    convert_to_final_tsv(run_informations, output_file, "")

def merge_vcf_files(run_informations):
    log.info("Merging all vcfs into one for CuteVariant analysis")
    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))
    output_file = osj(run_informations["tmp_analysis_folder"], f"merged_VANNOT_{run_informations["run_name"]}.design.vcf.gz")
    
    for vcf_file in vcf_files:
        subprocess.call([tabix, vcf_file], universal_newlines=True)

    cmd = [bcftools, "merge", "-o", output_file, "-O", "z"]
    for vcf_file in vcf_files:
        cmd.append(vcf_file)

    subprocess.call(cmd, universal_newlines=True)
    

def convert_to_final_tsv(run_informations, input_file, panel_name):
    log.info("Converting output file into readable tsv")
    container_name = f"VANNOT_convert_{run_informations['run_name']}_{os.path.basename(input_file).split('.')[0]}"
    if panel_name != "":
        output_file = osj(run_informations["tmp_analysis_folder"], f"{os.path.basename(input_file).split(".")[0]}.{panel_name}.tsv")
    else:
        output_file = osj(run_informations["tmp_analysis_folder"], f"{os.path.basename(input_file).split(".")[0]}.design.tsv")
    
    if not input_file.startswith("merged"):
        launch_convert_arguments = ["convert", "--input", input_file, "--output", output_file, "--explode_infos"]
        howard_launcher.launch(container_name, launch_convert_arguments)

def cleaner(run_informations):
    log.info("Moving results from temporary folder")
    results_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*tsv")) + glob.glob(osj(run_informations["tmp_analysis_folder"], "*vcf.gz"))

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
    
def panel_filtering(run_informations):
    tmp_vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*VANNOT_*vcf.gz"))
    panels = run_informations["run_panels"]
    
    for panel in panels:
        is_empty = True
        subprocess.call(["rsync", "-rvt", panel, run_informations["tmp_analysis_folder"]], universal_newlines=True)
        panel = os.path.join(run_informations["tmp_analysis_folder"], os.path.basename(panel))
        if os.path.basename(panel).split(".")[2] != "vcf":
            panel_name = os.path.basename(panel).split(".")[2]
        else:
            panel_name = "_".join(os.path.basename(panel).split(".")[1].split("_")[1:])
        for tmp_vcf_file in tmp_vcf_files:
            sample_name = tmp_vcf_file.split(".")[0]
            filtered_vcf = osj(run_informations["tmp_analysis_folder"], sample_name + "." + panel_name + ".vcf")
            command_list = [intersectBed, "-a", tmp_vcf_file, "-b", panel,"-header"]
            log.info(" ".join(command_list))
            with open(filtered_vcf, "a") as f : 
                subprocess.call(command_list, stdout=f, stderr=subprocess.STDOUT, universal_newlines=True)
            with open(filtered_vcf, "r") as readfile :
                lines = readfile.readlines()
                for line in lines:
                    if line.startswith("chr"):
                        is_empty = False
            
            if is_empty is True:
                log.info("Filtered VCF is empty, please check your .gene file, conversion to tsv is aborted")
                
            subprocess.call([bgzip, filtered_vcf])
            filtered_vcf = filtered_vcf + ".gz"
            if is_empty is False:
                convert_to_final_tsv(run_informations, filtered_vcf, panel_name)

if __name__ == "__main__":
    pass
