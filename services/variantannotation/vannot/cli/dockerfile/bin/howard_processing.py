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
        if os.path.basename(vcf_file).split(".")[1] == "POOL":
            fixed_vcf_file = info_to_format_script(vcf_file, run_informations)
        # exomiser_annotation()
        else:
            fixed_vcf_file = vcf_file
        cleaning_annotations(fixed_vcf_file, run_informations)

        # unmerge_vcf()
        
    # os.remove(cleaned_merged_vcf)
    
    
    
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
        subprocess.call(["bgzip", vcf_file])
        vcf_file = vcf_file + ".gz"

    cleaned_vcf = osj(os.path.dirname(vcf_file), "cleaned_" + os.path.basename(vcf_file)[:-3])
    info_to_keep = []
    for i in annotations_to_keep:
        for j in actual_info_fields:
            if re.search(i, j):
                info_to_keep.append("INFO/" + j.split(",")[0].split("=")[-1])
                
    if len(info_to_keep) == 0:
        log.info(f"No annotations to keep were found in {os.path.basename(vcf_file)}, deleting all annotations")
        info_to_keep = "INFO"
    else:
        log.info(f"Keeping following annotations: {' '.join(info_to_keep)} for sample {os.path.basename(vcf_file)}")
        info_to_keep = "^" + ",".join(info_to_keep)

    cmd = ["bcftools", "annotate", "-x"]
    cmd.append(info_to_keep)
    cmd.append(vcf_file)

    with open(cleaned_vcf, "w") as output:
        subprocess.call(cmd, stdout=output, universal_newlines=True)
    os.remove(vcf_file)
        
    subprocess.call(["bgzip", cleaned_vcf], universal_newlines=True)
    renamed_clean = osj(os.path.dirname(cleaned_vcf), os.path.basename(cleaned_vcf).replace("cleaned_", "") + ".gz")
    cleaned_vcf = cleaned_vcf + ".gz"
    os.rename(cleaned_vcf, renamed_clean)

    return renamed_clean

def merge_vcf(run_informations):
    vcf_file_to_merge = glob.glob(osj(run_informations["tmp_analysis_folder"], "*vcf*"))
    for i in vcf_file_to_merge:
        subprocess.call(["tabix", i], universal_newlines=True)
    
    output_merged = osj(run_informations["tmp_analysis_folder"], f"merged_{run_informations["run_name"]}.vcf.gz")
    cmd = ["bcftools", "merge"] + vcf_file_to_merge
    cmd_args = ["-m", "none", "-O", "z", "-o", output_merged]
    cmd = cmd + cmd_args
    log.debug(" ".join(cmd))
    subprocess.call(cmd, universal_newlines=True)
    for i in vcf_file_to_merge:
        os.remove(i)
        os.remove(i + ".tbi")
    return(output_merged)

def unmerge_vcf(run_informations, input):
    vcf_file_to_unmerge = input
    print(vcf_file_to_unmerge)


def info_to_format_script(vcf_file, run_informations):
    log.info(f"Moving POOL INFO sample specific columns to FORMAT column for {os.path.basename(vcf_file)}")
    module_config = osj(os.environ["DOCKER_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    output_file = osj(os.path.dirname(vcf_file), f"fixed_{os.path.basename(vcf_file)[:-3]}")
    sample = os.path.basename(vcf_file).split(".")[0]
    tmp_annot = osj(os.path.dirname(vcf_file), "annot.txt.tmp")
    tmp_annot_fixed = osj(os.path.dirname(vcf_file), "annot.fixed.txt.tmp")
    tmp_hdr = osj(os.path.dirname(vcf_file), "hdr.txt.tmp")
        
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        info_to_format_columns = data["info_to_format"][run_informations["run_platform_application"]]
    info_to_format_columns_query = "\\t%".join(info_to_format_columns)
    info_to_format_columns_query = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%" + info_to_format_columns_query + "\n"

    cmd = ["bcftools", "query", "-f", info_to_format_columns_query, vcf_file]
    
    with open(tmp_annot, "a") as writefile:
        subprocess.call(cmd, universal_newlines=True, stdout=writefile)

    with open(tmp_annot, "r") as readfile:
        with open(tmp_annot_fixed, "a") as writefile:
            lines = readfile.readlines()
            for line in lines:
                line = line.strip()
                writefile.write(line.replace(":", ",") + "\n")

    subprocess.call(["bgzip", tmp_annot_fixed], universal_newlines=True)    
    tmp_annot_fixed = tmp_annot_fixed + ".gz"
    cmd = ["tabix", "-s1", "-b2", "-e2", tmp_annot_fixed]
    subprocess.call(cmd, universal_newlines=True)

    vcf_file_gunzip = vcf_file[:-3]
    subprocess.call(["gunzip", vcf_file])
    with open(vcf_file_gunzip, "r") as readfile:
        with open(tmp_hdr, "a") as writefile:
            lines = readfile.readlines()
            for line in lines:
                for column in info_to_format_columns:
                    if line.startswith("##INFO=<ID=" + column):
                        writefile.write(line.replace("INFO", "FORMAT"))
    subprocess.call(["bgzip", vcf_file_gunzip], universal_newlines=True)

    info_to_format_columns_annotate = ",FORMAT/".join(info_to_format_columns)
    info_to_format_columns_annotate = "CHROM,POS,REF,ALT,FORMAT/" + info_to_format_columns_annotate

    cmd = ["bcftools", "annotate", "-s", sample, "-a", tmp_annot_fixed, "-h", tmp_hdr, "-c", info_to_format_columns_annotate, vcf_file]
    log.debug(" ".join(cmd))
    with open(output_file, "a") as writefile:
        subprocess.call(cmd, universal_newlines=True, stdout=writefile)

    subprocess.call(["bgzip", output_file], universal_newlines=True)
    output_file = output_file + ".gz"
    
    os.remove(tmp_annot)
    os.remove(tmp_annot_fixed)
    os.remove(tmp_annot_fixed + ".tbi")
    os.remove(tmp_hdr)
    os.remove(vcf_file)
    
    return(output_file)

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

def merge_vcf_files(run_informations):
    log.info("Merging all vcfs into one for CuteVariant analysis")
    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))
    output_file = osj(run_informations["tmp_analysis_folder"], f"merged_VANNOT_{run_informations["run_name"]}.design.vcf.gz")
    
    for vcf_file in vcf_files:
        subprocess.call(["tabix", vcf_file], universal_newlines=True)

    cmd = ["bcftools", "merge", "-o", output_file, "-O", "z"]
    for vcf_file in vcf_files:
        cmd.append(vcf_file)

    subprocess.call(cmd, universal_newlines=True)
    

def convert_to_final_tsv(run_informations, input_file, panel_name):
    log.info("Converting output file into readable tsv")
    container_name = f"VANNOT_convert_{run_informations['run_name']}_{os.path.basename(input_file).split('.')[0]}"
    if panel_name != "":
        output_file = osj(run_informations["tmp_analysis_folder"], f"{os.path.basename(input_file).split(".")[0]}.panel.{panel_name}.tsv")
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
        subprocess.call(["rsync", "-rvt", panel, run_informations["tmp_analysis_folder"]], universal_newlines=True)
        panel = os.path.join(run_informations["tmp_analysis_folder"], os.path.basename(panel))
        if os.path.basename(panel).split(".")[3] != "genes":
            panel_name = os.path.basename(panel).split(".")[1]
        else:
            panel_name = "_".join(os.path.basename(panel).split(".")[1].split("_")[1:])
        for tmp_vcf_file in tmp_vcf_files:
            sample_name = tmp_vcf_file.split(".")[0]
            filtered_vcf = osj(run_informations["tmp_analysis_folder"], sample_name + ".panel." + panel_name + ".vcf")
            command_list = ["intersectBed", "-a", tmp_vcf_file, "-b", panel,"-header"]
            log.info(" ".join(command_list))
            with open(filtered_vcf, "a") as f : 
                subprocess.call(command_list, stdout=f, stderr=subprocess.STDOUT, universal_newlines=True)
                
            subprocess.call(["bgzip", filtered_vcf])
            filtered_vcf = filtered_vcf + ".gz"

if __name__ == "__main__":
    pass
