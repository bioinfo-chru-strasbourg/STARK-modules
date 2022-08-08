# -*- coding: utf-8 -*-
import glob
from os.path import join as osj
import os
import commons
import subprocess
import logging as log
import re
import shutil

def initialisation(run_informations):
    archives_vcf_files = glob.glob(osj(run_informations["run_processing_folder"], "VCF", "*", "*.final.vcf*"))
    archives_run_name = glob.glob(osj(run_informations["run_processing_folder"], "VCF", "*"))
    other_vcfs = glob.glob(osj(run_informations["run_processing_folder"], "*.final.vcf*")) + glob.glob(osj(run_informations["run_processing_folder"], "VCF", "*.final.vcf*"))

    deleted_samples_file = osj(os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG"], "variantannotation", "varank", "configfiles", "deleted_samples_config.tsv")

    with open(deleted_samples_file, "r") as read_file:
        next(read_file)
        for line in read_file.readlines():
            line = line.split("\t")
            if line[4] == os.path.basename(run_informations["run_processing_folder"]):
                vcf_file_to_delete = osj(run_informations["run_processing_folder"], "VCF", line[1], line[0])
                if os.path.isfile(vcf_file_to_delete):
                    os.remove(vcf_file_to_delete)

    if len(archives_vcf_files) != 0:
        archives_vcf_files = sorted(archives_vcf_files)
        for run_name in archives_run_name:
            vcf_file_list = glob.glob(osj(run_name, "*.final.vcf*"))
            for processed_vcf_file in vcf_file_list:
                subprocess.run(["rsync","-rp", processed_vcf_file, run_informations["run_processing_folder"]])
    
    if len(other_vcfs) != 0:
        if not os.path.isdir(osj(run_informations["run_processing_folder"], "VCF", "DEFAULT")):
            os.makedirs(osj(run_informations["run_processing_folder"], "VCF", "DEFAULT", ""), 0o775)
        archives_default_vcf_folder = osj(run_informations["run_processing_folder"], "VCF", "DEFAULT", "")

        for archives_other_processing_vcf_files in other_vcfs:
            shutil.move(archives_other_processing_vcf_files, archives_default_vcf_folder)

        archives_default_vcf_files = glob.glob(osj(archives_default_vcf_folder, "*.final.vcf*"))

        for archives_default_processed_vcf_files in archives_default_vcf_files:
            subprocess.run(["rsync","-rp",archives_default_processed_vcf_files, run_informations["run_processing_folder"]])

    if os.path.isdir(run_informations["run_processing_folder_tsv"]):
        log.info(f"Cleaning old files in {run_informations['run_processing_folder_tsv']}")
        with open(osj(run_informations["run_processing_folder"], "TSVDeleting.txt"), mode="a"):
            pass
        shutil.rmtree(run_informations["run_processing_folder_tsv"])
        os.mkdir(run_informations["run_processing_folder_tsv"])
        os.remove(osj(run_informations["run_processing_folder"], "TSVDeleting.txt"))
    else:
        os.mkdir(run_informations["run_processing_folder_tsv"])
    
    