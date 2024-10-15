import glob
import logging as log
from os.path import join as osj
import os
import commons
import subprocess
import re
import json
import shutil

def find_samplesheet(run_informations):
    samples = glob.glob(osj(run_informations["run_repository"], "*", ""))
    samplesheet = glob.glob(osj(samples[0], "STARK", "*.SampleSheet.csv"))
    return samplesheet[0]
    
def find_controls(samplesheet):
    word = "CQI#"
    controls_samples = []
    with open(samplesheet, "r") as read_file:
        for line in read_file:
            line = line.strip()
            if word in line:
                line = line.split(",")
                controls_samples.append(line[0])
            
    return controls_samples

def vcf_synchronizer(run_informations):
    run_repository = run_informations["run_repository"]
    pattern = run_informations["vcf_pattern"]
    
    module_config = osj(os.environ["DOCKER_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        ignored_samples = data["ignored_samples"]
    samplesheet = find_samplesheet(run_informations)
    control_samples = find_controls(samplesheet)
    ignored_samples = ignored_samples + control_samples
    log.info("Ignoring following sample patterns for the analysis and dejavu generation : " + ", ".join(ignored_samples))

    if os.path.isdir(run_informations["archives_run_folder"]):
        shutil.rmtree(run_informations["archives_run_folder"])
    elif os.path.isfile(run_informations["archives_run_folder"]):
        os.remove(run_informations["archives_run_folder"])
    elif not os.path.isdir(run_informations["archives_run_folder"]):
        os.makedirs(run_informations["archives_run_folder"])
        os.chmod(run_informations["archives_run_folder"], 0o777)      

    kept_vcf = []
    treated_samples = []
    for element in reversed(pattern):
        vcf_files = glob.glob(osj(run_repository, element))
        for vcf_file in vcf_files:
            sample = vcf_file.split("/")[-1].split(".")[0]
            dated_stark_vcf = run_repository + "/" + sample + "\\/STARK\\/" + sample + ".reports\\/" + sample + ".\\d{8}-\\d{6}.final.vcf.gz"
            if not re.match(dated_stark_vcf, vcf_file) and sample not in treated_samples and sample not in ignored_samples:
                kept_vcf.append(vcf_file)

            if element != commons.get_default_pattern():
                log.info(f"Keeping the sample vcf with {pattern} pattern")
                treated_samples.append(sample)

    for ignored_sample in ignored_samples:
            for sample_vcf in kept_vcf:
                if ignored_sample in sample_vcf:
                    kept_vcf.remove(sample_vcf)

    for vcf_file in kept_vcf:
        log.info(f"Synchronizing {vcf_file}")
        subprocess.run(["rsync", "-rp", vcf_file, run_informations["archives_run_folder"]])
