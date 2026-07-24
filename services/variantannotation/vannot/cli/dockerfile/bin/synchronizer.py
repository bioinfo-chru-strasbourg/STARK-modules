import glob
import logging as log
from os.path import join as osj
import os
import commons
import subprocess
import re
import json
import shutil


def find_samplesheet(run_informations, archives):
    if archives is True :
        samples = glob.glob(osj(run_informations["run_archives"], "*", ""))
        for sample in samples:
            samplesheet = osj(sample, f"{os.path.basename(sample.strip('/'))}.SampleSheet.csv")
            is_ss = os.path.isfile(samplesheet)
            if is_ss is True:
                break
    else: 
        samples = glob.glob(osj(run_informations["run_repository"], "*", ""))
        for sample in samples:
            samplesheet = osj(sample, f"{os.path.basename(sample.strip('/'))}.SampleSheet.csv")
            is_ss = os.path.isfile(samplesheet)
            if is_ss is True:
                break
    return samplesheet


def find_tag(samplesheet, word):
    tagged_samples = []
    with open(samplesheet, "r") as read_file:
        for line in read_file:
            line = line.strip()
            if word in line:
                line = line.split(",")
                tagged_samples.append(line[0])

    return tagged_samples


def vcf_synchronizer(run_informations):
    if run_informations["type"] == "run_dejavu" and run_informations["run_repository"] == "":
        archives = True
        run_path = run_informations["run_archives"]
    else:
        archives = False
        run_path = run_informations["run_repository"]

    pattern = run_informations["vcf_pattern"]

    module_config = osj(
        os.environ["HOST_MODULE_CONFIG"],
        f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
    )
    ignored_samples = []
    if run_informations["onco"] == False:
        with open(module_config, "r") as read_file:
            data = json.load(read_file)
            ignored_samples = data["ignored_samples"]
        samplesheet = find_samplesheet(run_informations, archives)
        control_samples = find_tag(samplesheet, "CQI#", archives)
        platform_application = run_informations["run_platform_application"]
        pool_tag = f"APP#{platform_application}#POOL"
        pool_samples = find_tag(samplesheet, pool_tag, archives)

        ignored_samples = ignored_samples + control_samples + pool_samples
    else:
        ignored_samples = []
    log.info(
        "Ignoring following sample patterns for the analysis and dejavu generation : "
        + ", ".join(ignored_samples)
    )

    if os.path.isdir(run_informations["archives_run_folder"]):
        shutil.rmtree(run_informations["archives_run_folder"])
    elif os.path.isfile(run_informations["archives_run_folder"]):
        os.remove(run_informations["archives_run_folder"])

    os.makedirs(run_informations["archives_run_folder"])
    os.chmod(run_informations["archives_run_folder"], 0o777)

    kept_vcf = []
    treated_samples = []

    for element in reversed(pattern):
        vcf_files = glob.glob(osj(run_path, element))
        for vcf_file in vcf_files:
            sample = vcf_file.split("/")[-1].split(".")[0]
            if archives is False:
                stark_vcf = (
                    run_path
                    + "/"
                    + sample
                    + "\\/STARK\\/"
                    + sample
                    + ".reports\\/"
                    + sample
                    + ".\\d{8}-\\d{6}.final.vcf.gz"
                )
            elif archives is True:
                stark_vcf = (
                    run_path
                    + "/"
                    + sample
                    + "/"
                    + sample
                    + ".final.vcf.gz"
                )

            if (
                not re.match(stark_vcf, vcf_file)
                and sample not in treated_samples
                and sample not in ignored_samples
            ):
                kept_vcf.append(vcf_file)
            elif re.match(stark_vcf, vcf_file) and sample not in treated_samples and sample not in ignored_samples and archives is True:
                kept_vcf.append(stark_vcf)

            if element != commons.get_default_pattern(run_informations):
                log.info(f"Keeping the sample vcf with {element} pattern")
                treated_samples.append(sample)

    for ignored_sample in ignored_samples:
        for sample_vcf in kept_vcf:
            if ignored_sample in sample_vcf:
                kept_vcf.remove(sample_vcf)

    for vcf_file in kept_vcf:
        output = subprocess.check_output(f'zgrep -v \"#\" {vcf_file} | wc -l', shell=True, text=True)
        if int(output) == 0:
            kept_vcf.remove(vcf_file)

    for vcf_file in kept_vcf:
        vcf_file_output = os.path.basename(vcf_file).split(".")[0] + ".vcf.gz"
        log.info(f"Synchronizing {vcf_file}")
        subprocess.run(
            [
                "rsync",
                "-rp",
                vcf_file,
                osj(run_informations["archives_run_folder"], vcf_file_output),
            ]
        )
