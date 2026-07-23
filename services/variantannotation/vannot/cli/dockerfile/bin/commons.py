import os
import logging as log
from os.path import join as osj
from threading import local
import sys
import time
import json
import glob
import shutil
import gzip
import subprocess

def get_threads(threads_type):
    module_config = osj(os.environ["HOST_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    if not os.path.isfile(module_config):
        log.error(f"{module_config} do not exist, primordial file, check its existence")
        raise ValueError(module_config)
      
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        threads = data[threads_type]

    return threads

def create_listbypipeline_vcf(run_informations):
    run_repository = run_informations["run_repository"]
    tmp_analysis_folder = run_informations["tmp_analysis_folder"]
    sample_list = [i.split("/")[-1].split(".")[0] for i in glob.glob(osj(tmp_analysis_folder, "*")) if i.endswith(".vcf.gz")]
    if not os.path.isdir(osj(tmp_analysis_folder, "fullvcf")):
        os.makedirs(osj(tmp_analysis_folder, "fullvcf"))
    for i in sample_list:
        full_vcf = osj(run_repository, i, f"{i}.full.Design.vcf.gz")
        shutil.copy(full_vcf, osj(tmp_analysis_folder, "fullvcf", f"{i}.full.Design.vcf.gz"))

    for vcf_file in glob.glob(osj(tmp_analysis_folder, "fullvcf", "*.full.Design.vcf.gz")):
        annotation_file = vcf_file.replace(".vcf.gz", ".annotation.tsv")
        header_file = f"{annotation_file}.hdr"

        with gzip.open(vcf_file, "rt") as read_vcf, open(annotation_file, "w") as write_annot:
            write_annot.write("#CHROM\tPOS\tREF\tALT\tLISTBYPIPELINE\tGQ_list\n")
            callers = []
            for line in read_vcf:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header_fields = line.strip().split("\t")
                    callers = header_fields[9:]
                    continue

                fields = line.strip().split("\t")
                chrom, pos, _, ref, alt = fields[0:5]
                format_fields = fields[8].split(":")
                gt_index = format_fields.index("GT") if "GT" in format_fields else None
                gq_index = format_fields.index("GQ") if "GQ" in format_fields else None
                samples_data = fields[9:]

                called_callers = []
                gq_values = []
                for caller, sample_data in zip(callers, samples_data):
                    sample_values = sample_data.split(":")
                    if gt_index is None or gt_index >= len(sample_values):
                        gq_values.append(".")
                        continue
                    gt = sample_values[gt_index].replace("|", "/")
                    if gt in (".", "./.", "./0", "0/."):
                        gq_values.append(".")
                        continue

                    called_callers.append(caller)
                    if gq_index is not None and gq_index < len(sample_values):
                        gq_values.append(sample_values[gq_index])
                    else:
                        gq_values.append(".")

                listbypipeline = ",".join(called_callers) if called_callers else "."
                gq_list = ",".join(gq_values) if gq_values else "."
                write_annot.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{listbypipeline}\t{gq_list}\n")

        with open(header_file, "w") as write_hdr:
            write_hdr.write(
                '##FORMAT=<ID=ListByPipeline,Number=.,Type=String,Description="List of callers/pipelines which found the variant">\n'
                '##FORMAT=<ID=GQ_list,Number=.,Type=String,Description="List of GQ values for each caller which found the variant">\n'
            )
        os.remove(vcf_file)

    for sample in sample_list:
        annot_tsv = osj(tmp_analysis_folder, "fullvcf", f"{sample}.full.Design.annotation.tsv")
        original_vcf = osj(tmp_analysis_folder, f"{sample}.vcf.gz")
        header_file = osj(tmp_analysis_folder, "fullvcf", f"{sample}.full.Design.annotation.tsv.hdr")

        if not os.path.isfile(annot_tsv) or not os.path.isfile(original_vcf):
            continue

        annot_tsv_gz = f"{annot_tsv}.gz"
        subprocess.run(
            f"bgzip -f -c {annot_tsv} > {annot_tsv_gz} && tabix -f -s1 -b2 -e2 {annot_tsv_gz}",
            shell=True,
            check=True,
        )

        tmp_vcf = f"{original_vcf}.tmp.vcf.gz"
        cmd=[
            "bcftools", "annotate",
            "-a", annot_tsv_gz,
            "-h", header_file,
            "-c", "CHROM,POS,REF,ALT,FORMAT/LISTBYPIPELINE,FORMAT/GQ_list",
            "-Oz", "-o", tmp_vcf,
            original_vcf,
        ]
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)

        shutil.move(tmp_vcf, original_vcf)
        subprocess.run(["tabix", "-f", "-p", "vcf", original_vcf], check=True)

    shutil.rmtree(osj(tmp_analysis_folder, "fullvcf"), ignore_errors=True)
    tbi_files = glob.glob(osj(tmp_analysis_folder, "*.tbi"))
    for tbi_file in tbi_files:
        os.remove(tbi_file)


def get_memory(memory_type):
    module_config = osj(os.environ["HOST_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    if not os.path.isfile(module_config):
        log.error(f"{module_config} do not exist, primordial file, check its existence")
        raise ValueError(module_config)
      
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        memory = data[memory_type]

    return memory

def get_default_pattern(value):
    if value == "default":
        default_pattern = "*/STARK/*.reports/*.final.vcf.gz"
    if isinstance(value, dict):
        if value["run_repository"] is not None:
            default_pattern = "*/STARK/*.reports/*.final.vcf.gz"
        elif value["run_archives"] is not None:
            default_pattern = "*/*.final.vcf.gz"
    return default_pattern


def set_logger_info():
    mylog = log.getLogger()
    log_file = mylog.handlers[0].baseFilename
    info_file = log_file.replace(".log", ".info")

    if not os.path.isfile(info_file):
        python_command = " ".join(sys.argv)
        python_version = f"{sys.version.split(' ')[0].split('.')[0]}.{sys.version.split(' ')[0].split('.')[1]}"
        command = f"python{python_version} {python_command}"
        global time_seconds_start
        time_seconds_start = time.time() + 7200
        local_time = time.localtime(time_seconds_start)
        actual_time = time.strftime("%a %b %d %H:%M:%S %Y", local_time)
        start = actual_time

        logging = f"""Command: {command}
Start time: {start}\n"""

        with open(info_file, "a") as write_file:
            write_file.write(logging)

    else:
        time_seconds_end = time.time() + 7200
        local_time = time.localtime(time_seconds_end)
        actual_time = time.strftime("%a %b %d %H:%M:%S %Y", local_time)
        end = actual_time
        delta = time_seconds_end - time_seconds_start

        logging = f"""End time: {end}
Time run: {delta}"""

        with open(info_file, "a") as write_file:
            write_file.write(logging)


def logger_header(log_file):
    logging = f"""#########################
#         vAnnot        #
# Author: Mateusz Rauch #
#########################

####################
#   Release: {os.environ["DOCKER_SERVICE_CLI_RELEASE"]} #
####################

"""
    with open(log_file, "a") as write_file:
        write_file.write(logging)


def set_log_level(args):
    verbosity = args.verbosity
    time_seconds = time.time() + 7200
    local_time = time.localtime(time_seconds)
    actual_time = time.strftime("%Y%m%d_%H%M%S", local_time)

    if "run" in args:
        run = args.run
        if run.endswith("/"):
            run = run[:-1]
        mode = args.launchmode
        run_name = run.split("/")[-1]
        run_application = run.split("/")[-2]
        run_platform = run.split("/")[-3]
        log_file = f"{actual_time}_{run_platform}_{run_application}_{run_name}_{mode}.log"

    elif "folder" in args:
        folder = args.folder
        if folder.endswith("/"):
            folder = folder[:-1]
        folder_name = folder.split("/")[-1]
        log_file = f"{actual_time}_{folder_name}.log"
        
    elif "dejavu" in args:
        if args.run_dejavu is not None:
            dejavu = args.run_dejavu
        else:
            dejavu = args.dejavu
        if dejavu.endswith("/"):
            dejavu = dejavu[:-1]
        dejavu_name = dejavu.split("/")[-1]
        log_file = f"{actual_time}_dejavuonly_{dejavu_name}.log"

    log_file = osj(
        os.environ["HOST_SERVICES"],
        "logs",
        log_file,
    )
    if not os.path.exists(os.path.dirname(log_file)):
        os.makedirs(os.path.dirname(log_file))

    logger_header(log_file)

    configs = {
        "debug": log.DEBUG,
        "info": log.INFO,
        "warning": log.WARNING,
        "error": log.ERROR,
        "critical": log.CRITICAL,
    }
    if verbosity not in configs.keys():
        raise ValueError(
            "Unknown verbosity level:"
            + verbosity
            + "\nPlease use any in:"
            + configs.keys()
        )

    log.basicConfig(
        filename=log_file,
        force=True,
        filemode="a",
        format="vAnnot %(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )


if __name__ == "__main__":
    pass
