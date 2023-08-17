# -*- coding: utf-8 -*-
import glob
from os.path import join as osj
import os
import subprocess
import logging as log
import shutil
import json
import multiprocessing


def root_vcf_initialisation(run_informations):
    other_vcfs = glob.glob(
        osj(run_informations["run_processing_folder"], "*.final.vcf*")
    ) + glob.glob(osj(run_informations["run_processing_folder"], "*", "*.final.vcf*"))

    if len(other_vcfs) != 0:
        if not os.path.isdir(
            osj(run_informations["run_processing_folder"], "VCF", "DEFAULT")
        ):
            os.makedirs(
                osj(run_informations["run_processing_folder"], "VCF", "DEFAULT", ""),
                0o775,
            )
        archives_default_vcf_folder = osj(
            run_informations["run_processing_folder"], "VCF", "DEFAULT", ""
        )

        for archives_other_processing_vcf_files in other_vcfs:
            subprocess.run(
                [
                    "rsync",
                    "-rp",
                    archives_other_processing_vcf_files,
                    archives_default_vcf_folder,
                ]
            )
            os.remove(archives_other_processing_vcf_files)

    vcf_files = glob.glob(
        osj(run_informations["run_processing_folder"], "VCF", "*", "*.final.vcf*")
    )

    tmp_folders = glob.glob("/tmp_*")
    for tmp_folder in tmp_folders:
        shutil.rmtree(tmp_folder)

    os.mkdir(run_informations["analysis_folder"])

    for vcf_file in vcf_files:
        subprocess.run(
            [
                "rsync",
                "-rp",
                vcf_file,
                run_informations["analysis_folder"],
            ]
        )

    if os.path.isdir(osj(run_informations["run_processing_folder"], "Alamut")):
        subprocess.run(
            [
                "rsync",
                "-rp",
                osj(run_informations["run_processing_folder"], "Alamut"),
                run_informations["analysis_folder"],
            ]
        )


def initialisation(run_informations):
    archives_vcf_files = glob.glob(
        osj(run_informations["run_processing_folder"], "VCF", "*", "*.final.vcf*")
    )
    archives_run_name = glob.glob(
        osj(run_informations["run_processing_folder"], "VCF", "*")
    )

    deleted_samples_file = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG"],
        "variantannotation",
        "varank",
        "configfiles",
        "deleted_samples_config.tsv",
    )

    with open(deleted_samples_file, "r") as read_file:
        next(read_file)
        for line in read_file.readlines():
            line = line.split("\t")
            if line[4] == os.path.basename(run_informations["run_processing_folder"]):
                vcf_file_to_delete = osj(
                    run_informations["run_processing_folder"], "VCF", line[1], line[0]
                )
                if os.path.isfile(vcf_file_to_delete):
                    os.remove(vcf_file_to_delete)

    if len(archives_vcf_files) != 0:
        archives_vcf_files = sorted(archives_vcf_files)
        for run_name in archives_run_name:
            vcf_file_list = glob.glob(osj(run_name, "*.final.vcf*"))
            for processed_vcf_file in vcf_file_list:
                subprocess.run(
                    [
                        "rsync",
                        "-rp",
                        processed_vcf_file,
                        run_informations["analysis_folder"],
                    ]
                )

    if os.path.isdir(osj(run_informations["run_processing_folder"], "Alamut")):
        subprocess.run(
            [
                "rsync",
                "-rp",
                osj(run_informations["run_processing_folder"], "Alamut"),
                run_informations["analysis_folder"],
            ]
        )

    if os.path.isdir(run_informations["run_processing_folder_tsv"]):
        log.info(
            f"Cleaning old files in {run_informations['run_processing_folder_tsv']}"
        )
        with open(
            osj(run_informations["run_processing_folder"], "TSVDeleting.txt"), mode="a"
        ):
            pass
        for tsv_file in glob.glob(
            osj(run_informations["run_processing_folder_tsv"], "*")
        ):
            os.remove(tsv_file)
        os.remove(osj(run_informations["run_processing_folder"], "TSVDeleting.txt"))
    else:
        os.mkdir(run_informations["run_processing_folder_tsv"])


def configfile_manager(run_informations):
    configfile_shelter = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG"],
        "variantannotation",
        "varank",
        "configfiles",
    )
    used_configfile = osj(run_informations["run_processing_folder"], "configfile")

    if run_informations["run_platform_application"] == "default":
        configfile = osj(configfile_shelter, "configfile.default")
    else:
        configfile = osj(
            configfile_shelter,
            f"configfile.{run_informations['run_platform_application']}",
        )

    if (
        os.path.isfile(configfile)
        and not os.path.isfile(used_configfile)
        and not run_informations["run_platform_application"] == "default"
    ):
        log.info(
            f"Configfile {configfile} in {run_informations['run_processing_folder']} was synced from the configfile shelter"
        )
        subprocess.run(["rsync", "-rp", configfile, used_configfile])

    elif (
        os.path.isfile(configfile)
        and os.path.isfile(used_configfile)
        and not run_informations["run_platform_application"] == "default"
    ):
        log.info(
            f"Configfile {configfile} in {run_informations['run_processing_folder']} was updated from the configfile shelter"
        )
        subprocess.run(["rsync", "-rp", configfile, used_configfile])

    elif (
        os.path.isfile(used_configfile)
        and run_informations["run_platform_application"] == "default"
    ):
        log.info(
            f"Using existing configfile in {run_informations['run_processing_folder']}"
        )

    elif (
        not os.path.isfile(used_configfile)
        and run_informations["run_platform_application"] == "default"
    ):
        log.info(f"Using default configfile in configfile shelter")
        subprocess.run(["rsync", "-rp", configfile, used_configfile])

    subprocess.run(
        ["rsync", "-rp", used_configfile, run_informations["analysis_folder"]]
    )


def varank_launcher(run_informations):
    logfile = osj(run_informations["analysis_folder"], "VaRank.log")
    varank_bin = osj(os.environ["VARANK"], "bin", "VaRank")

    varank_json = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG"],
        "variantannotation",
        "varank",
        "varank_config.json",
    )

    with open(varank_json, "r") as read_file:
        data = json.load(read_file)
        threads = data.get("threads")
        if threads == "all":
            used_threads = multiprocessing.cpu_count() - 1
        else:
            used_threads = threads

    log.info(f"Generating TSV files, you can follow the progress in {logfile}")
    with open(logfile, "w") as f:
        subprocess.call(
            [
                varank_bin,
                "-vcfdir",
                run_informations["analysis_folder"],
                "-alamutHumanDB",
                "hg19",
                "-SamVa",
                '"yes"',
                "-AlamutProcesses",
                str(used_threads),
            ],
            stdout=f,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
        )
    os.chmod(logfile, 0o777)


def cleaner(run_informations):
    tsv_files = glob.glob(osj(run_informations["analysis_folder"], "*allVariants*.tsv"))
    statistics_files = glob.glob(osj(run_informations["analysis_folder"], "*statistics.tsv"))
    vcf_files = glob.glob(osj(run_informations["analysis_folder"], "*final.vcf.gz"))
    log.info(f"Cleaning the analysis folder")

    if os.path.isdir(osj(run_informations["run_processing_folder"], "TSV")):
        shutil.rmtree(osj(run_informations["run_processing_folder"], "TSV"))
        os.mkdir(osj(run_informations["run_processing_folder"], "TSV"))

    for file in tsv_files:
        os.chmod(file, 0o777)
        shutil.move(file, run_informations["run_processing_folder_tsv"])
    for file in statistics_files:
        os.chmod(file, 0o777)
        shutil.move(file, run_informations["run_processing_folder_tsv"])
    log.info(f"Moved generated tsv files")

    if os.path.isdir(osj(run_informations["run_processing_folder"], "Alamut")):
        shutil.rmtree(osj(run_informations["run_processing_folder"], "Alamut"))
        log.info(f"Removed old Alamut folder")

    subprocess.run(
        [
            "rsync",
            "-rp",
            osj(run_informations["analysis_folder"], "Alamut"),
            run_informations["run_processing_folder"],
        ]
    )
    log.info(f"Copied new Alamut folder")

    if os.path.isfile(osj(run_informations["run_processing_folder"], "VaRank.log")):
        os.remove(osj(run_informations["run_processing_folder"], "VaRank.log"))
        log.info(f"Removed old VaRank.log")

    subprocess.run(
        [
            "rsync",
            "-rp",
            osj(run_informations["analysis_folder"], "VaRank.log"),
            run_informations["run_processing_folder"],
        ]
    )
    log.info(f"Copied new VaRank.log")

    with open(
        osj(run_informations["run_processing_folder_tsv"], "TSVCopyComplete.txt"),
        mode="a",
    ):
        pass

    for file in vcf_files:
        os.remove(file)
    log.info(f"Deleted used vcf file")

    shutil.rmtree(run_informations["analysis_folder"])
    log.info(f"Deleted temporary analysis folder")


if __name__ == "__main__":
    pass
