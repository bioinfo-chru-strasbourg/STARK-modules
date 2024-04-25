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


def initialisation(run_informations):
    archives_vcf_files = glob.glob(
        osj(run_informations["run_processing_folder"], "VCF", "*", "*.final.vcf*")
    )
    deleted_samples_file = osj(
        os.environ[
            "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARIANTANNOTATION_FOLDER_CONFIG"
        ],
        "variantannotation",
        "variantannotation",
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
        vcf_file_list = glob.glob(
            osj(run_informations["run_processing_folder_vcf_run"], "*.final.vcf*")
        )
        for processed_vcf_file in vcf_file_list:
            subprocess.run(
                [
                    "rsync",
                    "-rp",
                    processed_vcf_file,
                    run_informations["analysis_folder"],
                ]
            )

    if os.path.isdir(run_informations["run_processing_folder"]):
        log.info(f"Cleaning old files in {run_informations['run_processing_folder']}")
        with open(
            osj(run_informations["run_processing_folder"], "RESULTSDeleting.txt"),
            mode="a",
        ):
            pass
        # for results_file in glob.glob(
        #     osj(run_informations["run_processing_folder"], "*")
        # ):
        #     os.remove(results_file)
        os.remove(osj(run_informations["run_processing_folder"], "RESULTSDeleting.txt"))
    else:
        os.mkdir(run_informations["run_processing_folder"])
    vcf_file_to_analyse = glob.glob(
        osj(run_informations["analysis_folder"], "*.final.vcf*")
    )
    for vcf_file in vcf_file_to_analyse:
        howard_launcher(run_informations, vcf_file)


def howard_launcher(run_informations, vcf_file):
    print(f"launching for {vcf_file}")
    logfile = osj(run_informations["analysis_folder"], "variantannotation.log")
    output_file = osj(
        run_informations["analysis_folder"], f"output_{os.path.basename(vcf_file)}"
    )
    log.info(f"Generating results files")
    with open(logfile, "w") as f:
        subprocess.call(
            [
                "howard",
                "annotation",
                "--input",
                vcf_file,
                "--output",
                output_file,
                "--annotations",
                "/HOWARD/databases/dbsnp/current/hg19/b156/dbsnp.parquet",
            ],
            stdout=f,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
        )
    os.chmod(logfile, 0o777)


def cleaner(run_informations):
    tsv_files = glob.glob(osj(run_informations["analysis_folder"], "*rankingByVar.tsv"))
    statistics_files = glob.glob(
        osj(run_informations["analysis_folder"], "*statistics.tsv")
    )
    vcf_files = glob.glob(osj(run_informations["analysis_folder"], "*final.vcf.gz"))
    log.info(f"Cleaning the analysis folder")

    if os.path.isdir(osj(run_informations["run_processing_folder"], "TSV")):
        shutil.rmtree(osj(run_informations["run_processing_folder"], "TSV"))
        os.mkdir(osj(run_informations["run_processing_folder"], "TSV"))
    else:
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
