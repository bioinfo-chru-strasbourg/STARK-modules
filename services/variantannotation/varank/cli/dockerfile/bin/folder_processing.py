# -*- coding: utf-8 -*-
import logging as log
import checker
import os
from os.path import join as osj
import archives_processing
import non_redundant_generator
import glob
import hgmd_processing


def launch_folder(args):
    folder_path = args.folder
    if folder_path.endswith("/"):
        folder_path = folder_path[:-1]
    folder_path_list = folder_path.split("/")
    in_application = False

    if folder_path.startswith(
        osj(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES"
            ],
            "Archives",
        )
    ):
        for i in glob.glob(
            osj(
                os.environ[
                    "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES"
                ],
                "Archives",
                "*",
                "*",
            )
        ):
            if folder_path_list[-1] == i.split("/")[-1]:
                in_application = True
                run_informations = {
                    "run_name": "none",
                    "run_application": folder_path_list[-1],
                    "run_platform": folder_path_list[-2],
                    "run_platform_application": f"{folder_path_list[-2]}.{folder_path_list[-1]}",
                    "run_depository": "none",
                    "run_pattern": "none",
                    "run_processing_folder": folder_path,
                    "run_repository": "none",
                    "run_processing_folder_tsv": osj(folder_path, "TSV"),
                    "run_processing_folder_vcf_run": "none",
                    "analysis_folder": f"/tmp_{folder_path_list[-2]}.{folder_path_list[-1]}",
                }
        if in_application is False:
            log.error(
                f"You can't start a folder analysis in {folder_path}, however you can start it in one of existing application"
            )
            raise ValueError(folder_path)
        elif in_application is True:
            log.info(
                f"Starting folder analysis from {folder_path_list[-1]} application"
            )

    else:
        run_informations = {
            "run_name": "none",
            "run_application": "none",
            "run_platform": "none",
            "run_platform_application": "default",
            "run_depository": "none",
            "run_pattern": "none",
            "run_processing_folder": folder_path,
            "run_repository": "none",
            "run_processing_folder_tsv": osj(folder_path, "TSV"),
            "run_processing_folder_vcf_run": "none",
            "analysis_folder": f"/tmp_{folder_path_list[-2]}.{folder_path_list[-1]}",
        }
        log.info(f"Starting folder analysis from {folder_path}")

    archives_processing.root_vcf_initialisation(run_informations)
    hgmd_processing.annotate(run_informations)
    checker.configfile(run_informations)
    archives_processing.configfile_manager(run_informations)
    archives_processing.varank_launcher(run_informations)
    archives_processing.cleaner(run_informations)
    checker.logfile(run_informations)
    non_redundant_generator.generate(run_informations)

    log.info(f"VaRank analysis for folder {folder_path} ended well")


if __name__ == "__main__":
    pass
