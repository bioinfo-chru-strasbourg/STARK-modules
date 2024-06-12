# -*- coding: utf-8 -*-
import logging as log
import checker
import os
from os.path import join as osj
import howard_processing

# import non_redundant_generator
import glob


def launch_folder(args):
    folder_path = args.folder
    if folder_path.endswith("/"):
        folder_path = folder_path[:-1]
    folder_path_list = folder_path.split("/")
    in_application = False

    if folder_path.startswith(osj(os.environ["DOCKER_SERVICES"], "Archives")):
        for i in glob.glob(osj(os.environ["DOCKER_SERVICES"], "Archives", "*", "*")):
            if folder_path_list[-1] == i.split("/")[-1]:
                in_application = True
                run_informations = {
                    "assembly": args.assembly,
                    "parameters_file": args.param,
                    "output_format": args.output_format,
                    "run_name": "none",
                    "run_application": folder_path_list[-1],
                    "run_platform": folder_path_list[-2],
                    "run_platform_application": f"{folder_path_list[-2]}.{folder_path_list[-1]}",
                    "run_depository": "none",
                    "run_repository": "none",
                    "vcf_pattern": "none",
                    "archives_project_folder": folder_path,
                    "archives_results_folder": osj(folder_path, "results"),
                    "archives_run_folder": "none",
                    "tmp_analysis_folder": osj(os.environ["DOCKER_TMP"], folder_path_list[-1]),
                    "module_config": osj(os.environ["DOCKER_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json"),
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
            "assembly": args.assembly,
            "parameters_file": args.param,
            "output_format": args.output_format,
            "run_name": "none",
            "run_application": "none",
            "run_platform": "none",
            "run_platform_application": "default",
            "run_depository": "none",
            "vcf_pattern": "none",
            "archives_project_folder": folder_path,
            "run_repository": "none",
            "archives_results_folder": osj(folder_path, "results"),
            "archives_run_folder": "none",
            "tmp_analysis_folder": osj(os.environ["DOCKER_TMP"], folder_path_list[-1]),
            "module_config": osj(os.environ["DOCKER_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json"),
        }
        log.info(f"Starting folder analysis from {folder_path}")

    howard_processing.folder_initialisation(run_informations)
    howard_processing.cleaner(run_informations)
    # non_redundant_generator.generate(run_informations)

    log.info(f"VariantAnnotation analysis for folder {folder_path} ended well")


if __name__ == "__main__":
    pass
