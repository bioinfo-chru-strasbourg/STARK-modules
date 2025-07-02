# -*- coding: utf-8 -*-
import logging as log
import checker
import os
from os.path import join as osj
import synchronizer
import archives_processing
import non_redundant_generator
import results_provider
import subprocess
import hgmd_processing


def launch_run(args):
    run_repository = args.run
    if run_repository.endswith("/"):
        run_repository = run_repository[:-1]
    run_repository_list = run_repository.split("/")

    run_informations = {
        "run_name": run_repository_list[-1],
        "run_application": run_repository_list[-2],
        "run_platform": run_repository_list[-3],
        "run_platform_application": f"{run_repository_list[-3]}.{run_repository_list[-2]}",
        "run_depository": osj(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DEPOSITORY"
            ],
            run_repository_list[-3],
            run_repository_list[-2],
            run_repository_list[-1],
        ),
        "vcf_pattern": args.pattern,
        "archives_project_folder": osj(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES"
            ],
            "Archives",
            run_repository_list[-3],
            run_repository_list[-2],
        ),
        "run_repository": run_repository,
        "run_processing_folder_tsv": osj(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES"
            ],
            "Archives",
            run_repository_list[-3],
            run_repository_list[-2],
            "TSV",
        ),
        "archives_run_folder": osj(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES"
            ],
            "Archives",
            run_repository_list[-3],
            run_repository_list[-2],
            "VCF",
            run_repository_list[-1],
        ),
        "tmp_analysis_folder": f"/tmp_{run_repository_list[-1]}/",
    }
    checker.depository_checker(run_informations)
    checker.pattern_checker(run_informations)
    varank_running_log = osj(run_repository, "VARANKRunning.txt")

    with open(varank_running_log, "w") as write_file:
        pass

    synchronizer.processing_folder_vcf_synchronizer(run_informations)
    archives_processing.run_initialisation(run_informations)
    hgmd_processing.annotate(run_informations)
    checker.configfile(run_informations)
    archives_processing.configfile_manager(run_informations)
    archives_processing.varank_launcher(run_informations)
    archives_processing.cleaner(run_informations)
    checker.logfile(run_informations)
    non_redundant_generator.generate(run_informations)
    results_provider.distribute(run_informations)

    log.info(f"VaRank analysis for run {run_informations['run_name']} ended well")


if __name__ == "__main__":
    pass
