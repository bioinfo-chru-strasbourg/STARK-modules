# -*- coding: utf-8 -*-
import logging as log
import checker
import os
from os.path import join as osj
import synchronizer
import archives_processing
import non_redundant_generator


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
        "run_pattern": args.pattern,
        "run_processing_folder": osj(
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
        "run_repository": run_repository,
        "run_processing_folder_vcf_run": osj(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES"
            ],
            "Archives",
            run_repository_list[-3],
            run_repository_list[-2],
            "VCF",
            run_repository_list[-1],
        ),
    }
    checker.depository_checker(run_informations)
    checker.pattern_checker(run_informations)
    varank_running_log = osj(run_repository, "VARANKRunning.txt")

    with open(varank_running_log, "w") as write_file:
        pass

    synchronizer.processing_folder_vcf_synchronizer(run_informations)
    archives_processing.initialisation(run_informations)
    checker.configfile(run_informations)
    archives_processing.configfile_manager(run_informations)
    archives_processing.varank_launcher(run_informations)
    checker.logfile(run_informations)
    archives_processing.cleaner(run_informations)
    non_redundant_generator.generate(run_informations)


if __name__ == "__main__":
    pass
