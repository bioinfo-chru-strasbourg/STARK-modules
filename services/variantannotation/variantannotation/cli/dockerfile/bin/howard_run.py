# -*- coding: utf-8 -*-
import logging as log
import checker
import os
from os.path import join as osj
import synchronizer
import howard_processing
import results_provider

# import non_redundant_generator
import subprocess


def launch_run(args):
    run_repository = args.run
    if run_repository.endswith("/"):
        run_repository = run_repository[:-1]
    run_repository_list = run_repository.split("/")

    run_informations = {
        "assembly": args.assembly,
        "parameters_file": args.param,
        "output_format": args.output_format,
        "run_name": run_repository_list[-1],
        "run_application": run_repository_list[-2],
        "run_platform": run_repository_list[-3],
        "run_platform_application": f"{run_repository_list[-3]}.{run_repository_list[-2]}",
        "run_depository": osj(
            os.environ["DOCKER_DEPOSITORY"],
            run_repository_list[-3],
            run_repository_list[-2],
            run_repository_list[-1],
        ),
        "run_repository": run_repository,
        "vcf_pattern": args.pattern,
        "archives_project_folder": osj(
            os.environ["DOCKER_SERVICES"],
            "Archives",
            run_repository_list[-3],
            run_repository_list[-2],
        ),
        "archives_results_folder": osj(
            os.environ["DOCKER_SERVICES"],
            "Archives",
            run_repository_list[-3],
            run_repository_list[-2],
            "results",
        ),
        "archives_run_folder": osj(
            os.environ["DOCKER_SERVICES"],
            "Archives",
            run_repository_list[-3],
            run_repository_list[-2],
            "VCF",
            run_repository_list[-1],
        ),
        "tmp_analysis_folder": osj(
            os.environ["DOCKER_TMP"],
            f"tmp_{run_repository_list[-1]}/",
        ),
        "module_config": osj(os.environ["DOCKER_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json"),
    }
    checker.depository_checker(run_informations)
    checker.pattern_checker(run_informations)
    variantannotation_running_log = osj(run_repository, "VARunning.txt")

    with open(variantannotation_running_log, "w") as write_file:
        pass

    synchronizer.processing_folder_vcf_synchronizer(run_informations)
    howard_processing.run_initialisation(run_informations)
    howard_processing.cleaner(run_informations)
    # non_redundant_generator.generate(run_informations)
    results_provider.distribute(run_informations)

    log.info(
        f"VariantAnnotation analysis for run {run_informations['run_name']} ended well"
    )


if __name__ == "__main__":
    pass
