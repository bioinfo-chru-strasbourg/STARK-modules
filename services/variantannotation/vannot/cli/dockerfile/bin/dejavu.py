# -*- coding: utf-8 -*-
import logging as log
import checker
import os
from os.path import join as osj
import synchronizer
import howard_processing
import results_provider
import dejavu_processing
import non_redundant


def launch_dejavu(args):
    analysis_project = args.dejavu
    if analysis_project.endswith("/"):
        analysis_project = analysis_project[:-1]
    analysis_project_name = analysis_project.split("/")

    run_informations = {
        "assembly": args.assembly,
        "parameters_file": "",
        "output_format": "",
        "type": "dejavu",
        "run_name": analysis_project_name[-1],
        "run_application": analysis_project_name[-1],
        "run_platform": analysis_project_name[-2],
        "run_platform_application": f"{analysis_project_name[-2]}.{analysis_project_name[-1]}",
        "run_depository": "",
        "run_repository": "",
        "vcf_pattern": "",
        "archives_project_folder": osj(
            os.environ["DOCKER_SERVICES"],
            "Archives",
            args.assembly,
            analysis_project_name[-2],
            analysis_project_name[-1],
        ),
        "archives_results_folder": "",
        "archives_run_folder": "",
        "parquet_db_run_folder": "",
        "parquet_db_project_folder": osj(
            os.environ["DOCKER_DEJAVU_DATABASE"],
            "current",
            args.assembly,
            "dejavu.partition.parquet",
            f"GROUP={analysis_project_name[-2]}",
            f"PROJECT={analysis_project_name[-1]}",
        ),
        "parquet_db_howard_folder": osj(
            "/",
            "databases",
            "dejavu",
            "current",
            args.assembly,
            "dejavu.partition.parquet",
        ),
        "parquet_db_folder": osj(
            os.environ["DOCKER_DEJAVU_DATABASE"],
            "current",
            args.assembly,
            "dejavu.partition.parquet",
        ),
        "tmp_analysis_folder": osj(
            os.environ["DOCKER_TMP"],
            f"tmp_dejavu_{analysis_project_name[-1]}/",
        ),
        "module_config": "",
    }

    dejavu_processing.convert_vcf_parquet(run_informations, args)
    dejavu_processing.calculate_dejavu(run_informations)
    os.remove(run_informations["tmp_analysis_folder"])

    log.info(
        f"vannot dejavu analysis for folder {run_informations['run_name']} ended well"
    )


if __name__ == "__main__":
    pass
