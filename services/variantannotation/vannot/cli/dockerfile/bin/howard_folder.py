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

def launch_folder(args):
    analysis_folder = args.folder
    if analysis_folder.endswith("/"):
        analysis_folder = analysis_folder[:-1]
    run_repository_list = analysis_folder.split("/")

    run_informations = {
        "assembly": args.assembly,
        "parameters_file": args.param,
        "output_format": args.output_format,
        "run_name": run_repository_list[-1],
        "run_application": "",
        "run_platform": "",
        "run_platform_application": "",
        "run_depository": "",
        "run_repository": "",
        "vcf_pattern": "",
        "archives_project_folder": "",
        "archives_results_folder": osj(analysis_folder, "results"),
        "archives_run_folder": "",
        "parquet_db_run_folder": "",
        "parquet_db_project_folder": "",
        "parquet_db_howard_folder": "",
        "parquet_db_folder": "",
        "tmp_analysis_folder": osj(
            os.environ["DOCKER_TMP"],
            f"tmp_{run_repository_list[-1]}/",
        ),
        "module_config": osj(os.environ["DOCKER_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json"),
    }
    howard_processing.run_initialisation(run_informations)
    merged_vcf = howard_processing.merge_vcf(run_informations)
    annotated_merged_vcf = howard_processing.howard_proc(run_informations, merged_vcf)
    howard_processing.unmerge_vcf(annotated_merged_vcf)
    howard_processing.convert_to_final_tsv(run_informations)
    non_redundant.generate(run_informations)

    howard_processing.cleaner(run_informations)

    log.info(
        f"vannot analysis for folder {run_informations['run_name']} ended well"
    )


if __name__ == "__main__":
    pass
