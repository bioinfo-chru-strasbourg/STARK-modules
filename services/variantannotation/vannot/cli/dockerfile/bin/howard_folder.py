import logging as log
import os
from os.path import join as osj
import howard_processing
import non_redundant


def launch_folder(args):
    analysis_folder = args.folder
    if analysis_folder.endswith("/"):
        analysis_folder = analysis_folder[:-1]
    analysis_folder_name = analysis_folder.split("/")

    run_informations = {
        "assembly": args.assembly,
        "parameters_file": args.param,
        "output_format": args.output_format,
        "type": "folder",
        "run_name": analysis_folder_name[-1],
        "run_application": "",
        "run_platform": "",
        "run_platform_application": "default",
        "run_depository": "",
        "run_repository": "",
        "vcf_pattern": "",
        "archives_project_folder": "",
        "archives_results_folder": osj(analysis_folder, "results"),
        "archives_run_folder": analysis_folder,
        "parquet_db_run_folder": "",
        "parquet_db_project_folder": "",
        "parquet_db_howard_folder": "",
        "parquet_db_folder": "",
        "tmp_analysis_folder": osj(
            os.environ["HOST_TMP"],
            f"tmp_{analysis_folder_name[-1]}/",
        ),
        "module_config": osj(
            os.environ["DOCKER_CONFIG"],
            f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
        ),
    }
    if os.environ["DOCKER_SERVICES"] in analysis_folder:
        run_informations["run_application"] = analysis_folder_name[-1]
        run_informations["run_platform"] = analysis_folder_name[-2]
        run_informations["run_platform_application"] = (
            f"{analysis_folder_name[-2]}.{analysis_folder_name[-1]}"
        )

        howard_processing.project_folder_initialisation(run_informations)
        merged_vcf = howard_processing.merge_vcf(run_informations)
        annotated_merged_vcf = howard_processing.howard_proc(
            run_informations, merged_vcf
        )
        howard_processing.unmerge_vcf(annotated_merged_vcf, run_informations)
        howard_processing.howard_score_transcripts(run_informations)
        howard_processing.convert_to_final_tsv(run_informations)
        non_redundant.generate(run_informations)
        howard_processing.cleaner(run_informations)

    else:
        howard_processing.folder_initialisation(run_informations)
        merged_vcf = howard_processing.merge_vcf(run_informations)
        annotated_merged_vcf = howard_processing.howard_proc(
            run_informations, merged_vcf
        )
        howard_processing.unmerge_vcf(annotated_merged_vcf, run_informations)
        howard_processing.howard_score_transcripts(run_informations)
        howard_processing.convert_to_final_tsv(run_informations)
        non_redundant.generate(run_informations)
        howard_processing.cleaner(run_informations)

    log.info(f"vannot analysis for folder {run_informations['run_name']} ended well")


if __name__ == "__main__":
    pass
