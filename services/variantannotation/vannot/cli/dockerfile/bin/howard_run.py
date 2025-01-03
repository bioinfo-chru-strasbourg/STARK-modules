import logging as log
import checker
import os
from os.path import join as osj
import synchronizer
import howard_processing
import results_provider
import dejavu_processing
import non_redundant


def launch_run(args):
    run_repository = args.run
    if run_repository.endswith("/"):
        run_repository = run_repository[:-1]
    run_repository_list = run_repository.split("/")

    run_informations = {
        "assembly": args.assembly,
        "parameters_file": args.param,
        "output_format": args.output_format,
        "type": "run",
        "run_name": run_repository_list[-1],
        "run_application": run_repository_list[-2],
        "run_platform": run_repository_list[-3],
        "run_platform_application": f"{run_repository_list[-3]}.{run_repository_list[-2]}",
        "run_depository": osj(
            os.environ["HOST_DEPOSITORY"],
            run_repository_list[-3],
            run_repository_list[-2],
            run_repository_list[-1],
        ),
        "run_repository": run_repository,
        "vcf_pattern": args.pattern,
        "archives_project_folder": osj(
            os.environ["HOST_SERVICES"],
            "Archives",
            args.assembly,
            run_repository_list[-3],
            run_repository_list[-2],
        ),
        "archives_results_folder": osj(
            os.environ["HOST_SERVICES"],
            "Archives",
            args.assembly,
            run_repository_list[-3],
            run_repository_list[-2],
            "results",
        ),
        "archives_run_folder": osj(
            os.environ["HOST_SERVICES"],
            "Archives",
            args.assembly,
            run_repository_list[-3],
            run_repository_list[-2],
            "VCF",
            run_repository_list[-1],
            "",
        ),
        "parquet_db_run_folder": osj(
            os.environ["HOST_DATABASES"],
            "dejavu",
            "current",
            args.assembly,
            "dejavu.partition.parquet",
            f"GROUP={run_repository_list[-3]}",
            f"PROJECT={run_repository_list[-2]}",
            f"RUN={run_repository_list[-1]}",
        ),
        "parquet_db_project_folder": osj(
            os.environ["HOST_DATABASES"],
            "dejavu",
            "current",
            args.assembly,
            "dejavu.partition.parquet",
            f"GROUP={run_repository_list[-3]}",
            f"PROJECT={run_repository_list[-2]}",
        ),
        "parquet_db_howard_folder": osj(
            os.environ["HOST_DATABASES"],
            "dejavu",
            "current",
            args.assembly,
            "dejavu.partition.parquet",
        ),
        "parquet_db_folder": osj(
            os.environ["HOST_DATABASES"],
            "dejavu",
            "current",
            args.assembly,
            "dejavu.partition.parquet",
        ),
        "tmp_analysis_folder": osj(
            os.environ["HOST_TMP"],
            f"tmp_{run_repository_list[-1]}/",
        ),
        "module_config": osj(
            os.environ["HOST_CONFIG"],
            f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
        ),
    }
    checker.depository_checker(run_informations)
    checker.pattern_checker(run_informations)
    run_informations = checker.panel_checker(run_informations)
    variantannotation_running_log = osj(run_repository, "VANNOTRunning.txt")

    with open(variantannotation_running_log, "w") as write_file:
        pass

    synchronizer.vcf_synchronizer(run_informations)
    dejavu_processing.convert_vcf_parquet(run_informations, args)
    dejavu_processing.calculate_dejavu(run_informations)
    howard_processing.run_initialisation(run_informations)
    merged_vcf = howard_processing.merge_vcf(run_informations, "1", "")
    fambarcode_vcf = howard_processing.fambarcode_vcf(
        run_informations,
        merged_vcf,
    )
    annotated_merged_vcf = howard_processing.howard_proc(
        run_informations, fambarcode_vcf
    )
    howard_processing.unmerge_vcf(annotated_merged_vcf, run_informations)
    howard_processing.howard_score_transcripts(run_informations)
    howard_processing.gmc_score(run_informations)
    print(howard_processing.merge_vcf(run_informations, "2", ""))
    if run_informations["run_panels"] != "":
        howard_processing.panel_filtering(run_informations)
    howard_processing.convert_to_final_tsv(run_informations)

    # non_redundant.generate(run_informations)
    howard_processing.cleaner(run_informations)
    results_provider.distribute(run_informations)

    lock_file = osj(run_repository, "VANNOTComplete.txt")
    with open(lock_file, "w") as write_file:
        pass

    log.info(f"vannot analysis for run {run_informations['run_name']} ended well")


if __name__ == "__main__":
    pass
