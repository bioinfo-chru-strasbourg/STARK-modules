import logging as log
import os
from os.path import join as osj
import commons
import dejavu_processing
import checker
import synchronizer

def launch_dejavu(args):
    if args.dejavu is not None:
        analysis_project = args.dejavu
        if analysis_project.endswith("/"):
            analysis_project = analysis_project[:-1]
        analysis_project_name = analysis_project.split("/")

        run_informations = {
            "assembly": args.assembly,
            "onco": args.onco,
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
            "force": args.force,
            "archives_project_folder": osj(
                os.environ["HOST_SERVICES"],
                "Archives",
                args.assembly,
                analysis_project_name[-2],
                analysis_project_name[-1],
            ),
            "archives_results_folder": "",
            "archives_run_folder": "",
            "parquet_db_run_folder": "",
            "parquet_db_project_folder": osj(
                os.environ["HOST_DATABASES"],
                "dejavu",
                "current",
                args.assembly,
                "dejavu.partition.parquet",
                f"GROUP={analysis_project_name[-2]}",
                f"PROJECT={analysis_project_name[-1]}",
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
                f"tmp_dejavu_{analysis_project_name[-1]}/",
            ),
            "module_config": "",
        }

    if args.run_dejavu is not None:
        if args.run_dejavu.endswith("/"):
            run_path = args.run_dejavu[:-1]
        else:
            run_path = args.run_dejavu
        run_path_list = run_path.split("/")

        run_informations = {
            "assembly": args.assembly,
            "onco": args.onco,
            "parameters_file": "",
            "output_format": "",
            "type": "run_dejavu",
            "force": args.force,
            "run_name": run_path_list[-1],
            "run_application": run_path_list[-2],
            "run_platform": run_path_list[-3],
            "run_platform_application": f"{run_path_list[-3]}.{run_path_list[-2]}",
            "run_repository": "",
            "run_archives": "",
            "run_depository": "",
            "vcf_pattern": args.pattern,
            "archives_project_folder": osj(
                os.environ["HOST_SERVICES"],
                "Archives",
                args.assembly,
                run_path_list[-2],
                run_path_list[-1],
            ),
            "archives_results_folder": "",
            "archives_run_folder": osj(
                os.environ["HOST_SERVICES"],
                "Archives",
                args.assembly,
                run_path_list[-3],
                run_path_list[-2],
                "VCF",
                run_path_list[-1],
                "",
            ),
            "parquet_db_run_folder": osj(
                os.environ["HOST_DATABASES"],
                "dejavu",
                "current",
                args.assembly,
                "dejavu.partition.parquet",
                f"GROUP={run_path_list[-3]}",
                f"PROJECT={run_path_list[-2]}",
                f"RUN={run_path_list[-1]}",
            ),
            "parquet_db_project_folder": osj(
                os.environ["HOST_DATABASES"],
                "dejavu",
                "current",
                args.assembly,
                "dejavu.partition.parquet",
                f"GROUP={run_path_list[-3]}",
                f"PROJECT={run_path_list[-2]}",
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
                f"tmp_{run_path_list[-1]}/",
            ),
            "module_config": osj(
                os.environ["HOST_CONFIG"],
                f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
            ),
        }

        if run_path_list[-4] == "archives":
            run_informations["run_archives"] = osj(
                os.environ["HOST_ARCHIVES"],
                run_path_list[-3],
                run_path_list[-2],
                run_path_list[-1],
            )
            if len(run_informations["vcf_pattern"]) == 1:
                run_informations["vcf_pattern"] = ["*/*.final.vcf.gz"]
            else:
                for i in range(len(run_informations["vcf_pattern"])):
                    if run_informations["vcf_pattern"][i] == commons.get_default_pattern("default"):
                        run_informations["vcf_pattern"][i] = "*/*.final.vcf.gz"

        elif run_path_list[-4] == "repository":
            run_informations["run_repository"] = osj(
                os.environ["HOST_REPOSITORY"],
                run_path_list[-3],
                run_path_list[-2],
                run_path_list[-1],
            )

        checker.platform_mapping(run_informations)

        if os.path.isdir(run_informations["parquet_db_run_folder"]) and args.force is False:
            log.info(
                f"Parquet database for run {run_informations['run_name']} already exists, skipping dejavu generation"
            )
        elif (os.path.isdir(run_informations["parquet_db_run_folder"]) and args.force is True) or not os.path.isdir(run_informations["parquet_db_run_folder"]):
            checker.pattern_checker(run_informations)
            synchronizer.vcf_synchronizer(run_informations)
            dejavu_processing.convert_vcf_parquet(run_informations, args)

    if args.dejavu is not None:
        dejavu_processing.convert_vcf_parquet(run_informations, args)

    if args.global_dejavu is True:
        dejavu_processing.calculate_dejavu(run_informations)

    log.info(
        f"vannot dejavu analysis for folder {run_informations['run_name']} ended well"
    )


if __name__ == "__main__":
    pass
