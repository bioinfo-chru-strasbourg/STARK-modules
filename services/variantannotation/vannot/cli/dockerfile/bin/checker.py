import os
import json
import logging as log
import glob
from os.path import join as osj
import commons


def absolute_folder_path(path):
    if (
        os.path.isabs(path)
        and os.path.isdir(path)
    ):
        return path
    elif (
        not os.path.isabs(path)
        and os.path.isdir(path)
    ):
        return os.path.abspath(path)
    else:
        raise ValueError(path)

def platform_mapping(run_informations):
    run_platform = run_informations["run_platform"]
    module_config = osj(
        os.environ["HOST_MODULE_CONFIG"],
        f"{os.environ['DOCKER_SUBMODULE_NAME']}_config.json",
    )
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        if run_platform in data["platform_mapping"]:
            mapped_platform = data["platform_mapping"][run_platform]
            log.info(
                f"Mapping {run_platform} to {mapped_platform} for the analysis"
            )
            run_informations.update({"run_platform": mapped_platform})
            run_informations.update({"run_platform_application": f"{mapped_platform}.{run_informations['run_application']}"})

            run_application = run_informations["run_application"]
            run_name = run_informations["run_name"]
            assembly = run_informations["assembly"]
            run_informations.update({
                "archives_project_folder": osj(
                    os.environ["HOST_SERVICES"],
                    "Archives",
                    assembly,
                    mapped_platform,
                    run_application,
                ),
                "archives_results_folder": osj(
                    os.environ["HOST_SERVICES"],
                    "Archives",
                    assembly,
                    mapped_platform,
                    run_application,
                    "results",
                ),
                "archives_run_folder": osj(
                    os.environ["HOST_SERVICES"],
                    "Archives",
                    assembly,
                    mapped_platform,
                    run_application,
                    "VCF",
                    run_name,
                    "",
                ),
                "parquet_db_run_folder": osj(
                    os.environ["HOST_DATABASES"],
                    "dejavu",
                    "current",
                    assembly,
                    "dejavu.partition.parquet",
                    f"GROUP={mapped_platform}",
                    f"PROJECT={run_application}",
                    f"RUN={run_name}",
                ),
                "parquet_db_project_folder": osj(
                    os.environ["HOST_DATABASES"],
                    "dejavu",
                    "current",
                    assembly,
                    "dejavu.partition.parquet",
                    f"GROUP={mapped_platform}",
                    f"PROJECT={run_application}",
                ),
                "parquet_db_howard_folder": osj(
                    os.environ["HOST_DATABASES"],
                    "dejavu",
                    "current",
                    assembly,
                    "dejavu.partition.parquet",
                ),
            })
        else:
            log.info(
                f"No mapping found for {run_platform}, using it as is for the analysis"
            )
        
def absolute_run_path(path):
    if (
        os.path.isabs(path)
        and os.path.isdir(path)
        and path.startswith(os.environ["HOST_REPOSITORY"])
        and "." not in os.path.basename(path)
    ):
        return path
    elif (
        not os.path.isabs(path)
        and os.path.isdir(path)
        and os.path.abspath(path).startswith(os.environ["HOST_REPOSITORY"])
        and "." not in os.path.basename(path)
    ):
        return os.path.abspath(path)
    else:
        raise ValueError(path)


def depository_checker(run_informations):
    run_repository = run_informations["run_repository"]
    run_depository = run_informations["run_depository"]

    if not os.path.isdir(run_depository):
        log.warning(
            "Specified depository folder doesn't exists, maybe it was sent to the Archives ? Creating a new folder with all subfolders tree."
        )
        run_repository_sample_folders = glob.glob(osj(run_repository, "*", ""))
        os.makedirs(run_depository, 0o775)

        for sample_folder in run_repository_sample_folders:
            sample_folder = osj(run_depository, os.path.basename(sample_folder[:-1]))
            if not os.path.isdir(sample_folder):
                os.mkdir(sample_folder, 0o755)


def pattern_checker(run_informations):
    if run_informations["type"] == "run_dejavu" and run_informations["run_repository"] == "":
        run_path = run_informations["run_archives"]
    else:
        run_path = run_informations["run_repository"]

    pattern = run_informations["vcf_pattern"]
    for element in pattern:
        vcf_files = glob.glob(osj(run_path, element))
        if len(vcf_files) == 0 and element != commons.get_default_pattern(run_informations):
            log.error(
                f"There is no vcf files with the specified pattern {element}, please check your command-line"
            )
            raise ValueError(element)

        elif len(vcf_files) == 0 and element == commons.get_default_pattern(run_informations):
            log.error(
                f"There is no vcf files with the default STARK analysis pattern {element}, please check the analysis integrity"
            )
            raise ValueError(element)


def panel_checker(run_informations):
    panels = find_panel(run_informations)
    run_informations.update({"run_panels": panels})
    log.info(f"Found {len(run_informations["run_panels"])} panels in your run, using them to filter results at the end of the process")
    return run_informations


def find_panel(run_informations):
    samples = glob.glob(osj(run_informations["run_repository"], "*", ""))
    panels = glob.glob(osj(samples[0], "STARK", "*.manifest.*genes"))
    return panels


if __name__ == "__main__":
    pass
