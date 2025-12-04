import glob
from operator import truediv
from os.path import join as osj
import os
import logging as log
import datetime
import subprocess


def distribute(run_informations):
    wrapper_log_file = log.getLogger().handlers[0].baseFilename
    run_repository_sample_folder_list = glob.glob(
        osj(run_informations["run_repository"], "*", "")
    )
    d = datetime.datetime.now()
    date = f"{d.strftime('%Y%m%d')}-{d.strftime('%H%M%S')}"

    renamed_wrapper_log_file = f"VANNOT.{date}.analysis.log"
    output_folders = [
        run_informations["run_repository"],
        run_informations["run_depository"],
    ]
    
    sample_count = len(run_repository_sample_folder_list)

    for sample in run_repository_sample_folder_list:
        sample = os.path.basename(os.path.dirname(sample))

        for output_folder in output_folders:
            if not os.path.isdir(osj(output_folder, sample, "VANNOT", "")):
                os.makedirs(osj(output_folder, sample, "VANNOT", ""))
                log.info(
                    f"Created VANNOT folder for sample {sample} into {output_folder}"
                )

        existing_results_files = glob.glob(
            osj(run_informations["run_repository"], sample, "VANNOT", "*")
        )

        for existing_results_file in existing_results_files:
            if os.path.isfile(existing_results_file):
                os.remove(existing_results_file)

        variantannotation_folder_name = f"{os.path.basename(sample)}_{date}"
        results_files = (
            glob.glob(
                osj(
                    run_informations["archives_results_folder"],
                    f"VANNOT_*{sample}*.tsv",
                )
            )
            + glob.glob(
                osj(
                    run_informations["archives_results_folder"], f"VANNOT_*{sample}.log"
                )
            )
            + glob.glob(
                osj(
                    run_informations["archives_results_folder"],
                    f"VANNOT_*{sample}*.vcf.gz",
                )
            )
        )
        if sample_count > 1:
            merged_vcfs = glob.glob(
                osj(
                    run_informations["archives_results_folder"],
                    f"VANNOT_merged_{run_informations["run_name"]}*.vcf.gz",
                )
            )
            for merged_vcf in merged_vcfs:
                if os.path.basename(merged_vcf).split(".")[1] != "design":
                    renamed_merged_vcf = f"VANNOT.{date}.panel.{os.path.basename(merged_vcf).split(".")[2]}.vcf.gz"
                else:
                    renamed_merged_vcf = f"VANNOT.{date}.design.vcf.gz"
                for output_folder in output_folders:
                    subprocess.run(
                        ["rsync", "-rp", merged_vcf, osj(output_folder, renamed_merged_vcf)]
                    )
                    log.info(f"Copied {renamed_merged_vcf} into {output_folder}")
                    
        elif sample_count == 1:
            sample_name = os.path.basename(sample)
            vcfs = glob.glob(
                osj(
                    run_informations["archives_results_folder"],
                    f"VANNOT_{sample_name}*.vcf.gz",
                )
            )
            for vcf in vcfs:
                if os.path.basename(vcf).split(".")[1] != "design":
                    renamed_vcf = f"VANNOT.{date}.panel.{os.path.basename(vcf).split(".")[2]}.vcf.gz"
                else:
                    renamed_vcf = f"VANNOT.{date}.design.vcf.gz"
                for output_folder in output_folders:
                    subprocess.run(
                        ["rsync", "-rp", vcf, osj(output_folder, renamed_vcf)]
                    )
                    log.info(f"Copied {renamed_vcf} into {output_folder}")

        non_redundants = glob.glob(
            osj(run_informations["archives_results_folder"], "Non_Redondant.*.tsv")
        )
        for non_redundant in non_redundants:
            if os.path.basename(non_redundant).split(".")[1] != "design":
                renamed_non_redundant = f"VANNOT.{date}.Non_Redondant.panel.{os.path.basename(non_redundant).split(".")[1]}.tsv"
            else:
                renamed_non_redundant = f"VANNOT.{date}.Non_Redondant.design.tsv"
            for output_folder in output_folders:
                subprocess.run(
                    [
                        "rsync",
                        "-rp",
                        non_redundant,
                        osj(output_folder, renamed_non_redundant),
                    ]
                )
                log.info(f"Copied {renamed_non_redundant} into {output_folder}")

        for results_file in results_files:
            for output_folder in output_folders:
                subprocess.run(
                    [
                        "rsync",
                        "-rp",
                        results_file,
                        osj(
                            output_folder,
                            os.path.basename(sample),
                            "VANNOT",
                            variantannotation_folder_name,
                            "",
                        ),
                    ]
                )
            log.info(
                f"Copied {results_file} into {variantannotation_folder_name} in repository and depository"
            )

            subprocess.run(
                [
                    "rsync",
                    "-rp",
                    results_file,
                    osj(
                        run_informations["run_repository"],
                        os.path.basename(sample),
                        "VANNOT",
                        "",
                    ),
                ]
            )
            for log_file in glob.glob(
                osj(
                    run_informations["run_repository"],
                    os.path.basename(sample),
                    "VANNOT",
                    "*.log",
                )
            ):
                os.remove(log_file)
            log.info(
                f"Copied {results_file} into easy-available results in repository and removed log files"
            )

    for output_folder in output_folders:
        subprocess.run(
            [
                "rsync",
                "-rp",
                wrapper_log_file,
                osj(output_folder, renamed_wrapper_log_file),
            ]
        )
        log.info(f"Copied {renamed_wrapper_log_file} into {output_folder}")

    if os.path.isfile(osj(run_informations["run_repository"], "VANNOTRunning.txt")):
        os.remove(osj(run_informations["run_repository"], "VANNOTRunning.txt"))
        log.info("Deleted VANNOTRunning.txt")


if __name__ == "__main__":
    pass
