# -*- coding: utf-8 -*-
import glob
from operator import truediv
from os.path import join as osj
import os
import re
import logging as log
import datetime
import subprocess


def distribute(run_informations):
    varank_log_file = osj(run_informations["run_processing_folder"], "VaRank.log")
    wrapper_log_file = log.getLogger().handlers[0].baseFilename
    non_redundant = osj(
        run_informations["run_processing_folder_tsv"], "Non_Redondant.tsv"
    )
    run_repository_sample_folder_list = glob.glob(
        osj(run_informations["run_repository"], "*", "")
    )
    d = datetime.datetime.now()
    date = f"{d.strftime('%Y%m%d')}-{d.strftime('%H%M%S')}"

    renamed_varank_log_file = f"VARANK.{date}.VaRank.log"
    renamed_wrapper_log_file = f"VARANKComplete.txt"
    renamed_non_redundant = f"VARANK.{date}.Non_Redondant.tsv"

    subprocess.run(
        [
            "rsync",
            "-rp",
            varank_log_file,
            osj(run_informations["run_repository"], renamed_varank_log_file),
        ]
    )
    subprocess.run(
        [
            "rsync",
            "-rp",
            varank_log_file,
            osj(run_informations["run_depository"], renamed_varank_log_file),
        ]
    )
    subprocess.run(
        [
            "rsync",
            "-rp",
            non_redundant,
            osj(run_informations["run_repository"], renamed_non_redundant),
        ]
    )
    subprocess.run(
        [
            "rsync",
            "-rp",
            non_redundant,
            osj(run_informations["run_depository"], renamed_non_redundant),
        ]
    )

    for sample in run_repository_sample_folder_list:
        sample = os.path.basename(os.path.dirname(sample))

        if not os.path.isdir(
            osj(run_informations["run_repository"], sample, "VARANK", "")
        ):
            os.mkdir(osj(run_informations["run_repository"], sample, "VARANK", ""))
        if not os.path.isdir(
            osj(run_informations["run_depository"], sample, "VARANK", "")
        ):
            os.mkdir(osj(run_informations["run_depository"], sample, "VARANK", ""))

        varank_folder_name = f"{os.path.basename(sample)}_{date}"
        tsv_files = glob.glob(
            osj(run_informations["run_processing_folder_tsv"], f"*{sample}*.tsv")
        )

        for tsv_file in tsv_files:
            subprocess.run(
                [
                    "rsync",
                    "-rp",
                    tsv_file,
                    osj(
                        run_informations["run_repository"],
                        os.path.basename(sample),
                        "VARANK",
                        varank_folder_name,
                        "",
                    ),
                ]
            )
            subprocess.run(
                [
                    "rsync",
                    "-rp",
                    tsv_file,
                    osj(
                        run_informations["run_depository"],
                        os.path.basename(sample),
                        "VARANK",
                        varank_folder_name,
                        "",
                    ),
                ]
            )

        subprocess.run(
            [
                "rsync",
                "-rp",
                wrapper_log_file,
                osj(run_informations["run_repository"], renamed_wrapper_log_file),
            ]
        )
        subprocess.run(
            [
                "rsync",
                "-rp",
                wrapper_log_file,
                osj(run_informations["run_depository"], renamed_wrapper_log_file),
            ]
        )
        os.remove(osj(run_informations["run_repository"], "VARANKRunning.txt"))


if __name__ == "__main__":
    pass
