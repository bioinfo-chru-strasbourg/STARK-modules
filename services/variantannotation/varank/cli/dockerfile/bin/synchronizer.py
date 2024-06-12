# -*- coding: utf-8 -*-
import glob
import logging as log
from os.path import join as osj
import os
import commons
import subprocess
import re


def processing_folder_vcf_synchronizer(run_informations):
    run_repository = run_informations["run_repository"]
    pattern = run_informations["vcf_pattern"]
    for i in pattern:
        log.info(f"Syncronizing vcf files according to pattern {i}")

    if not os.path.isdir(run_informations["archives_run_folder"]):
        os.makedirs(run_informations["archives_run_folder"])
        os.chmod(run_informations["archives_run_folder"], 0o777)

    for element in pattern:
        vcf_files = glob.glob(osj(run_repository, element))
        if element == commons.default_pattern:
            for stark_vcf_file in vcf_files:
                if re.match(
                    run_repository
                    + "\\/.+\\/STARK\\/.+\\.reports\\/[^.]+.final.vcf.gz",
                    stark_vcf_file,
                ):
                    log.info(
                        f'Syncronizing {os.path.basename(stark_vcf_file)} from "{element}" pattern'
                    )
                    subprocess.run(
                        [
                            "rsync",
                            "-rp",
                            stark_vcf_file,
                            run_informations["archives_run_folder"],
                        ]
                    )
        else:
            for pattern_vcf_file in vcf_files:
                subprocess.run(
                    [
                        "rsync",
                        "-rp",
                        pattern_vcf_file,
                        run_informations["archives_run_folder"],
                    ]
                )
                log.info(
                    f'Overwriting {os.path.basename(pattern_vcf_file)} with vcf from "{element}" pattern'
                )
