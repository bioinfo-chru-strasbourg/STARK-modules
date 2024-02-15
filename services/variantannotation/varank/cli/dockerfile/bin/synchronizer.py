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
    pattern = run_informations["run_pattern"]
    for i in pattern:
        log.info(f"Syncronizing vcf files according to pattern {i}")

    if not os.path.isdir(run_informations["run_processing_folder_vcf_run"]):
        os.makedirs(run_informations["run_processing_folder_vcf_run"])
        os.chmod(run_informations["run_processing_folder_vcf_run"], 0o777)

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
                            run_informations["run_processing_folder_vcf_run"],
                        ]
                    )
        else:
            for pattern_vcf_file in vcf_files:
                subprocess.run(
                    [
                        "rsync",
                        "-rp",
                        pattern_vcf_file,
                        run_informations["run_processing_folder_vcf_run"],
                    ]
                )
                log.info(
                    f'Overwriting {os.path.basename(pattern_vcf_file)} with vcf from "{element}" pattern'
                )

    archives_vcf_files = glob.glob(
        osj(run_informations["run_processing_folder_vcf_run"], "*")
    )
    hgmd_db = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DATABASES"],
        "HGMD",
        "current",
        "hg19",
        "HGMD_Pro_2023.4_hg19.clean2.vcf.gz",
    )
    for vcf_file in archives_vcf_files:
        vcf_name = os.path.basename(vcf_file)
        vcfout = osj(
            run_informations["run_processing_folder_vcf_run"], f"hgmd_{vcf_name}"
        )
        subprocess.run(
            [
                "python3",
                "/tools/varank/VaRank/scripts/HGMD_annotations.py",
                "-i",
                vcf_file,
                "-o",
                vcfout,
                "-a",
                "hg19",
                "--hgmd",
                hgmd_db,
                "--bcftools",
                "/tools/bin/bcftools",
            ]
        )
        os.remove(vcf_file)
        os.remove(vcf_file + ".tbi")
        os.chmod(vcfout, 0o777)
