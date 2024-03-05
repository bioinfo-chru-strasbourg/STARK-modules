# -*- coding: utf-8 -*-
import glob
import logging as log
from os.path import join as osj
import os
import commons
import subprocess
import re



def annotate(run_informations):

    analysis_folder_vcf_files = glob.glob(
        osj(run_informations["analysis_folder"], "*.vcf.gz")
    )
    hgmd_db = osj(
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_DATABASES"],
        "HGMD",
        "current",
        "hg19",
        "HGMD_Pro_2023.4_hg19.clean2.vcf.gz",
    )

    for vcf_file in analysis_folder_vcf_files:
        vcf_name = os.path.basename(vcf_file)
        vcfout = osj(os.path.dirname(vcf_file), f"hgmd_{vcf_name}")
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