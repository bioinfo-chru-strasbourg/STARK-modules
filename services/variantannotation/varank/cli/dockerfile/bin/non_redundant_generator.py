# -*- coding: utf-8 -*-
import glob
from operator import truediv
from os.path import join as osj
import os
import re
import logging as log
import datetime


def generate(run_informations):
    input_dir = run_informations["run_processing_folder_tsv"]
    columns_to_remove = []
    d = datetime.datetime.now()
    output_file = osj(
        input_dir,
        "Non_Redondant"
        + "_"
        + d.strftime("%d")
        + "_"
        + d.strftime("%m")
        + "_"
        + d.strftime("%Y")
        + ".tsv",
    )

    if "WES" in input_dir:
        filter_mode = True
        log.warning(
            f"WES Folder detected, some variants will be filtered in the non-redundant : gnomadAltFreq_popmax_index < 0.005 & gnomadHomCount_all_index < 2"
        )
    else:
        filter_mode = False

    with open(
        osj(
            os.environ[
                "DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_CONFIG"
            ],
            "variantannotation",
            "varank",
            "configfiles",
            "non_redundant_config.txt",
        ),
        "r",
    ) as read_file:
        for line in read_file.readlines():
            line = line.strip()
            if re.match("UNWANTED_VARANK_COLUMN.default", line):
                line = re.sub(r"^.*?=", "", line)
                columns_to_remove = line.split(",")
            elif re.match(
                "UNWANTED_VARANK_COLUMN."
                + run_informations["run_platform_application"],
                line,
            ):
                line = re.sub(r"^.*?=", "", line)
                columns_to_remove = line.split(",")

    first_file = True
    index_to_keep = []
    input_files_list = glob.glob(
        osj(
            input_dir,
            "*_allVariants.rankingByVar.tsv",
        )
    )

    already_seen_variants_dic = {}

    with open(output_file + "_unsorted", "w") as write_file:
        for i, f in enumerate(input_files_list):
            with open(f, "r") as read_file:
                for l in read_file:
                    if l.startswith("#"):
                        continue
                    elif l in ["\n", "\r\n"]:
                        continue
                    elif l.startswith("variantID"):
                        if first_file:
                            header = l
                            l = l.rstrip("\r\n").split("\t")
                            split_header_length = len(l)
                            for i in range(len(l)):
                                if l[i] not in columns_to_remove:
                                    index_to_keep.append(i)
                            write_file.write(
                                "\t".join(l[i] for i in index_to_keep) + "\r\n"
                            )
                            gnomadAltFreq_popmax_index = l.index("gnomadAltFreq_popmax")
                            gnomadHomCount_all_index = l.index("gnomadHomCount_all")
                            first_file = False
                        else:
                            assert l == header


if __name__ == "__main__":
    pass
