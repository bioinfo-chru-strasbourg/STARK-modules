# -*- coding: utf-8 -*-
import glob
from operator import truediv
from os.path import join as osj
import os
import re
import logging as log
import datetime
import subprocess


def generate(run_informations):
    # Adapted from the Samuel's non_redondant.py script
    input_dir = run_informations["run_processing_folder_tsv"]
    columns_to_remove = []
    d = datetime.datetime.now()
    output_file = osj(input_dir, "Non_Redondant.tsv")

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
                    else:
                        l = l.rstrip("\r\n").split("\t")
                        if "," in l[gnomadAltFreq_popmax_index]:
                            l[gnomadAltFreq_popmax_index] = l[
                                gnomadAltFreq_popmax_index
                            ].replace(",", ".")

                        if filter_mode:
                            if float(l[gnomadAltFreq_popmax_index]) < 0.005:
                                continue
                            if int(l[gnomadHomCount_all_index]) < 2:
                                continue

                        if l[0] not in already_seen_variants_dic:
                            already_seen_variants_dic[l[0]] = ""

                            if len(l) != split_header_length:
                                if len(l) == split_header_length - 1:
                                    l.append("")
                                else:
                                    log.error(
                                        f"Anormal number of tabulations on current lines compared to header. Exiting..."
                                    )
                                    log.error(f"{len(l)}")
                                    log.error(f"{l}")
                                    raise ValueError(len(l))
                            write_file.write(
                                "\t".join(l[i] for i in index_to_keep) + "\r\n"
                            )
    log.info(f"Sorting non redundant file...")
    subprocess.call(
        "head -n1 " + output_file + "_unsorted > " + output_file, shell=True
    )
    subprocess.call(
        "tail -n+2 " + output_file + "_unsorted | sort -t $'\t' -k 1 >> " + output_file,
        shell=True,
    )
    os.remove(output_file + "_unsorted")
    log.info(f"{output_file} generated !")


if __name__ == "__main__":
    pass
