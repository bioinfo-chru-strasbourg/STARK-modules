import glob
from operator import truediv
from os.path import join as osj
import os
import json
import logging as log
import subprocess


def generate(run_informations):
    # Adapted from the Samuel's non_redondant.py script
    input_dir = run_informations["tmp_analysis_folder"]
    columns_to_remove = []

    if "WES" in input_dir:
        filter_mode = True
        log.warning(
            "WES Folder detected, some variants will be filtered in the non-redundant : gnomadAltFreq_popmax_index < 0.005 & gnomadHomCount_all_index < 2"
        )
    else:
        filter_mode = False

    module_config = osj(
        os.environ["HOST_MODULE_CONFIG"],
        f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
    )
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        columns_to_remove = data["non_redundant_unwanted_columns"][
            run_informations["run_platform_application"]
        ]

    input_files_list = glob.glob(
        osj(
            input_dir,
            "VANNOT_*.tsv",
        )
    )
    input_files_dict = {}
    if run_informations["type"] == "run":
        for input_file in input_files_list:
            if "transcripts" in input_file:
                continue
            elif input_file.split(".")[2] != "tsv":
                panel = input_file.split(".")[2]
                if panel in input_files_dict.keys():
                    previous_value = input_files_dict[panel]
                    previous_value.append(input_file)
                    input_files_dict[panel] = previous_value
                else:
                    input_files_dict[panel] = [input_file]
            else:
                if "design" in input_files_dict.keys():
                    previous_value = input_files_dict["design"]
                    previous_value.append(input_file)
                    input_files_dict["design"] = previous_value
                else:
                    input_files_dict["design"] = [input_file]
    else:
        for input_file in input_files_list:
            if "transcripts" in input_file:
                continue
            if "folder" in input_files_dict.keys():
                previous_value = input_files_dict["folder"]
                previous_value.append(input_file)
                input_files_dict["folder"] = previous_value
            else:
                input_files_dict["folder"] = [input_file]

    samples = []
    for i in input_files_list:
        samples.append(os.path.basename(i).split(".")[0][7:])

    for key, input_files_list in input_files_dict.items():
        index_to_keep = []
        already_seen_variants_dic = {}
        first_file = True
        cleaned_header = []
        output_file = osj(input_dir, f"Non_Redondant.{key}.tsv")

        with open(output_file + "_unsorted", "w") as write_file:
            for i, f in enumerate(input_files_list):
                with open(f, "r") as read_file:
                    for l in read_file:
                        if l.startswith("##"):
                            continue
                        elif l in ["\n", "\r\n"]:
                            continue
                        elif l.startswith("#CHROM"):
                            if first_file:
                                l = l.rstrip("\r\n").split("\t")
                                split_header_length = len(l)
                                for i in range(len(l)):
                                    if (
                                        l[i] not in columns_to_remove
                                        and l[i] not in samples
                                    ):
                                        index_to_keep.append(i)
                                    if l[i] not in samples:
                                        cleaned_header.append(l[i])
                                header = "\t".join(cleaned_header)

                                write_file.write(
                                    "\t".join(l[i] for i in index_to_keep) + "\r\n"
                                )
                                if (
                                    "gnomadAltFreq_popmax" in l
                                    and "gnomadHomCount_all" in l
                                ):
                                    gnomadAltFreq_popmax_index = l.index(
                                        "gnomadAltFreq_popmax"
                                    )
                                    gnomadHomCount_all_index = l.index(
                                        "gnomadHomCount_all"
                                    )
                                first_file = False
                            else:
                                cleaned_line = []
                                l = l.rstrip("\r\n").split("\t")
                                for i in range(len(l)):
                                    if l[i] not in samples:
                                        cleaned_line.append(l[i])
                                l = "\t".join(cleaned_line)
                                assert l == header
                        else:
                            l = l.rstrip("\r\n").split("\t")
                            if (
                                "gnomadAltFreq_popmax" in l
                                and "gnomadHomCount_all" in l
                            ):
                                if "," in l[gnomadAltFreq_popmax_index]:
                                    l[gnomadAltFreq_popmax_index] = l[
                                        gnomadAltFreq_popmax_index
                                    ].replace(",", ".")

                                if filter_mode:
                                    if float(l[gnomadAltFreq_popmax_index]) < 0.005:
                                        continue
                                    if int(l[gnomadHomCount_all_index]) < 2:
                                        continue

                                if "." in l[gnomadAltFreq_popmax_index]:
                                    l[gnomadAltFreq_popmax_index] = l[
                                        gnomadAltFreq_popmax_index
                                    ].replace(".", ",")

                            if (
                                f"{l[0]}_{l[1]}_{l[2]}_{l[3]}"
                                not in already_seen_variants_dic
                            ):
                                already_seen_variants_dic[
                                    f"{l[0]}_{l[1]}_{l[2]}_{l[3]}"
                                ] = ""

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
            "tail -n+2 "
            + output_file
            + "_unsorted | sort -t $'\t' -k 1 >> "
            + output_file,
            shell=True,
        )
        os.remove(output_file + "_unsorted")
        log.info(f"{output_file} generated !")


if __name__ == "__main__":
    pass
