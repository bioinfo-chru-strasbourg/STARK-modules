import glob
from os.path import join as osj
import os
import subprocess
import logging as log
import shutil
import json
import re
import time
import gzip
from multiprocessing import Pool

from vannotplus.family.barcode import main_barcode_fast
from vannotplus.exomiser.exomiser import main_exomiser
from vannotplus.annot.score import main_annot
from vannotplus.__main__ import load_config, main_config

import howard_launcher
import commons

def project_folder_initialisation(run_informations):
    print("sam:project_folder_initialisation")
    vcf_file_list = glob.glob(osj(run_informations["archives_run_folder"], "VCF", "*", "*.vcf*"))
    if os.path.isdir(run_informations["tmp_analysis_folder"]):
        log.info("Cleaning temporary analysis folder")
        shutil.rmtree(run_informations["tmp_analysis_folder"])
        os.mkdir(run_informations["tmp_analysis_folder"])
    else:
        os.mkdir(run_informations["tmp_analysis_folder"])

    log.info("Copying vcf files to temporary folder")
    for processed_vcf_file in vcf_file_list:
        subprocess.run(
            [
                "rsync",
                "-rp",
                processed_vcf_file,
                osj(run_informations["tmp_analysis_folder"], ""),
            ]
        )

    vcf_file_to_analyse = glob.glob(
        osj(run_informations["tmp_analysis_folder"], "*vcf*")
    )

    for vcf_file in vcf_file_to_analyse:
        unmerge_vcf(vcf_file, run_informations)
        unmerged_vcfs = glob.glob(osj(run_informations["tmp_analysis_folder"], "unmerged_*vcf*"))
        if len(unmerged_vcfs) > 1:
            for unmerged_vcf in unmerged_vcfs:
                info_to_format_script(unmerged_vcf, run_informations)
            merge_vcf(run_informations, "0", os.path.basename(vcf_file))
        else:
            info_to_format_script(vcf_file, run_informations)

        cleaned_vcf = cleaning_annotations(vcf_file, run_informations)
        
        sample_list = subprocess.run(["bcftools", "query", "-l", cleaned_vcf],universal_newlines=True,stdout=subprocess.PIPE,).stdout.strip().split("\n")

        if run_informations["run_platform_application"] != None and len(sample_list) >= 1:
            output_exomiser = osj(
                run_informations["tmp_analysis_folder"],
                "exomized_" + os.path.basename(cleaned_vcf),
            )
            vannotplus_config = osj(os.environ["HOST_MODULE_CONFIG"], "vannotplus.yml")
            main_exomiser(
                cleaned_vcf,
                output_exomiser,
                run_informations["run_application"],
                load_config(vannotplus_config),
            )
            os.remove(cleaned_vcf)
            os.rename(output_exomiser, vcf_file)

def folder_initialisation(run_informations):
    vcf_file_list = glob.glob(osj(run_informations["archives_run_folder"], "*.vcf*"))
    if os.path.isdir(run_informations["tmp_analysis_folder"]):
        log.info("Cleaning temporary analysis folder")
        shutil.rmtree(run_informations["tmp_analysis_folder"])
        os.mkdir(run_informations["tmp_analysis_folder"])
    else:
        os.mkdir(run_informations["tmp_analysis_folder"])

    log.info("Copying vcf files to temporary folder")
    for processed_vcf_file in vcf_file_list:
        subprocess.run(
            [
                "rsync",
                "-rp",
                processed_vcf_file,
                osj(run_informations["tmp_analysis_folder"], ""),
            ]
        )

    vcf_file_to_analyse = glob.glob(
        osj(run_informations["tmp_analysis_folder"], "*vcf*")
    )
    for vcf_file in vcf_file_to_analyse:
        # unmerge_vcf(vcf_file, run_informations)
        # unmerged_vcfs = glob.glob(osj(run_informations["tmp_analysis_folder"], "unmerged_*vcf*"))
        # if len(unmerged_vcfs) > 1:
        #     for unmerged_vcf in unmerged_vcfs:
        #         info_to_format_script(unmerged_vcf, run_informations)
        #     vcf_file = merge_vcf(run_informations, "0", os.path.basename(vcf_file))
        # else:
        #     info_to_format_script(vcf_file, run_informations)

        cleaned_vcf = cleaning_annotations(vcf_file, run_informations)

        sample_list = subprocess.run(["bcftools", "query", "-l", cleaned_vcf],universal_newlines=True,stdout=subprocess.PIPE,).stdout.strip().split("\n")

        if run_informations["run_platform_application"] != None and len(sample_list) >= 1:
            output_exomiser = osj(
                run_informations["tmp_analysis_folder"],
                "exomized_" + os.path.basename(cleaned_vcf),
            )
            vannotplus_config = osj(os.environ["HOST_MODULE_CONFIG"], "vannotplus.yml")
            main_exomiser(
                cleaned_vcf,
                output_exomiser,
                run_informations["run_application"],
                load_config(vannotplus_config),
            )
            os.remove(cleaned_vcf)
            os.rename(output_exomiser, vcf_file)


def run_initialisation(run_informations):
    vcf_file_list = glob.glob(osj(run_informations["archives_run_folder"], "*.vcf*"))
    if os.path.isdir(run_informations["tmp_analysis_folder"]):
        log.info("Cleaning temporary analysis folder")
        shutil.rmtree(run_informations["tmp_analysis_folder"])
        os.mkdir(run_informations["tmp_analysis_folder"])
    else:
        os.mkdir(osj(run_informations["tmp_analysis_folder"]))

    log.info("Copying vcf files to temporary folder")
    for processed_vcf_file in vcf_file_list:
        subprocess.run(
            [
                "rsync",
                "-rp",
                processed_vcf_file,
                osj(run_informations["tmp_analysis_folder"], ""),
            ]
        )

    vcf_file_to_analyse = glob.glob(
        osj(run_informations["tmp_analysis_folder"], "*.vcf*")
    )

    for vcf_file in vcf_file_to_analyse:
        info_to_format_script(vcf_file, run_informations)
        cleaned_vcf = cleaning_annotations(vcf_file, run_informations)
        
        sample_list = subprocess.run(["bcftools", "query", "-l", cleaned_vcf],universal_newlines=True,stdout=subprocess.PIPE,).stdout.strip().split("\n")
        print(sample_list)
        print(cleaned_vcf)
        if run_informations["run_platform_application"] != None and len(sample_list) >= 1:
            output_exomiser = osj(
                run_informations["tmp_analysis_folder"],
                "exomized_" + os.path.basename(cleaned_vcf),
            )
            vannotplus_config = osj(os.environ["HOST_MODULE_CONFIG"], "vannotplus.yml")
            main_exomiser(
                cleaned_vcf,
                output_exomiser,
                run_informations["run_application"],
                load_config(vannotplus_config),
            )
            os.remove(cleaned_vcf)
            os.rename(output_exomiser, vcf_file)

def cleaning_annotations(vcf_file, run_informations):
    module_config = osj(
        os.environ["HOST_MODULE_CONFIG"],
        f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
    )
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        annotations_to_keep = data["keep_vcf_info"][
            run_informations["run_platform_application"]
        ]
    log.info("Cleaning INFO column in the provided vcfs")
    log.info(f"Kept informations : {", ".join(annotations_to_keep)}")

    actual_info_fields = subprocess.run(
        ["zgrep", "##INFO", vcf_file], capture_output=True, text=True
    )
    actual_info_fields = actual_info_fields.stdout.strip().split("##")
    if not vcf_file.endswith(".vcf.gz"):
        subprocess.call(["bgzip", vcf_file])
        vcf_file = vcf_file + ".gz"

    cleaned_vcf = osj(
        os.path.dirname(vcf_file), "cleaned_" + os.path.basename(vcf_file)[:-3]
    )
    info_to_keep = []
    for i in annotations_to_keep:
        for j in actual_info_fields:
            if re.search(i, j):
                info_to_keep.append("INFO/" + j.split(",")[0].split("=")[-1])

    if len(info_to_keep) == 0:
        log.info(
            f"No annotations to keep were found in {os.path.basename(vcf_file)}, deleting all annotations"
        )
        info_to_keep = "INFO"
    else:
        log.info(
            f"Keeping following annotations: {' '.join(info_to_keep)} for sample {os.path.basename(vcf_file)}"
        )
        info_to_keep = "^" + ",".join(info_to_keep)

    cmd = ["bcftools", "annotate", "-x"]
    cmd.append(info_to_keep)
    print(vcf_file)
    cmd.append(vcf_file)
    print(" ".join(cmd))
    with open(cleaned_vcf, "w") as output:
        subprocess.call(cmd, stdout=output, universal_newlines=True)
    os.remove(vcf_file)

    subprocess.call(["bgzip", cleaned_vcf], universal_newlines=True)
    renamed_clean = osj(
        os.path.dirname(cleaned_vcf),
        os.path.basename(cleaned_vcf).replace("cleaned_", "") + ".gz",
    )
    cleaned_vcf = cleaned_vcf + ".gz"
    os.rename(cleaned_vcf, renamed_clean)
    return renamed_clean


def fambarcode_vcf(run_informations, input_vcf):
    output = osj(
        run_informations["tmp_analysis_folder"],
        "fambarcode_" + os.path.basename(input_vcf),
    )
    module_config = osj(
        os.environ["HOST_MODULE_CONFIG"],
        f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
    )
    howard_config_container = osj(
        os.environ["HOST_MODULE_CONFIG"], "howard", "howard_config.json"
    )
    howard_config_host = osj(os.environ["HOST_CONFIG"], "howard", "howard_config.json")

    if not os.path.isfile(module_config):
        log.error(f"{module_config} do not exist, primordial file, check its existence")
        raise ValueError(module_config)
    elif not os.path.isfile(howard_config_container):
        log.error(
            f"{howard_config_container} do not exist, primordial file, check its existence"
        )
        raise ValueError(howard_config_container)

    vannotplus_config = osj(os.environ["HOST_MODULE_CONFIG"], "vannotplus.yml")
    fambarcode_config = load_config(vannotplus_config)

    exact_time = time.time() + 7200
    local_time = time.localtime(exact_time)
    actual_time = time.strftime("%H%M%S", local_time)
    start = actual_time
    container_name = f"VANNOT_fambarcode_{start}_{run_informations['run_name']}"

    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        howard_image = data["howard_image"]
    log.info(f"Using {howard_image}")

    howard_bin = f"docker run --rm --name {container_name} --env http_proxy={os.environ["http_proxy"]} --env ftp_proxy={os.environ["ftp_proxy"]} --env https_proxy={os.environ["https_proxy"]} -v /tmp/:/tmp/ -v {os.environ["HOST_TMP"]}:{os.environ["HOST_TMP"]} -v {os.environ["HOST_DATABASES"]}:/databases/ -v {os.environ["HOST_SERVICES"]}:{os.environ["HOST_SERVICES"]} -v {os.environ["HOST_CONFIG"]}:{os.environ["HOST_CONFIG"]} -v {howard_config_host}:/tools/howard/current/config/config.json -v /var/run/docker.sock:/var/run/docker.sock {howard_image} calculation"
    fambarcode_config["howard"]["bin"] = howard_bin
    main_barcode_fast(
        input_vcf,
        output,
        run_informations["run_application"],
        fambarcode_config,
    )

    os.rename(
        output,
        osj(os.path.dirname(output), "_".join(os.path.basename(output).split("_")[1:])),
    )
    output = osj(
        os.path.dirname(output), "_".join(os.path.basename(output).split("_")[1:])
    )
    print("sam:fambarcode_vcf output", output)
    return output

    # tmp_output = osj(os.path.dirname(output), "tmp_" + os.path.basename(output)[:-3])
    # tmp_header = osj(os.path.dirname(output), "tmp_header_" + os.path.basename(output)[:-3])

    # with gzip.open(output, "rt") as read_file:
    #     with open(tmp_output, "w") as write_file:
    #         with open(tmp_header, "w") as write_header:
    #             lines = read_file.readlines()
    #             for line in lines:
    #                 if line.startswith("##"):
    #                     write_header.write(line)
    #                 elif line.startswith("#CHROM"):
    #                     write_file.write(line)
    #                     line = line.rstrip("\n").split("\t")
    #                     if "FORMAT" in line:
    #                         format_index = line.index("FORMAT")
    #                 else:
    #                     line = line.rstrip("\n").split("\t")
    #                     format_values = line[format_index].split(":")
    #                     if format_values.count("BCF") > 1:
    #                         bcf_wanted = ["BCF", "BCFS"]
    #                         bcf_original = ["BCF", "BCFS"]
    #                         bcf_count = format_values.count("BCF")
    #                         for i in range(1, bcf_count):
    #                             bcf_wanted.append(f"BCF_{i}")
    #                             bcf_wanted.append(f"BCFS_{i}")
    #                             bcf_original.append("BCF")
    #                             bcf_original.append("BCFS")
    #                         bcf_wanted = ":".join(bcf_wanted)
    #                         bcf_original = ":".join(bcf_original)
    #                         format_values = ":".join(format_values)
    #                         line[format_index] = format_values.replace(bcf_original, bcf_wanted)
    #                     else:
    #                         bcf_count = 0
    #                     write_file.write("\t".join(line) + "\n")

    # os.remove(output)
    # output = output[:-3]

    # with open(tmp_header, "r") as read_file:
    #     with open(tmp_output, "r") as read_file2:
    #         with open(output, "w") as write_file:
    #             lines = read_file.readlines()
    #             for line in lines:
    #                 if bcf_count > 1:
    #                     if line.startswith("##FORMAT=<ID=BCF,"):
    #                         write_file.write(line)
    #                         for i in range(1, bcf_count):
    #                             write_file.write(f"##FORMAT=<ID=BCF_{i},Number=.,Type=String,Description=\"barcode family calculation\">\n")
    #                     elif line.startswith("##FORMAT=<ID=BCFS,"):
    #                         write_file.write(line)
    #                         for i in range(1, bcf_count):
    #                             write_file.write(f"##FORMAT=<ID=BCFS_{i},Number=.,Type=String,Description=\"barcode family samples\">\n")
    #                     else:
    #                         write_file.write(line)
    #                 else:
    #                     write_file.write(line)
    #             lines = read_file2.readlines()
    #             for line in lines:
    #                 write_file.write(line)
    
    # os.remove(tmp_output)
    # os.remove(tmp_header)
    # subprocess.call("bgzip " + output, shell=True)  
    # output = output + ".gz" 
    # return output

def merge_vcf(run_informations, step, base_vcf):
    if step == "0":
        vcf_file_to_merge = glob.glob(osj(run_informations["tmp_analysis_folder"], "fixed_unmerged_*.vcf*"))
        if len(vcf_file_to_merge) == 0:
            vcf_file_to_merge = glob.glob(osj(run_informations["tmp_analysis_folder"], "unmerged_*.vcf*"))
    else:
        print("sam: running merge_vcf step=1")
        print("run_informations['tmp_analysis_folder']", run_informations["tmp_analysis_folder"])
        vcf_file_to_merge = glob.glob(
            osj(run_informations["tmp_analysis_folder"], "*.vcf*")
        )
    if len(vcf_file_to_merge) > 1:
        for vcf_file in vcf_file_to_merge:
            tmp_output = osj(os.path.dirname(vcf_file), "tmp_" + os.path.basename(vcf_file))
            with gzip.open(vcf_file, "rt") as read_file:
                with open(tmp_output, "w") as write_file:
                    lines = read_file.readlines()
                    for line in lines:
                        if line.startswith("##FORMAT=<ID=VAF,Number="):
                            log.info(f"Fixing VAF for sample {vcf_file}")
                            write_file.write("##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"VAF Variant Frequency, calculated from quality\">\n")
                        elif line.startswith("##FORMAT=<ID=AD,Number="):
                            log.info(f"Fixing AD for sample {vcf_file}")
                            write_file.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n")
                        elif line.startswith("##FORMAT=<ID=FT,Number="):
                            log.info(f"Fixing FT for sample {vcf_file}")
                            write_file.write("##FORMAT=<ID=FT,Number=.,Type=String,Description=\"Genotype-level filter\">\n")
                        elif line.startswith("##reference=file:"):
                            log.info(f"Fixing reference for sample {vcf_file}")
                            write_file.write("##reference=file:///STARK/databases/genomes/current/hg19.fa\n")
                        elif line.startswith("##FORMAT=<ID=PL,Number="):
                            log.info(f"Fixing PL for sample {vcf_file}")
                            write_file.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n")
                        else:
                            write_file.write(line)
            os.remove(vcf_file)
            subprocess.call("bgzip " + tmp_output, shell=True)
            os.rename(tmp_output + ".gz", vcf_file)
        
        log.info(f"Merging {len(vcf_file_to_merge)} vcf files")
        for i in vcf_file_to_merge:
            subprocess.call(["tabix", i], universal_newlines=True)

        output_merged = osj(
            run_informations["tmp_analysis_folder"],
            f"VANNOT_merged_{run_informations["run_name"]}.design.vcf.gz",
        )
        cmd = ["bcftools", "merge"] + vcf_file_to_merge
        cmd_args = ["-m", "none", "-O", "z", "-o", output_merged]
        cmd = cmd + cmd_args
        log.debug(" ".join(cmd))
        subprocess.call(cmd, universal_newlines=True)
        if step == "1":
            for i in vcf_file_to_merge:
                os.remove(i)
                os.remove(i + ".tbi")
            os.rename(output_merged, osj(os.path.dirname(output_merged), os.path.basename(output_merged).removeprefix("VANNOT_merged_")))
            return osj(os.path.dirname(output_merged), os.path.basename(output_merged).removeprefix("VANNOT_merged_"))
        if step == "0":
            for i in vcf_file_to_merge:
                os.remove(i)
                os.remove(i + ".tbi")
            if base_vcf != "":
                os.rename(output_merged, osj(os.path.dirname(output_merged), base_vcf))
            else:
                os.rename(output_merged, osj(os.path.dirname(output_merged), os.path.basename(output_merged).removeprefix("VANNOT_").replace(".design", "")))
            return osj(os.path.dirname(output_merged), os.path.basename(output_merged).removeprefix("VANNOT_").replace(".design", ""))
        else:
            return output_merged
    elif len(vcf_file_to_merge) == 0:
        log.error("No vcf files to merge after transcript score calculation")
        raise ValueError(vcf_file_to_merge)
    else:
        return vcf_file_to_merge[0]


def unmerge_vcf(input, run_informations):
    vcf_file_to_unmerge = input
    print("sam: unmerge_vcf", vcf_file_to_unmerge)
    sample_list = (
        subprocess.run(
            ["bcftools", "query", "-l", vcf_file_to_unmerge],
            universal_newlines=True,
            stdout=subprocess.PIPE,
        )
        .stdout.strip()
        .split("\n")
    )
    if len(sample_list) > 1:
        for sample in sample_list:
            if run_informations["type"] == "run":
                output_file = osj(
                    os.path.dirname(vcf_file_to_unmerge), f"unmerged_VANNOT_{sample}.design.vcf"
                )
            elif run_informations["type"] == "dejavu":
                output_file = osj(os.path.dirname(vcf_file_to_unmerge), f"{sample}.vcf")
            else:
                output_file = osj(
                    os.path.dirname(vcf_file_to_unmerge), f"unmerged_VANNOT_{sample}.vcf"
                )
            cmd = ["bcftools", "view", "-c1", "-s", sample, vcf_file_to_unmerge]
            with open(output_file, "w") as writefile:
                subprocess.call(cmd, universal_newlines=True, stdout=writefile)
            subprocess.call(["bgzip", output_file], universal_newlines=True)
            print("sam:unmerge_vcf output_file:", output_file)
        os.remove(vcf_file_to_unmerge)


def info_to_format_script(vcf_file, run_informations):
    log.info(
        f"Moving desired columns to FORMAT column for {os.path.basename(vcf_file)}"
    )
    module_config = osj(
        os.environ["HOST_MODULE_CONFIG"],
        f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
    )
    output_file = osj(
        os.path.dirname(vcf_file), f"fixed_{os.path.basename(vcf_file)[:-3]}"
    )
    fixed_file = osj(
        os.path.dirname(vcf_file), f"fixed2_{os.path.basename(vcf_file)[:-3]}"
    )
    is_pool = False
    with gzip.open(vcf_file, "rt") as readfile:
        with open(fixed_file, "wt") as writefile:
            lines = readfile.readlines()
            for line in lines:
                if line.startswith("##INFO=<ID=POOL_F_Depth"):
                    is_pool = True
                else:
                    continue
            for line in lines:
                if line.startswith("##INFO=<ID=BARCODE") and is_pool == True:
                    line = line.replace("##INFO=<ID=BARCODE", "##INFO=<ID=POOL_BARCODE")
                    writefile.write(line)
                elif line.startswith("##"):
                    writefile.write(line)
                elif line.startswith("#CHROM"):
                    line = line.rstrip("\n").split("\t")
                    sample = line[-1]
                    writefile.write("\t".join(line) + "\n")
                else:
                    line = line.rstrip("\n").split("\t")
                    info_line = line[7].split(";")
                    info_line_check = line[7].replace("=",";").split(";")
                    if "POOL_F_Depth" in info_line_check and "POOL_M_Depth" in info_line_check:
                        for i,info_value in enumerate(info_line):
                            if info_value.startswith("BARCODE="):
                                info_value = info_value.replace("BARCODE=", "POOL_BARCODE=")
                            info_line[i] = info_value
                    line[7] = ";".join(info_line)
                    writefile.write("\t".join(line) + "\n")
    
    subprocess.call(["bgzip", fixed_file], universal_newlines=True)    
    os.rename(fixed_file + ".gz", vcf_file)
    
    sample = subprocess.run(["bcftools", "query", "-l", vcf_file],universal_newlines=True,stdout=subprocess.PIPE,).stdout.strip().split("\n")[0]
    tmp_annot = osj(os.path.dirname(vcf_file), "annot.txt.tmp")
    tmp_annot_fixed = osj(os.path.dirname(vcf_file), "annot.fixed.txt.tmp")
    tmp_hdr = osj(os.path.dirname(vcf_file), "hdr.txt.tmp")
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        info_to_format_columns = data["info_to_format"][
            run_informations["run_platform_application"]
        ]

    info_to_format_columns_query = "\\t%".join(info_to_format_columns)
    info_to_format_columns_query = (
        "%CHROM\\t%POS\\t%REF\\t%ALT\\t%" + info_to_format_columns_query + "\n"
    )

    cmd = ["bcftools", "query", "-f", info_to_format_columns_query, vcf_file]

    with open(tmp_annot, "a") as writefile:
        subprocess.call(cmd, universal_newlines=True, stdout=writefile)

    with open(tmp_annot, "r") as readfile:
        with open(tmp_annot_fixed, "a") as writefile:
            lines = readfile.readlines()
            for line in lines:
                line = line.strip()
                writefile.write(line.replace(":", ",") + "\n")
    subprocess.call(["bgzip", tmp_annot_fixed], universal_newlines=True)
    tmp_annot_fixed = tmp_annot_fixed + ".gz"
    cmd = ["tabix", "-s1", "-b2", "-e2", tmp_annot_fixed]
    subprocess.call(cmd, universal_newlines=True)
    vcf_file_gunzip = vcf_file[:-3]
    subprocess.call(["gunzip", vcf_file])
    with open(vcf_file_gunzip, "r") as readfile:
        with open(tmp_hdr, "a") as writefile:
            lines = readfile.readlines()
            for line in lines:
                for column in info_to_format_columns:
                    if line.startswith("##INFO=<ID=" + column):
                        writefile.write(line.replace("INFO", "FORMAT"))
    subprocess.call(["bgzip", vcf_file_gunzip], universal_newlines=True)
    info_to_format_columns_annotate = ",FORMAT/".join(info_to_format_columns)
    info_to_format_columns_annotate = (
        "CHROM,POS,REF,ALT,FORMAT/" + info_to_format_columns_annotate
    )

    cmd = [
        "bcftools",
        "annotate",
        "-s",
        sample,
        "-a",
        tmp_annot_fixed,
        "-h",
        tmp_hdr,
        "-c",
        info_to_format_columns_annotate,
        vcf_file,
    ]
    log.debug(" ".join(cmd))
    with open(output_file, "a") as writefile:
        subprocess.call(cmd, universal_newlines=True, stdout=writefile)
    shutil.copy(output_file, output_file + ".unzipped")
    unzipped_output_file = output_file + ".unzipped"
    subprocess.call(["bgzip", output_file], universal_newlines=True)
    output_file = output_file + ".gz"

    os.remove(tmp_annot)
    os.remove(tmp_annot_fixed)
    os.remove(tmp_annot_fixed + ".tbi")
    os.remove(tmp_hdr)
    if os.stat(unzipped_output_file).st_size == 0:
        os.remove(unzipped_output_file)
        os.remove(output_file)
        return vcf_file
    else:
        os.remove(unzipped_output_file)
        os.remove(vcf_file)
        os.rename(output_file, vcf_file)
        return vcf_file


def howard_proc(run_informations, vcf_file):
    log.info(f"Launching HOWARD analysis for {vcf_file}")
    if run_informations["type"] == "run":
        if run_informations["output_format"] != None:
            output_file = osj(
                run_informations["tmp_analysis_folder"],
                f"VANNOT_{os.path.basename(vcf_file).split(".")[0]}.design.{run_informations["output_format"]}",
            )
        else:
            output_file = osj(
                run_informations["tmp_analysis_folder"],
                f"VANNOT_{os.path.basename(vcf_file).split(".")[0]}.design.{".".join(os.path.basename(vcf_file).split(".")[1:])}",
            )
    elif run_informations["type"] == "folder":
        if run_informations["output_format"] != None:
            output_file = osj(
                run_informations["tmp_analysis_folder"],
                f"VANNOT_{os.path.basename(vcf_file).split(".")[0]}.{run_informations["output_format"]}",
            )
        else:
            output_file = osj(
                run_informations["tmp_analysis_folder"],
                f"VANNOT_{os.path.basename(vcf_file).split(".")[0]}.{".".join(os.path.basename(vcf_file).split(".")[1:])}",
            )
    if run_informations["parameters_file"] == None:
        for option in [
            run_informations["run_platform_application"],
            run_informations["run_platform"],
            "default",
            "none",
        ]:
            configfile = osj(
                os.environ["HOST_MODULE_CONFIG"],
                "paramfiles",
                f"param.{option}.json",
            )
            if os.path.isfile(configfile):
                break
            elif configfile == osj(
                os.environ["HOST_MODULE_CONFIG"], "paramfiles", "param.none.json"
            ):
                log.error(
                    "Missing parameter file for your analysis, didn't find application nor platform not default parameters"
                )
                raise ValueError(configfile)
            else:
                continue
    else:
        configfile = run_informations["parameters_file"]

    log.info(f"Using {configfile} as parameter for HOWARD analysis")
    if not os.path.isfile(configfile):
        log.error("param.default.json not found, please check your config directory")
        raise ValueError(configfile)

    howard_config = osj(
        os.environ["HOST_MODULE_CONFIG"], "howard", "howard_config.json"
    )
    threads = commons.get_threads("threads_annotation")
    memory = commons.get_memory("memory_annotation")

    exact_time = time.time() + 7200
    local_time = time.localtime(exact_time)
    actual_time = time.strftime("%H%M%S", local_time)
    start = actual_time

    container_name = f"VANNOT_annotate_{start}_{run_informations['run_name']}_{os.path.basename(vcf_file).split('.')[0]}"
    launch_annotate_arguments = [
        "process",
        "--input",
        vcf_file,
        "--output",
        output_file,
        "--param",
        configfile,
        "--assembly",
        run_informations["assembly"],
        "--memory",
        memory,
        "--threads",
        threads,
        "--config",
        howard_config,
    ]

    log.info("Annotating input files with HOWARD")
    howard_launcher.launch(container_name, launch_annotate_arguments)
    os.remove(vcf_file)

    return output_file


def gmc_score(run_informations):
    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))

    for vcf_file in vcf_files:
        output_gmc = osj(
            run_informations["tmp_analysis_folder"],
            "gmc_" + os.path.basename(vcf_file),
        )
        vannotplus_config = osj(os.environ["HOST_MODULE_CONFIG"], "vannotplus.yml")
        main_annot(
            vcf_file,
            output_gmc,
            load_config(vannotplus_config),
        )
        os.remove(vcf_file)
        os.rename(output_gmc, vcf_file)
        print("sam:gmc_score output", vcf_file)

# def divide_memory(threads: str, memory: str) -> str:
#     suffix = memory[-1]
#     memory = int(memory[:-1])
#     memory = memory/int(threads)
#     return str(memory) + suffix

# def prioritize_worker(vcf_file: str, run_informations: dict, memory: str, start: str, transcript_param: str, threads: str) -> None:



# def prioritize_worker_unpacker(args):
#     vcf_file, run_informations, memory, start, transcript_param, threads = args
#     prioritize_worker(vcf_file, run_informations, memory, start, transcript_param, threads)
#Fix parallelization prio
def howard_score_transcripts(run_informations):
    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))
    print(vcf_files)
    howard_config = osj(os.environ["HOST_CONFIG"], "howard", "howard_config.json")


    if len(vcf_files) > 1:
        for vcf_file in vcf_files:
            os.rename(vcf_file, osj(os.path.dirname(vcf_file), os.path.basename(vcf_file).removeprefix("unmerged_")))
        vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))

    transcript_param = osj(
        os.environ["HOST_MODULE_CONFIG"],
        "howard",
        "param.transcripts.json",
    )

    threads = commons.get_threads("threads_annotation")
    memory = commons.get_memory("memory_annotation")

    exact_time = time.time() + 7200
    local_time = time.localtime(exact_time)
    actual_time = time.strftime("%H%M%S", local_time)
    start = actual_time

#     with Pool(1) as pool:
#         pool.map(prioritize_worker_unpacker, [(vcf_file, run_informations, divide_memory(threads, memory), start, transcript_param, threads
# ) for vcf_file in vcf_files])
    for vcf_file in vcf_files:
        output_file_transcripts = osj(run_informations["tmp_analysis_folder"], f"{os.path.basename(vcf_file).split(".")[0]}_output_transcripts.vcf.gz")

        container_name = f"VANNOT_transcripts_{start}_{run_informations['run_name']}_{os.path.basename(vcf_file).split('.')[0]}"
        launch_annotate_arguments = [
            "process",
            "--input",
            vcf_file,
            "--output",
            output_file_transcripts,
            "--param",
            transcript_param,
            "--config",
            howard_config,
            "--memory",
            memory,
            "--threads",
            threads,
            "--debug"
        ]

        log.info("Prioritization of transcripts")
        howard_launcher.launch(container_name, launch_annotate_arguments)
        os.remove(vcf_file)
        os.rename(output_file_transcripts, vcf_file)

        with open(
            osj(os.environ["HOST_MODULE_CONFIG"], "howard", "param.transcripts.json"),
            "r",
        ) as read_file:
            data = json.load(read_file)
            transcripts_output = data["transcripts"]["export"]["output"]

        sample_name = (os.path.basename(vcf_file).split(".")[0]).removeprefix("VANNOT_")
        transcripts_output_renamed = f"VANNOT_transcripts_{sample_name}.tsv"
        shutil.copy(
            transcripts_output,
            osj(run_informations["tmp_analysis_folder"], transcripts_output_renamed),
        )
        # shutil.copy(
        #     transcripts_output,
        #     "/home1/data/STARK/data/samtranscripts/",
        # )
        # print("sam: howard_score_transcripts output", osj(run_informations["tmp_analysis_folder"], transcripts_output_renamed))


def convert_to_final_tsv(run_informations):
    log.info("Converting output file into readable tsv")
    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))
    threads = commons.get_threads("threads_conversion")
    memory = commons.get_memory("memory_conversion")
    module_config = osj(
        os.environ["HOST_MODULE_CONFIG"],
        f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
    )
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        ordered_fields = data["vcf_to_tsv_column_order"][
            run_informations["run_platform_application"]
        ]

    explode_infos_fields = ",".join(ordered_fields)
    explode_infos_fields = explode_infos_fields + ",*"
    print(explode_infos_fields)
    howard_config = osj(
        os.environ["HOST_MODULE_CONFIG"], "howard", "howard_config.json"
    )

    for vcf_file in vcf_files:
        
        panel_name = vcf_file.split(".")[-3]
        exact_time = time.time() + 7200
        local_time = time.localtime(exact_time)
        actual_time = time.strftime("%H%M%S", local_time)
        start = actual_time
        container_name = f"VANNOT_convert_{start}_{run_informations['run_name']}_{os.path.basename(vcf_file).split('.')[0]}"
        if run_informations["type"] == "run":
            if panel_name != "design":
                output_file = osj(
                    run_informations["tmp_analysis_folder"],
                    f"{os.path.basename(vcf_file).split(".")[0]}.panel.{panel_name}.tsv",
                )
            else:
                output_file = osj(
                    run_informations["tmp_analysis_folder"],
                    f"{os.path.basename(vcf_file).split(".")[0]}.design.tsv",
                )
        else:
            output_file = osj(
                run_informations["tmp_analysis_folder"],
                f"{os.path.basename(vcf_file).split(".")[0]}.tsv",
            )

        if "merged" not in vcf_file:
            launch_convert_arguments = [
                "query",
                "--input",
                vcf_file,
                "--output",
                output_file,
                "--explode_infos",
                "--explode_infos_fields",
                explode_infos_fields,
                "--query",
                "SELECT * FROM variants",
                "--threads",
                threads,
                "--memory",
                memory,
                "--config",
                howard_config,
            ]
            howard_launcher.launch(container_name, launch_convert_arguments)
            output_file = format_explode(vcf_file, output_file)
            tsv_modifier(output_file, run_informations)

def format_explode(vcf_file, tsv_file):
    output_file = osj(os.path.dirname(tsv_file), "tmp_" + os.path.basename(tsv_file))
    tsv_lines = []
    format_dict = {}
    index_dict = {}
    format_type = ""
    is_format = True

    with gzip.open(vcf_file, "rt") as vcf_read:
        for line in vcf_read:
            if line.startswith("#CHROM"):
                line = line.rstrip("\n").split("\t")
                for i in range(len(line)):
                    if line[i] == "FORMAT":
                        format_type = i
                        format_values = i+1
                if format_type == "":
                    is_format = False
            elif not line.startswith("#") and is_format is True:
                line = line.rstrip("\n").split("\t")
                variant_name = line[0] + "_" + line[1] + "_" + line[3] + "_" + line[4]
                format_dict[variant_name] = [line[format_type].split(":"), line[format_values].split(":")]
    if is_format is True:
        with open(tsv_file, "r") as tsv_read:
            for line in tsv_read:
                if line.startswith("#CHROM"):
                    header = line.rstrip("\n").split("\t")
                else:
                    tsv_lines.append(line)

        is_second_ad_alt = False
        is_varscan_ad = False
        for values in format_dict.values():
            if "AD" in values[0]:
                ad_index = values[0].index("AD")
                ad_value = values[1][ad_index]
                if ad_value.count(",") == 2:
                    is_second_ad_alt = True
                elif "," not in ad_value and ad_value != ".":
                    is_varscan_ad = True

        modified_header = []
        new_tsv_lines = []
        for column_name in header:
            if column_name == "#CHROM":
                modified_header.append("chr")
            elif column_name == "AD":
                if is_varscan_ad == False:
                    ad_index = header.index(column_name)
                    column_name = "AD_ref"
                    modified_header.append(column_name)
                    modified_header.append("AD_alt")
                    if is_second_ad_alt is True:
                        modified_header.append("AD_alt2")
                else:
                    modified_header.append(column_name)
            else:
                modified_header.append(column_name)

        for line in tsv_lines:
            modified_tsv_line = []
            line = line.rstrip("\n").split("\t")
            for count, content in enumerate(line):
                if count == ad_index:
                    modified_tsv_line.append(content)
                    modified_tsv_line.append("")
                    if content.count(",") == 2:
                        modified_tsv_line.append("")
                else:
                    modified_tsv_line.append(content)
            new_tsv_lines.append(modified_tsv_line)

        tsv_lines = new_tsv_lines

        vcf_header = set()
        for values in format_dict.values():
            vcf_header = vcf_header | set(values[0])
            
        vcf_header = list(vcf_header)
        vcf_header_cleaned = []
        for i in vcf_header:
            if i not in modified_header:
                vcf_header_cleaned.append(i)
        new_header = modified_header + vcf_header_cleaned
        
        with open(output_file, "w") as write_file:
            write_file.write("\t".join(new_header) + "\n")

        for i in new_header:
            index_dict[i] = new_header.index(i)

        for line in tsv_lines:
            for i in range(len(new_header)-len(line)):
                line.append("")
            chr_index = index_dict["chr"]
            pos_index = index_dict["POS"]
            alt_index = index_dict["ALT"]
            ref_index = index_dict["REF"]
            tsv_variant = line[chr_index] + "_" + line[pos_index] + "_" + line[ref_index] + "_" + line[alt_index]
            if tsv_variant in format_dict.keys():
                format_value_per_type = {}
                for count, format_type in enumerate(format_dict[tsv_variant][0]):
                    format_value_per_type[format_type] = format_dict[tsv_variant][1][count]
                for key, value in format_value_per_type.items():
                    ad_ref_index = index_dict["AD_ref"]
                    ad_alt_index = index_dict["AD_alt"]
                    if key == "AD":
                        if is_varscan_ad is False and value != ".":
                            ad_ref = value.split(",")[0]
                            line[ad_ref_index] = ad_ref
                            ad_alt = value.split(",")[1]
                            line[ad_alt_index] = ad_alt
                            if value.count(",") == 2:
                                ad_alt_bis_index = index_dict["AD_alt"]
                                ad_alt_bis = value.split(",")[2]
                                line[ad_alt_bis_index] = ad_alt_bis
                        else:
                            line[index_dict[key]] = value
                    else:
                        line[index_dict[key]] = value
                with open(output_file, "a") as write_file:
                    write_file.write("\t".join(line) + "\n")

        os.remove(tsv_file)
        os.rename(output_file, tsv_file)
        return tsv_file
    
    else:
        with open(tsv_file, "r") as tsv_read:
            for line in tsv_read:
                if line.startswith("#CHROM"):
                    header = line.rstrip("\n").split("\t")
                else:
                    tsv_lines.append(line.rstrip("\n"))

        modified_header = []
        for column_name in header:
            if column_name == "#CHROM":
                modified_header.append("chr")
            else:
                modified_header.append(column_name)      

        new_header = modified_header
        with open(output_file, "w") as write_file:
            write_file.write("\t".join(new_header) + "\n")
        for line in tsv_lines:
            with open(output_file, "a") as write_file:
                write_file.write(line + "\n")

        os.remove(tsv_file)
        os.rename(output_file, tsv_file)
        return tsv_file

def tsv_modifier(input_file, run_informations):
    sample = os.path.basename(input_file).removesuffix(".tsv").removeprefix("VANNOT_")
    if sample.endswith(".design"):
        sample = sample.removesuffix(".design")
    elif len(sample.split(".")) > 1 and sample.split(".")[1] == "panel":
        sample = sample.split(".")[0]
    output_file = osj(
        os.path.dirname(input_file), "tmp_" + os.path.basename(input_file)
    )
    module_config = osj(os.environ["HOST_MODULE_CONFIG"],f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        values_to_delete = data["tsv_columns_to_remove"]

    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        last_order = data["vcf_to_tsv_column_order"][f"{run_informations["run_platform_application"]}"][-1]

    values_to_delete.append(sample)
    dejavu_to_keep = []
    prefix_to_keep = [run_informations["run_application"], "WES_AGILENT", "WES_TWIST", "WES_ROCHE"]
    for i in prefix_to_keep:
        dejavu_to_keep.append(f"{i}_ALLELECOUNT")
        dejavu_to_keep.append(f"{i}_HETCOUNT")
        dejavu_to_keep.append(f"{i}_HOMCOUNT")
        dejavu_to_keep.append(f"{i}_ALLELEFREQ")
        dejavu_to_keep.append(f"{i}_SAMPLECOUNT")

    index_to_keep = []
    alphanumerical_list_index = []
    last_order_exec = False

    with open(output_file, "w") as write_file:
        with open(input_file, "r") as read_file:
            for line in read_file:
                line = line.rstrip("\n").split("\t")
                if line[0] == "chr":
                    for i in range(len(line)):
                        if last_order_exec == False and line[i] not in values_to_delete and not (line[i].endswith("_ALLELECOUNT") or line[i].endswith("_HETCOUNT") or line[i].endswith("_HOMCOUNT") or line[i].endswith("_ALLELEFREQ") or line[i].endswith("_SAMPLECOUNT")) :
                            if line[i] == last_order:
                                index_to_keep.append(i)
                                last_order_exec = True
                            else:
                                index_to_keep.append(i)
                        elif last_order_exec == False and line[i].endswith("_ALLELECOUNT") or line[i].endswith("_HETCOUNT") or line[i].endswith("_HOMCOUNT") or line[i].endswith("_ALLELEFREQ") or line[i].endswith("_SAMPLECOUNT"):
                            if line[i] in dejavu_to_keep:
                                index_to_keep.append(i)
                            else:
                                continue
                        elif last_order_exec == True:
                            if line[i].endswith("_ALLELECOUNT") or line[i].endswith("_HETCOUNT") or line[i].endswith("_HOMCOUNT") or line[i].endswith("_ALLELEFREQ") or line[i].endswith("_SAMPLECOUNT"):
                                if line[i] in dejavu_to_keep:
                                    alphanumerical_list_index.append(i)
                            else:
                                alphanumerical_list_index.append(i)

                    alphanumerical_list = [line[i] for i in alphanumerical_list_index]
                    alphanumerical_list.sort()
                    alphanumerical_list_new_index = []
                    for i in alphanumerical_list:
                        for j in range(len(line)):
                            if line[j] == i:
                                alphanumerical_list_new_index.append(j)
                            else:
                                continue
                    index_to_keep = index_to_keep + alphanumerical_list_new_index
                    write_file.write("\t".join([line[i] for i in index_to_keep]) + "\n")
                else:
                    # print(line)
                    for count, element in enumerate(line):
                        if "/" in element:
                            line[count] = f'="{element}"'
                        if "." in element:
                            element = element.split(".")
                            if element[0].isdigit() and element[1].isdigit():
                                line[count] = element[0] + "," + element[1]
                    # print("\t".join([line[i] for i in index_to_keep]) + "\n")
                    write_file.write("\t".join([line[i] for i in index_to_keep]) + "\n")
    os.remove(input_file)
    os.rename(output_file, input_file)

def cleaner(run_informations):
    log.info("Moving results from temporary folder")
    results_files = glob.glob(
        osj(run_informations["tmp_analysis_folder"], "*tsv")
    ) + glob.glob(osj(run_informations["tmp_analysis_folder"], "*vcf.gz"))

    if os.path.isdir(run_informations["archives_results_folder"]):
        log.info("Removing old results folder in archives")
        shutil.rmtree(run_informations["archives_results_folder"])
        os.mkdir(run_informations["archives_results_folder"])
    else:
        os.mkdir(run_informations["archives_results_folder"])

    for results_file in results_files:
        os.chmod(results_file, 0o777)
        log.info(
            f"Moving {results_file} to {run_informations["archives_results_folder"]}"
        )
        shutil.move(results_file, run_informations["archives_results_folder"])

    with open(
        osj(run_informations["archives_results_folder"], "VANNOTCopyComplete.txt"),
        mode="a",
    ):
        pass

    shutil.rmtree(run_informations["tmp_analysis_folder"])
    log.info("Deleted temporary analysis folder")


def panel_filtering(run_informations):
    tmp_vcf_files = glob.glob(
        osj(run_informations["tmp_analysis_folder"], "*VANNOT_*vcf.gz")
    )
    panels = run_informations["run_panels"]

    for panel in panels:
        subprocess.call(
            ["rsync", "-rvt", panel, run_informations["tmp_analysis_folder"]],
            universal_newlines=True,
        )
        panel = os.path.join(
            run_informations["tmp_analysis_folder"], os.path.basename(panel)
        )
        if os.path.basename(panel).split(".")[3] != "genes":
            panel_name = os.path.basename(panel).split(".")[1]
        else:
            panel_name = "_".join(os.path.basename(panel).split(".")[1].split("_")[1:])
        for tmp_vcf_file in tmp_vcf_files:
            sample_name = tmp_vcf_file.split(".")[0]
            filtered_vcf = osj(
                run_informations["tmp_analysis_folder"],
                sample_name + ".panel." + panel_name + ".vcf",
            )
            command_list = ["intersectBed", "-a", tmp_vcf_file, "-b", panel, "-header"]
            log.info(" ".join(command_list))
            with open(filtered_vcf, "a") as f:
                subprocess.call(
                    command_list,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    universal_newlines=True,
                )

            subprocess.call(["bgzip", filtered_vcf])
            filtered_vcf = filtered_vcf + ".gz"
                    

if __name__ == "__main__":
    pass
