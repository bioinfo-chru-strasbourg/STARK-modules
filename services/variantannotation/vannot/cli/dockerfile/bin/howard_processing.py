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
import synchronizer

def ignore_samples(run_informations):
    module_config = osj(
        os.environ["HOST_MODULE_CONFIG"],
        f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
    )
    ignored_samples = []
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        ignored_samples = data["ignored_samples"]
    if run_informations["type"] == "run":
        samplesheet = synchronizer.find_samplesheet(run_informations)
        control_samples = synchronizer.find_tag(samplesheet, "CQI#")
        ignored_samples = ignored_samples + control_samples
        
    log.info(
        "Ignoring following sample patterns for the merge vcf file : "
        + ", ".join(ignored_samples)
    )
    return ignored_samples

def format_to_qual_filter_id(run_informations):
    """
    Reverse of qual_filter_id_to_format.
    Take the VID, VQUAL and VFILTER values stored in the FORMAT column of each
    per-sample VCF and write them back into the site-level ID, QUAL and FILTER
    columns, then drop the VID/VQUAL/VFILTER FORMAT subfields and their headers.
    Replaces the backup_sample/restore_sample pair: QUAL/FILTER/ID travel with
    the variants through merge/annotate/unmerge instead of a side backup folder.
    """
    restore_fields = ("VID", "VQUAL", "VFILTER")
    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))
    for vcf_file in vcf_files:
        if "merged" in os.path.basename(vcf_file):
            continue
        log.info(
            f"Restoring QUAL, FILTER and ID from the FORMAT column for {os.path.basename(vcf_file)}"
        )

        is_gz = vcf_file.endswith(".gz")
        open_func = gzip.open if is_gz else open
        base = os.path.basename(vcf_file)[:-3] if is_gz else os.path.basename(vcf_file)
        tmp_output = osj(os.path.dirname(vcf_file), "qfirestore_" + base)

        with open_func(vcf_file, "rt") as read_file, open(tmp_output, "w") as write_file:
            for line in read_file:
                if line.startswith("##FORMAT=<ID="):
                    match = re.match(r"##FORMAT=<ID=([^,]+),", line)
                    if match and match.group(1) in restore_fields:
                        continue  # drop VID/VQUAL/VFILTER header definitions
                    write_file.write(line)
                elif line.startswith("#"):
                    write_file.write(line)
                else:
                    parts = line.rstrip("\n").split("\t")
                    fmt_keys = parts[8].split(":")
                    fmt_map = dict(zip(fmt_keys, parts[9].split(":")))

                    # restore the original values verbatim
                    if "VID" in fmt_map:
                        parts[2] = fmt_map["VID"]
                    if "VQUAL" in fmt_map:
                        parts[5] = fmt_map["VQUAL"]
                    if "VFILTER" in fmt_map:
                        parts[6] = fmt_map["VFILTER"]

                    # remove VID/VQUAL/VFILTER from FORMAT and every sample column
                    keep_idx = [i for i, k in enumerate(fmt_keys) if k not in restore_fields]
                    parts[8] = ":".join(fmt_keys[i] for i in keep_idx)
                    for col in range(9, len(parts)):
                        values = parts[col].split(":")
                        parts[col] = ":".join(values[i] for i in keep_idx if i < len(values))

                    write_file.write("\t".join(parts) + "\n")

        os.remove(vcf_file)
        if is_gz:
            subprocess.call(["bgzip", tmp_output], universal_newlines=True)
            os.rename(tmp_output + ".gz", vcf_file)
        else:
            os.rename(tmp_output, vcf_file)

def qual_filter_id_to_format(run_informations):
    """
    Copy the site-level QUAL, FILTER and ID values into the FORMAT column of
    each per-sample VCF as new FORMAT subfields (VID, VQUAL, VFILTER).
    HOWARD cannot move QUAL/FILTER/ID into FORMAT, so it is done with bcftools.
    Must be called while the per-sample VCFs are present in the tmp folder
    (e.g. right after restore_sample).
    """
    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))
    for vcf_file in vcf_files:
        if "merged" in os.path.basename(vcf_file):
            continue
        log.info(
            f"Copying QUAL, FILTER and ID into the FORMAT column for {os.path.basename(vcf_file)}"
        )
        sample = (
            subprocess.run(
                ["bcftools", "query", "-l", vcf_file],
                universal_newlines=True,
                stdout=subprocess.PIPE,
            )
            .stdout.strip()
            .split("\n")[0]
        )

        tmp_annot = osj(os.path.dirname(vcf_file), "qfi_annot.txt.tmp")
        tmp_annot_fixed = osj(os.path.dirname(vcf_file), "qfi_annot.fixed.txt.tmp")
        tmp_hdr = osj(os.path.dirname(vcf_file), "qfi_hdr.txt.tmp")

        # Extract site-level ID, QUAL and FILTER per variant
        query = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%ID\\t%QUAL\\t%FILTER\n"
        cmd = ["bcftools", "query", "-f", query, vcf_file]
        with open(tmp_annot, "w") as writefile:
            subprocess.call(cmd, universal_newlines=True, stdout=writefile)

        # ";" is not the FORMAT separator (":"), so ID/QUAL/FILTER can be
        # stored verbatim inside the FORMAT subfield without loss
        with open(tmp_annot, "r") as readfile, open(tmp_annot_fixed, "w") as writefile:
            for line in readfile:
                writefile.write(line)
        subprocess.call(["bgzip", tmp_annot_fixed], universal_newlines=True)
        tmp_annot_fixed = tmp_annot_fixed + ".gz"
        subprocess.call(
            ["tabix", "-s1", "-b2", "-e2", tmp_annot_fixed], universal_newlines=True
        )

        # Header defining the new FORMAT fields
        with open(tmp_hdr, "w") as writefile:
            writefile.write('##FORMAT=<ID=VID,Number=1,Type=String,Description="Original variant ID">\n')
            writefile.write('##FORMAT=<ID=VQUAL,Number=1,Type=String,Description="Original variant QUAL">\n')
            writefile.write('##FORMAT=<ID=VFILTER,Number=1,Type=String,Description="Original variant FILTER">\n')

        output_file = osj(
            os.path.dirname(vcf_file), "qfi_" + os.path.basename(vcf_file)
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
            "CHROM,POS,REF,ALT,FORMAT/VID,FORMAT/VQUAL,FORMAT/VFILTER",
            "-O",
            "z",
            "-o",
            output_file,
            vcf_file,
        ]
        log.debug(" ".join(cmd))
        subprocess.call(cmd, universal_newlines=True)

        os.remove(tmp_annot)
        os.remove(tmp_annot_fixed)
        os.remove(tmp_annot_fixed + ".tbi")
        os.remove(tmp_hdr)
        os.remove(vcf_file)
        os.rename(output_file, vcf_file)

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

        if run_informations["run_platform_application"] != None and len(sample_list) >= 1 and run_informations["onco"] == False:
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
            field_id_match = re.search(r"ID=([^,]+)", j)
            if field_id_match:
                field_id = field_id_match.group(1)
                if re.fullmatch(i, field_id):
                    info_to_keep.append("INFO/" + field_id)

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

    if run_informations["onco"] == True:
        log.info("Keeping sample based annotations from STARK into format column if any remaining")
        #if there are duplicates, keep both but renaming the incoming one adding a digit to the end of the name (BCF -> BCF2, BCFS -> BCFS2)
        #howard calculation --input=2506278.final.vcf.gz --output=test.vcf --calculations="INFO_TO_FORMAT" --param='{"calculation": {"calculations": {"INFO_TO_FORMAT": {"annotation_fields": {"ADP":null},"remove_info_fields": true}}}}'

        with open(module_config, "r") as read_file:
            data = json.load(read_file)
            annotation_fields = data["howard_ift"][
                run_informations["run_platform_application"]
            ]
        cleaned_vcf = convert_flag_to_integer(cleaned_vcf, annotation_fields)
        json_query = json.dumps(
            {
                "calculation": {
                    "calculations": {
                        "INFO_TO_FORMAT": {
                            "annotation_fields": annotation_fields,
                            "remove_info_fields": True,
                        }
                    }
                }
            }
        )

        howard_config = osj(
            os.environ["HOST_MODULE_CONFIG"], "howard", "howard_onco_config.json"
        )

        threads = commons.get_threads("threads_annotation")
        memory = commons.get_memory("memory_annotation")

        exact_time = time.time() + 7200
        local_time = time.localtime(exact_time)
        actual_time = time.strftime("%H%M%S", local_time)
        start = actual_time
        vcf_file = vcf_file.replace(".vcf.gz", ".vcf")
        container_name = f"VANNOT_itf_{start}_{run_informations['run_name']}_{os.path.basename(vcf_file).split('.')[0]}"
        launch_annotate_arguments = [
            "calculation",
            "--input",
            cleaned_vcf,
            "--output",
            vcf_file,
            "--calculations",
            "INFO_TO_FORMAT",
            "--param",
            json_query,
            "--memory",
            memory,
            "--threads",
            threads,
            "--config",
            howard_config,
            "--debug",
        ]

        log.info("ITF generation with HOWARD")
        howard_launcher.launch(container_name, launch_annotate_arguments)
        os.rename(vcf_file, cleaned_vcf)

    subprocess.call(["bgzip", cleaned_vcf], universal_newlines=True)
    renamed_clean = osj(
        os.path.dirname(cleaned_vcf),
        os.path.basename(cleaned_vcf).replace("cleaned_", "") + ".gz",
    )
    cleaned_vcf = cleaned_vcf + ".gz"
    os.rename(cleaned_vcf, renamed_clean)
    return renamed_clean

def convert_integer_to_flag(vcf_file, flag_fields):
    """
    Reverse of convert_flag_to_integer.
    FIELD=1 -> FLAG present (FIELD, no value)
    FIELD=0 -> FLAG absent (removed from INFO)
    Header: Number=1,Type=Integer -> Number=0,Type=Flag
    """
    flag_fields = set(flag_fields)
    if not flag_fields:
        return vcf_file

    is_gz = vcf_file.endswith(".gz")
    open_func = gzip.open if is_gz else open

    log.info(f"Converting Integer INFO fields back to FLAG: {', '.join(sorted(flag_fields))}")

    if is_gz:
        tmp_output = osj(os.path.dirname(vcf_file), "flagrestore_" + os.path.basename(vcf_file)[:-3])
    else:
        tmp_output = osj(os.path.dirname(vcf_file), "flagrestore_" + os.path.basename(vcf_file))

    info_pattern = re.compile(r"##INFO=<ID=([^,]+),")

    with open_func(vcf_file, "rt") as read_file, open(tmp_output, "w") as write_file:
        for line in read_file:
            if line.startswith("##INFO=<ID="):
                match = info_pattern.match(line)
                if match and match.group(1) in flag_fields:
                    line = re.sub(r"Number=[^,]+", "Number=0", line, count=1)
                    line = re.sub(r"Type=[^,]+", "Type=Flag", line, count=1)
                write_file.write(line)
            elif line.startswith("#"):
                write_file.write(line)
            else:
                parts = line.rstrip("\n").split("\t")
                if len(parts) > 7:
                    info_items = parts[7].split(";") if parts[7] != "." else []
                    new_info = []
                    for item in info_items:
                        key = item.split("=")[0]
                        if key in flag_fields:
                            value = item.split("=")[1] if "=" in item else ""
                            if value == "1":
                                new_info.append(key)  # flag present
                            # "0" / autre -> drop (flag absent)
                        else:
                            new_info.append(item)
                    parts[7] = ";".join(new_info) if new_info else "."
                    write_file.write("\t".join(parts) + "\n")
                else:
                    write_file.write(line)

    os.remove(vcf_file)
    if is_gz:
        subprocess.call(["bgzip", tmp_output], universal_newlines=True)
        os.rename(tmp_output + ".gz", vcf_file)
    else:
        os.rename(tmp_output, vcf_file)

    return vcf_file


def restore_flags_samples(run_informations):
    """
    Restore Integer -> Flag on every per-sample VCF in the tmp folder,
    using the converted_flags.json record produced by convert_flag_to_integer.
    Must be called while the per-sample VCFs still exist (before merge step 2).
    """
    record_file = osj(run_informations["tmp_analysis_folder"], "converted_flags.json")
    if not os.path.isfile(record_file):
        log.info("No converted_flags.json record found, skipping flag restoration on samples")
        return
    with open(record_file, "r") as rf:
        flag_fields = json.load(rf)
    if not flag_fields:
        log.info("Empty converted flags record, nothing to restore on samples")
        return

    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))
    for vcf_file in vcf_files:
        log.info(f"Restoring flags on sample VCF {os.path.basename(vcf_file)}")
        convert_integer_to_flag(vcf_file, flag_fields)
    os.remove(record_file)

def convert_flag_to_integer(vcf_file, annotation_fields):
    """
    Convert FLAG type INFO fields to Integer before HOWARD INFO_TO_FORMAT.
    FLAG present  -> FIELD=1
    FLAG absent   -> FIELD=0
    Header: Number=0,Type=Flag -> Number=1,Type=Integer
    annotation_fields: dict (keys = field names) or set/list
    """
    is_gz = vcf_file.endswith(".gz")
    open_func = gzip.open if is_gz else open

    field_names = set(annotation_fields.keys()) if isinstance(annotation_fields, dict) else set(annotation_fields)

    flag_pattern = re.compile(r"##INFO=<ID=([^,]+),[^>]*Type=Flag")
    flag_fields = set()

    with open_func(vcf_file, "rt") as read_file:
        for line in read_file:
            if not line.startswith("#"):
                break
            match = flag_pattern.match(line)
            if match and match.group(1) in field_names:
                flag_fields.add(match.group(1))

    if not flag_fields:
        log.info("No FLAG fields to convert")
        return vcf_file

    log.info(f"Converting FLAG INFO fields to Integer: {', '.join(sorted(flag_fields))}")

    # Persist the converted flag fields so they can be restored later
    record_file = osj(os.path.dirname(vcf_file), "converted_flags.json")
    existing = set()
    if os.path.isfile(record_file):
        with open(record_file, "r") as rf:
            existing = set(json.load(rf))
    existing |= flag_fields
    with open(record_file, "w") as wf:
        json.dump(sorted(existing), wf)

    if is_gz:
        tmp_output = osj(os.path.dirname(vcf_file), "flagfix_" + os.path.basename(vcf_file)[:-3])
    else:
        tmp_output = osj(os.path.dirname(vcf_file), "flagfix_" + os.path.basename(vcf_file))

    with open_func(vcf_file, "rt") as read_file, open(tmp_output, "w") as write_file:
        for line in read_file:
            if line.startswith("##INFO=<ID="):
                match = flag_pattern.match(line)
                if match and match.group(1) in flag_fields:
                    line = re.sub(r"Number=0", "Number=1", line, count=1)
                    line = re.sub(r"Type=Flag", "Type=Integer", line, count=1)
                write_file.write(line)
            elif line.startswith("#"):
                write_file.write(line)
            else:
                parts = line.rstrip("\n").split("\t")
                if len(parts) > 7:
                    info_items = parts[7].split(";") if parts[7] != "." else []
                    new_info = []
                    flags_found = set()

                    for item in info_items:
                        field_name = item.split("=")[0]
                        if field_name in flag_fields:
                            new_info.append(f"{field_name}=1")
                            flags_found.add(field_name)
                        else:
                            new_info.append(item)

                    # Absent flags -> 0
                    for flag in flag_fields:
                        if flag not in flags_found:
                            new_info.append(f"{flag}=0")

                    parts[7] = ";".join(new_info) if new_info else "."
                    write_file.write("\t".join(parts) + "\n")
                else:
                    log.warning("No FORMAT field")
                    write_file.write(line)

    os.remove(vcf_file)
    if is_gz:
        subprocess.call(["bgzip", tmp_output], universal_newlines=True)
        os.rename(tmp_output + ".gz", vcf_file)
    else:
        os.rename(tmp_output, vcf_file)

    return vcf_file

def restore_merged_header(run_informations, merged_vcf, reference_vcf=None):
    """
    Put the correct FORMAT/INFO definitions back onto the merged VCF using a
    corrected per-sample VCF as reference. The merged VCF body (and the FORMAT
    values of every sample) is left untouched: no FORMAT_TO_INFO is applied here.
    Must be called while the per-sample VCFs still exist.
    """
    if reference_vcf is None:
        candidates = glob.glob(
            osj(run_informations["tmp_analysis_folder"], "VANNOT_*.vcf.gz")
        )
        candidates = [
            v for v in candidates
            if os.path.realpath(v) != os.path.realpath(merged_vcf)
        ]
        if not candidates:
            log.warning("No corrected per-sample VCF found, merged header left unchanged")
            return merged_vcf
        reference_vcf = candidates[0]

    log.info(f"Restoring merged header from {os.path.basename(reference_vcf)}")

    # Collect correct FORMAT/INFO definitions from the reference sample VCF
    reference_header = subprocess.run(
        ["bcftools", "view", "-h", reference_vcf],
        universal_newlines=True, stdout=subprocess.PIPE,
    ).stdout
    correct_def = {}  # (FORMAT|INFO, ID) -> full meta line
    for line in reference_header.splitlines(keepends=True):
        match = re.match(r"##(FORMAT|INFO)=<ID=([^,]+),", line)
        if match:
            correct_def[(match.group(1), match.group(2))] = line

    # Rebuild the merged header, swapping only matching FORMAT/INFO lines and
    # keeping the merged #CHROM line (all sample columns) intact.
    merged_header = subprocess.run(
        ["bcftools", "view", "-h", merged_vcf],
        universal_newlines=True, stdout=subprocess.PIPE,
    ).stdout
    new_hdr = osj(os.path.dirname(merged_vcf), "merged_hdr.txt.tmp")
    with open(new_hdr, "w") as writefile:
        for line in merged_header.splitlines(keepends=True):
            match = re.match(r"##(FORMAT|INFO)=<ID=([^,]+),", line)
            key = (match.group(1), match.group(2)) if match else None
            writefile.write(correct_def.get(key, line) if key else line)

    restored_merged = osj(
        os.path.dirname(merged_vcf), "reheader_" + os.path.basename(merged_vcf)
    )
    subprocess.call(
        ["bcftools", "reheader", "-h", new_hdr, "-o", restored_merged, merged_vcf],
        universal_newlines=True,
    )
    os.remove(new_hdr)
    os.remove(merged_vcf)
    os.rename(restored_merged, merged_vcf)
    log.info(f"Correct header restored on merged VCF {os.path.basename(merged_vcf)}")
    return merged_vcf

def format_to_info(run_informations):
    module_config = run_informations["module_config"]
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        annotation_fields = data["format_to_info"]
        if run_informations["run_platform_application"] not in annotation_fields:
            log.info(
                f"No FORMAT_TO_INFO annotations defined for {run_informations['run_platform_application']}, skipping"
            )
            return

    format_to_info_columns_config = annotation_fields[
        run_informations["run_platform_application"]
    ]
    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))
    for vcf_file in vcf_files:
        if "merged" in os.path.basename(vcf_file):
            vcf_files.remove(vcf_file)

    for vcf_file in vcf_files:
        log.info(
            f"Copying FORMAT fields into the INFO column for {os.path.basename(vcf_file)}"
        )

        # Filter to FORMAT fields that actually exist in this VCF
        vcf_header = subprocess.run(
            ["bcftools", "view", "-h", vcf_file],
            universal_newlines=True,
            stdout=subprocess.PIPE,
        ).stdout
        existing_format_fields = set(re.findall(r"##FORMAT=<ID=([^,]+)", vcf_header))
        format_to_info_columns = [f for f in format_to_info_columns_config if f in existing_format_fields]

        if not format_to_info_columns:
            log.info(f"No matching FORMAT fields found in {os.path.basename(vcf_file)}, skipping")
            continue

        log.info(f"FORMAT fields to copy to INFO: {format_to_info_columns}")

        sample = (
            subprocess.run(
                ["bcftools", "query", "-l", vcf_file],
                universal_newlines=True,
                stdout=subprocess.PIPE,
            )
            .stdout.strip()
            .split("\n")[0]
        )
        tmp_annot = osj(os.path.dirname(vcf_file), "annot.txt.tmp")
        tmp_hdr = osj(os.path.dirname(vcf_file), "hdr.txt.tmp")

        # Extract the FORMAT values (per site) for the first sample
        format_to_info_columns_query = "\\t%".join(format_to_info_columns)
        format_to_info_columns_query = (
            "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%" + format_to_info_columns_query + "]\n"
        )
        cmd = [
            "bcftools",
            "query",
            "-s",
            sample,
            "-f",
            format_to_info_columns_query,
            vcf_file,
        ]
        with open(tmp_annot, "w") as writefile:
            subprocess.call(cmd, universal_newlines=True, stdout=writefile)
        subprocess.call(["bgzip", tmp_annot], universal_newlines=True)
        tmp_annot = tmp_annot + ".gz"
        subprocess.call(
            ["tabix", "-s1", "-b2", "-e2", tmp_annot], universal_newlines=True
        )

        # Build matching INFO header lines from the existing FORMAT definitions
        original_info_headers = {}  # column -> original INFO header line
        vcf_file_gunzip = vcf_file[:-3]
        subprocess.call(["gunzip", vcf_file])
        with open(vcf_file_gunzip, "r") as readfile, open(tmp_hdr, "w") as writefile:
            for line in readfile:
                if not line.startswith("#"):
                    break
                for column in format_to_info_columns:
                    if line.startswith("##FORMAT=<ID=" + column + ","):
                        # Save the correct INFO line for later restoration
                        original_info_headers[column] = line.replace("##FORMAT=<ID=", "##INFO=<ID=", 1)
                        # Write with forced Number=.,Type=String to avoid bcftools segfault
                        new_line = line.replace("##FORMAT=<ID=", "##INFO=<ID=", 1)
                        new_line = re.sub(r"Number=[^,]+", "Number=.", new_line, count=1)
                        new_line = re.sub(r"Type=[^,]+", "Type=String", new_line, count=1)
                        writefile.write(new_line)
        subprocess.call(["bgzip", vcf_file_gunzip], universal_newlines=True)

        format_to_info_columns_annotate = ",INFO/".join(format_to_info_columns)
        format_to_info_columns_annotate = (
            "CHROM,POS,REF,ALT,INFO/" + format_to_info_columns_annotate
        )

        # Add the new INFO fields without touching the FORMAT column
        output_file = osj(
            os.path.dirname(vcf_file), "fti_" + os.path.basename(vcf_file)
        )
        cmd = [
            "bcftools",
            "annotate",
            "-a",
            tmp_annot,
            "-h",
            tmp_hdr,
            "-c",
            format_to_info_columns_annotate,
            "-O",
            "z",
            "-o",
            output_file,
            vcf_file,
        ]
        log.debug(" ".join(cmd))
        subprocess.call(cmd, universal_newlines=True)

        # Restore the correct Number/Type in the output header
        if original_info_headers:
            restore_hdr = osj(os.path.dirname(vcf_file), "restore_hdr.txt.tmp")
            current_header = subprocess.run(
                ["bcftools", "view", "-h", output_file],
                universal_newlines=True,
                stdout=subprocess.PIPE,
            ).stdout
            with open(restore_hdr, "w") as writefile:
                for hdr_line in current_header.splitlines(keepends=True):
                    replaced = False
                    for column, original_line in original_info_headers.items():
                        if hdr_line.startswith("##INFO=<ID=" + column + ","):
                            writefile.write(original_line)
                            replaced = True
                            break
                    if not replaced:
                        writefile.write(hdr_line)

            restored_output = osj(os.path.dirname(output_file), "restored_" + os.path.basename(output_file))
            subprocess.call([
                "bcftools", "reheader",
                "-h", restore_hdr,
                "-o", restored_output,
                output_file,
            ], universal_newlines=True)
            os.remove(output_file)
            os.remove(restore_hdr)
            os.rename(restored_output, output_file)
            log.info("Restored original Number/Type in INFO header after FORMAT_TO_INFO")

        os.remove(tmp_annot)
        os.remove(tmp_annot + ".tbi")
        os.remove(tmp_hdr)
        os.remove(vcf_file)
        os.rename(output_file, vcf_file)
    
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

def normalize_merge_headers(vcf_file, header_backup):
    """
    Rewrite VAF/AD/FT/PL/AS_FilterStatus/reference header lines to a single
    canonical definition so every VCF agrees before bcftools merge.
    header_backup: a plain dict owned by the caller, mutated in place. The
    first original line seen for each field is kept; later calls on an
    already-normalized file won't overwrite it.
    """
    field_descriptions = {
        "VAF": "VAF Variant Frequency, calculated from quality",
        "AD": "Allelic depths for the ref and alt alleles in the order listed",
        "FT": "Genotype-level filter",
        "PL": (
            "Normalized, Phred-scaled likelihoods for genotypes "
            "as defined in the VCF specification"
        ),
        "AS_FilterStatus": (
            "Filter status for each allele, as assessed by ApplyVQSR. "
            "Note that the VCF filter field will reflect the most lenient/sensitive "
            "status across all alleles."
        ),
    }
    merge_format_types = {
        "VAF": ("1", "Float"),
        "AD": ("R", "Integer"),
        "FT": (".", "String"),
        "PL": ("G", "Integer"),
        "AS_FilterStatus": (".", "String"),
    }
    reference_prefix = "##reference=file:"
    reference_line = "##reference=file:///STARK/databases/genomes/current/hg19.fa\n"

    def format_prefix(field):
        return f"##FORMAT=<ID={field},Number="

    def format_line(field):
        number, type_ = merge_format_types[field]
        return (
            f'##FORMAT=<ID={field},Number={number},Type={type_},'
            f'Description="{field_descriptions[field]}">\n'
        )

    tmp_output = osj(
        os.path.dirname(vcf_file), "hdrfix_" + os.path.basename(vcf_file)[:-3]
    )
    with gzip.open(vcf_file, "rt") as read_file, open(tmp_output, "w") as write_file:
        for line in read_file:
            if line.startswith(reference_prefix):
                if reference_prefix not in header_backup:
                    header_backup[reference_prefix] = line
                    log.info("Saving original header line for reference")
                write_file.write(reference_line)
                continue

            field = next(
                (f for f in merge_format_types if line.startswith(format_prefix(f))),
                None,
            )
            if field is None:
                write_file.write(line)
                continue

            prefix = format_prefix(field)
            if prefix not in header_backup:
                header_backup[prefix] = line
                log.info(f"Saving original header line for {field}")
            write_file.write(format_line(field))

    os.remove(vcf_file)
    subprocess.call(["bgzip", tmp_output], universal_newlines=True)
    os.rename(tmp_output + ".gz", vcf_file)


def restore_merge_headers(run_informations, header_backup):
    """
    Reverse of normalize_merge_headers: put the original FORMAT/reference
    header lines back on every vcf.gz in the tmp folder (per-sample AND
    merged). header_backup is the same dict filled in earlier by merge_vcf/
    normalize_merge_headers during this run. Must run AFTER restore_merged_header.
    """
    if not header_backup:
        log.info("No merge header backup collected, skipping header restoration")
        return

    vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))
    for vcf_file in vcf_files:
        log.info(f"Restoring original header lines for {os.path.basename(vcf_file)}")
        tmp_output = osj(
            os.path.dirname(vcf_file), "hdrrestore_" + os.path.basename(vcf_file)[:-3]
        )
        with gzip.open(vcf_file, "rt") as read_file, open(tmp_output, "w") as write_file:
            for line in read_file:
                prefix = next((p for p in header_backup if line.startswith(p)), None)
                write_file.write(header_backup[prefix] if prefix else line)

        os.remove(vcf_file)
        subprocess.call(["bgzip", tmp_output], universal_newlines=True)
        os.rename(tmp_output + ".gz", vcf_file)

def merge_vcf(run_informations, step, base_vcf, header_backup=None):
    if header_backup is None:
        header_backup = {}
    if step == "0":
        vcf_file_to_merge = glob.glob(osj(run_informations["tmp_analysis_folder"], "fixed_unmerged_*.vcf*"))
        if len(vcf_file_to_merge) == 0:
            vcf_file_to_merge = glob.glob(osj(run_informations["tmp_analysis_folder"], "unmerged_*.vcf*"))
    elif step == "2":
        ignored_samples = ignore_samples(run_informations)
        vcf_file_to_merge = glob.glob(
            osj(run_informations["tmp_analysis_folder"], "*.vcf*")
        )
        for ignored_sample in ignored_samples:
            for vcf_file in vcf_file_to_merge:
                if ignored_sample in os.path.basename(vcf_file):
                    vcf_file_to_merge.remove(vcf_file)
                    log.info(f"Sample {ignored_sample} is ignored, not included in the merge")
    else:
        print("sam: running merge_vcf step=1")
        print("run_informations['tmp_analysis_folder']", run_informations["tmp_analysis_folder"])
        vcf_file_to_merge = glob.glob(
            osj(run_informations["tmp_analysis_folder"], "*.vcf*")
        )
        
    print(vcf_file_to_merge)
    if len(vcf_file_to_merge) > 1:
        for vcf_file in vcf_file_to_merge:
            normalize_merge_headers(vcf_file, header_backup)

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
            cmd = ["bcftools", "view", "-c1", "-I", "-s", sample, vcf_file_to_unmerge]
            
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

    if run_informations["onco"] == True:
        howard_config = osj(
            os.environ["HOST_MODULE_CONFIG"], "howard", "howard_onco_config.json"
        )
    else:
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
        "--debug",
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
    if run_informations["onco"] == True:
        howard_config = osj(
            os.environ["HOST_MODULE_CONFIG"], "howard", "howard_onco_config.json"
        )
    else:
        howard_config = osj(os.environ["HOST_MODULE_CONFIG"], "howard", "howard_config.json")
        
    transcript_param = osj(os.environ["HOST_MODULE_CONFIG"],"howard",f"param.transcripts.{run_informations["run_platform"].lower()}.json",)

    if len(vcf_files) > 1:
        for vcf_file in vcf_files:
            os.rename(vcf_file, osj(os.path.dirname(vcf_file), os.path.basename(vcf_file).removeprefix("unmerged_")))
        vcf_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.vcf.gz"))

    threads = commons.get_threads("threads_annotation")
    memory = commons.get_memory("memory_annotation")

    exact_time = time.time() + 7200
    local_time = time.localtime(exact_time)
    actual_time = time.strftime("%H%M%S", local_time)
    start = actual_time

    with open(transcript_param, "r") as read_file:
        data = json.load(read_file)
    try:
        transcripts_output = data["calculation"]["TRANSCRIPTS_EXPORT"]["export"]["output"]
    except (KeyError, TypeError):
        log.warning(
            "transcripts export output not found in howard param; "
            "transcript TSV will not be copied"
        )
        transcripts_output = ""

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

        sample_name = (os.path.basename(vcf_file).split(".")[0]).removeprefix("VANNOT_")
        transcripts_output_renamed = f"VANNOT_transcripts_{sample_name}.tsv"
        if transcripts_output:
            shutil.copy(
                transcripts_output,
                osj(run_informations["tmp_analysis_folder"], transcripts_output_renamed),
            )
        print("howard_score_transcripts output", osj(run_informations["tmp_analysis_folder"], transcripts_output_renamed))


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
        fields_to_keep_raw = data.get("fields_to_keep_raw", {}).get(
            run_informations["run_platform_application"], []
        )

    explode_infos_fields = ",".join(re.escape(field) for field in ordered_fields)
    explode_infos_fields = explode_infos_fields + ",.*"
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
        # FORMAT + sample column(s) must go to the END of the TSV.
        # --explode_infos only adds INFO columns, so FORMAT/sample stay in their
        # original VCF position. DuckDB's "* EXCLUDE (...)" keeps every column but
        # lets us re-append FORMAT and the sample(s) last.
        sample = subprocess.run(
            ["bcftools", "query", "-l", vcf_file],
            universal_newlines=True,
            stdout=subprocess.PIPE,
        ).stdout.strip().split("\n")

        trailing_cols = ["FORMAT"] + sample
        trailing_quoted = ", ".join(f'"{c}"' for c in trailing_cols)
        select_query = f"SELECT * EXCLUDE ({trailing_quoted}), {trailing_quoted} FROM variants"
        
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
            force_info_fields_as_string(vcf_file, fields_to_keep_raw)
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
                select_query,
                "--threads",
                threads,
                "--memory",
                memory,
                "--config",
                howard_config,
            ]
            howard_launcher.launch(container_name, launch_convert_arguments)

            header_fixed = osj(
                os.path.dirname(output_file), "hdrfix_" + os.path.basename(output_file)
            )
            with open(output_file, "r") as read_file, open(header_fixed, "w") as write_file:
                is_header = True
                for line in read_file:
                    if is_header:
                        write_file.write(line.replace("\\", ""))
                        is_header = False
                    else:
                        write_file.write(line)
            os.replace(header_fixed, output_file)

            if run_informations["onco"] == False:
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
        # is_varscan_ad = False
        for values in format_dict.values():
            if "AD" in values[0]:
                ad_index = values[0].index("AD")
                ad_value = values[1][ad_index]
                if ad_value.count(",") == 2:
                    is_second_ad_alt = True
                # elif "," not in ad_value and ad_value != ".":
                #     is_varscan_ad = True

        modified_header = []
        new_tsv_lines = []

        for column_name in header:
            if column_name == "#CHROM":
                modified_header.append("chr")
            elif column_name == "AD":
                # if is_varscan_ad == False:
                column_name = "AD_ref"
                modified_header.append(column_name)
                modified_header.append("AD_alt")
                if is_second_ad_alt is True:
                    modified_header.append("AD_alt2")
                # else:
                #     modified_header.append(column_name)
            else:
                modified_header.append(column_name)

        ad_index = modified_header.index("AD_ref")

        for line in tsv_lines:
            modified_tsv_line = []
            line = line.rstrip("\n").split("\t")
            for count, content in enumerate(line):
                # if count == ad_index and is_varscan_ad is False:
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
                    if key == "AD":
                        ad_ref_index = index_dict["AD_ref"]
                        ad_alt_index = index_dict["AD_alt"]
                        if value.count(",") == 0:
                            ad_ref = value
                            line[ad_ref_index] = ad_ref
                            ad_alt = ""
                            line[ad_alt_index] = ad_alt
                        elif value.count(",") == 1:
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
    
def check_if_tsv_empty(input_file, run_informations):
    def get_sample_name(file_path):
        name = os.path.basename(file_path).removesuffix(".tsv").removeprefix("VANNOT_")
        if name.endswith(".design"):
            name = name.removesuffix(".design")
        elif len(name.split(".")) > 1 and name.split(".")[1] == "panel":
            name = name.split(".")[0]
        return name

    with open(input_file, "r") as read_file:
        first_line = read_file.readline()

    if first_line.strip() != "":
        return

    log.info(
        f"{os.path.basename(input_file)} is empty, retrieving a header from another tsv"
    )
    empty_sample = get_sample_name(input_file)
    tsv_files = glob.glob(osj(run_informations["tmp_analysis_folder"], "*.tsv"))

    for tsv_file in tsv_files:
        if tsv_file == input_file or "transcripts" in os.path.basename(tsv_file):
            continue
        with open(tsv_file, "r") as read_file:
            candidate_header = read_file.readline()
        if candidate_header.strip() == "":
            continue
        donor_sample = get_sample_name(tsv_file)
        header = [
            empty_sample if column == donor_sample else column
            for column in candidate_header.rstrip("\n").split("\t")
        ]
        with open(input_file, "w") as write_file:
            write_file.write("\t".join(header) + "\n")
        log.info(
            f"Used header from {os.path.basename(tsv_file)} for {os.path.basename(input_file)}"
        )
        return

    log.warning(
        f"No non-empty tsv found to retrieve a header for {os.path.basename(input_file)}"
    )

def force_info_fields_as_string(vcf_file, fields_to_keep_raw):
    """
    Rewrite the VCF header so the given INFO fields become Number=1,Type=String.
    Prevents HOWARD explode_infos from re-typing comma lists (Float -> adds '.0'
    to ints and drops trailing missing '.' entries). Keeps the raw INFO value
    exactly as written, like STARK.
    """
    if not fields_to_keep_raw:
        return vcf_file
    tmp_output = osj(
        os.path.dirname(vcf_file), "rawstr_" + os.path.basename(vcf_file)[:-3]
    )
    id_pattern = re.compile(r"##INFO=<ID=([^,]+),")
    with gzip.open(vcf_file, "rt") as read_file, open(tmp_output, "w") as write_file:
        for line in read_file:
            if line.startswith("##INFO=<ID="):
                match = id_pattern.match(line)
                if match and match.group(1) in fields_to_keep_raw:
                    line = re.sub(r"Number=[^,]+", "Number=.", line, count=1)
                    line = re.sub(r"Type=[^,]+", "Type=String", line, count=1)
            write_file.write(line)
    os.remove(vcf_file)
    subprocess.call(["bgzip", tmp_output], universal_newlines=True)
    os.rename(tmp_output + ".gz", vcf_file)
    return vcf_file

def tsv_modifier(input_file, run_informations):
    check_if_tsv_empty(input_file, run_informations)
    sample = os.path.basename(input_file).removesuffix(".tsv").removeprefix("VANNOT_")
    if sample.endswith(".design"):
        sample = sample.removesuffix(".design")
    elif len(sample.split(".")) > 1 and sample.split(".")[1] == "panel":
        sample = sample.split(".")[0]
    output_file = osj(
        os.path.dirname(input_file), "tmp_" + os.path.basename(input_file)
    )
    module_config = osj(os.environ["HOST_MODULE_CONFIG"],f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")

    if run_informations["onco"] == False:
        with open(module_config, "r") as read_file:
            data = json.load(read_file)
            values_to_delete = data["tsv_columns_to_remove_diag"]
    else:
        with open(module_config, "r") as read_file:
            data = json.load(read_file)
            values_to_delete = data["tsv_columns_to_remove_onco"]    

    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        ordered_fields = data["vcf_to_tsv_column_order"][
            run_informations["run_platform_application"]
        ]
        last_order = ordered_fields[-1] if ordered_fields else None

    dejavu_to_keep = []
    if run_informations["onco"] == False:
        values_to_delete.append(sample)
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
                if line[0] == "chr" or line[0] == "#CHROM":
                    for i in range(len(line)):
                        if last_order is None:
                            last_order_exec = True
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
                    if run_informations["run_platform"] != "DIAGGEN":
                        for count, element in enumerate(line): 
                            if "/" in element and element.count("/") == 1 and element.split("/")[0].isdigit() and element.split("/")[1].isdigit():
                                line[count] = f'="{element}"' #fix the problem with excel reading dates instead of fractions
                            elif "," not in element and element.count(".") == 1: #replace decimales with . to , not with lists
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
        if os.path.basename(panel).endswith(".manifest.genes"):
            panel_name = "_".join(os.path.basename(panel).split(".")[1].split("_")[1:])
        else:
            panel_name = os.path.basename(panel).split(".")[-2]
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
