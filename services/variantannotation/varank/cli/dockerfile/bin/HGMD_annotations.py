# -*- coding: utf-8 -*-

import argparse
import logging
import argparse
import pprint
import subprocess
import re
from os.path import join as osj
import os
import time
from pprint import pprint
import gzip
import sys


"""
Aim: Annotate raw VCF with HGMD pro annotations
"""


def split_vcf(vcf: str) -> tuple:
    """
    From vcf file, return list of header and number of mutations
    """
    header = []
    fields = []
    if vcf.endswith(".gz"):
        with gzip.open(vcf, "rb") as f:
            for j, l in enumerate(f):
                l = l.decode().strip()
                if l.startswith("##"):
                    header.append(l)
                elif l.startswith("#"):
                    fields.append(l)
    else:
        with open(vcf, "r") as f:
            for j, l in enumerate(f):
                if l.startswith("##"):
                    header.append(l)
                elif l.startswith("#"):
                    fields.append(l)
    mutations_count = j + 1 - (len(header) + 1)
    return (header, fields, mutations_count)


def explode_header(header: list, log: logging.Logger, notparse=[]) -> dict:
    """
    From header of vcf file return huge dico of header
    notparse: avoid splitting in misformatted field especially tier tools
    """
    dico = {}
    error = []
    # for each row in header
    for lines in header:
        lines = lines.strip()
        effective = lines.split("##")[1]
        # if row is fill and need to be translate in dict
        if effective.endswith(">") and all(
            [items for items in notparse if not effective.startswith(items)]
        ):
            desc = effective.split("=", 1)
            # if key id is not already in dict create key
            if not desc[0] in dico.keys():
                dico[desc[0]] = {}
            # for each value in field ID, maximum split 3 to avoid split in description field (description should be the last) TODO
            tmp = {}
            r = re.search("<(.*)>", desc[1]).group(1)
            try:
                wde = re.search("(.*),Description", r).group(1)
            except AttributeError:
                wde = False
            # if row contain Description (it should be)
            if wde:
                hook = wde.split(",")
                # print(hook)
                hook.append("Description=" + re.search("Description=(.*)", r).group(1))
            else:
                hook = r.split(",")
            # Becarefull if only one field
            for it in hook:
                field = it.split("=", 1)
                try:
                    key = field[0]
                    value = field[1]
                    tmp[key] = value
                except IndexError:
                    log.warning("probably wrong header row ", field)
                    # exit()
                # Check if Description field is correct
                if key == "Description":
                    if not value.startswith('"') or not value.endswith('"'):
                        log.error(value)
                        error.append(lines)
            # STACK ID value as dict name and other values key pair in level -1
            if "ID" in tmp.keys():
                value = tmp["ID"]
                tmp.pop("ID")
                dico[desc[0]][value] = tmp
            else:
                dico[desc[0]] = tmp
            # print(tmp)
        # extra field values
        elif not effective.startswith("contig"):
            val = effective.split("=", 1)
            dico[val[0]] = val[1]
    # print(json.dumps(dico, indent=4))
    return dico


def systemcall(command: str, log, skip=None) -> list:
    """
    Call bash command
    In case of crash exit code 1 stop script
    """
    log.info(command)
    p = subprocess.Popen(
        [command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    out, err = p.communicate()
    if not err:
        return out.decode("utf8").strip().split("\n")
    else:
        issues = err.decode("utf8").strip()
        try:
            re.search(r"(Warning|WARNING)", issues).group()
            log.warning(err.decode("utf-8").strip())
            return out.decode("utf8").strip().split("\n")
        except AttributeError:
            try:
                re.search(r"INFO", issues).group()
                log.info(err.decode("utf-8").strip())
                return out.decode("utf8").strip().split("\n")
            except AttributeError:
                log.error(err.decode("utf-8").strip())
                if skip:
                    pass
                else:
                    exit()


class PPrintFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None):
        super().__init__(fmt, datefmt)

    def format(self, record):
        if isinstance(record.msg, dict):
            formatted_msg = pprint.pformat(record.msg, sort_dicts=False, depth=2)
            record.msg = formatted_msg
        elif isinstance(record.msg, list):
            formatted_msg = pprint.pformat(
                [str(item) if isinstance(item, int) else item for item in record.msg],
                depth=1,
            )
            record.msg = formatted_msg
        return super().format(record)


class ShowHgmd(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        # print("list of fonts here")
        parser.exit()  # exits the program with no more arg parsing and checking


def show_hgmd():
    if sys.argv[1] == "--show":
        log = setup_logger(False, "INFO", "hgmd-annotations.log")
        header, fields, mutations_count = split_vcf(sys.argv[2])
        header_hgmd = explode_header(header, log, "GATK")
        pprint(header_hgmd, sort_dicts=False)
        sys.exit(1)


def setup_logger(
    log_to_file: bool, log_level: str, logfile_path=None
) -> logging.Logger:
    """
    Logging lib handling
    """
    logger = logging.getLogger("my_logger")
    logger.setLevel(log_level)

    formatter = PPrintFormatter("%(asctime)s - %(levelname)s - %(message)s")

    if log_to_file:
        file_handler = logging.FileHandler(logfile_path)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    else:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    return logger


def parseargs() -> argparse.Namespace:
    """
    argparse Lib
    """
    parser = argparse.ArgumentParser(description="HGMD annotator")
    parser.add_argument(
        "--logfile",
        action="store_true",
        help="Log messages to a file in the same output folder of vcf",
    )
    parser.add_argument(
        "--verbosity",
        "-v",
        choices=["debug", "info", "warning", "error", "critical"],
        default="info",
        help="Set verbosity level",
    )
    parser.add_argument("-i", "--input", type=str, help="Input raw vcf", required=True)
    parser.add_argument(
        "-o", "--output", type=str, help="Annotated vcf with HGMD", required=True
    )
    parser.add_argument(
        "--hgmd", type=str, help="Path of HGMD vcf database", required=True
    )
    parser.add_argument(
        "-a", "--assembly", type=str, choices=["hg19", "hg38"], help="Assembly version"
    )
    parser.add_argument(
        "--fields",
        type=str,
        help="List of fields to annotate, example --fields DM,PHEN default all fields",
        default="all",
    )
    parser.add_argument(
        "--bgzip", type=str, default="bgzip", help="Path of bgzip binary"
    )
    parser.add_argument(
        "--tabix", type=str, default="tabix", help="Path of tabix binary"
    )
    parser.add_argument(
        "--bcftools", type=str, default="bcftools", help="Path of bcftools binary"
    )
    parser.add_argument(
        "--show",
        nargs=1,
        action=show_hgmd(),
        help="--show <HGMD-vcf> to pprint hgmd header",
    )
    return parser.parse_args()


def bcftools_annotate(
    vcfin: str,
    vcfout: str,
    header_field: str,
    log,
    hgmd: str,
    bcftools,
    bgzip,
    tabix,
) -> str:
    """
    Arguments:

    variants: an object storing variant data to be annotated.
    vcfin: a string representing the input VCF file.
    vcfout: a string representing the output VCF file.
    header_field: a string indicating which field in the VCF file to annotate.
    vcftoannot=None: an optional string representing an additional VCF file to be used for carrying new annotations.

    Returns:
    vcfout: a string representing the path to the output VCF file.
    """
    if os.path.isfile(vcfout):
        os.remove(vcfout)
    if os.path.isfile(os.path.join(vcfout, ".tbi")):
        os.remove(os.path.join(vcfout, ".tbi"))

    log.info("bcftools annotate to " + vcfout)
    if not vcfin.endswith(".gz"):
        log.info("Compress and index input")
        systemcall(bgzip + " --force " + vcfin, log)
        systemcall(tabix + " -p vcf " + vcfin + ".gz", log)
        vcfin = vcfin + ".gz"
    else:
        if not os.path.isfile(os.path.join(vcfin + ".tbi")):
            systemcall(tabix + " -p vcf " + vcfin, log)     

    if not vcfout.endswith(".gz"):
        vcfout = vcfout + ".gz"
    systemcall(
        bcftools
        + " annotate -O 'z' -o "
        + vcfout
        + " -a "
        + hgmd
        + " -c CHROM,POS,REF,ALT,"
        + ",".join(header_field)
        + " "
        + vcfin,
        log,
    )
    systemcall(tabix + " -p vcf " + vcfout, log)
    return vcfout


def main():
    args = parseargs()
    output_folder = os.path.dirname(args.output)
    log = setup_logger(
        args.logfile, args.verbosity.upper(), osj(output_folder, "hgmd-annotations.log")
    )
    log.info("HGMD annotations")

    header, fields, mutations_count = split_vcf(args.hgmd)
    header_hgmd = explode_header(header, log, "GATK")

    log.info(header_hgmd["reference"])
    log.info(header_hgmd["source"])
    log.info("Variants HGMD " + str(mutations_count))
    log.info("\t")

    header_hgmd_info = list(header_hgmd["INFO"].keys())

    log.info("Available fields " + " ".join(header_hgmd_info))

    if args.fields == "all":
        chosen = header_hgmd_info
    else:
        chosen = args.fields.split(",")

    header_field = []
    header_field_missing = []
    for values in chosen:
        if values not in header_hgmd_info:
            header_field_missing.append(values)
        else:
            header_field.append(values)
    if not header_field:
        log.error("Fields values empty EXIT")
        sys.exit(1)
    log.info("Annotations " + ",".join(header_field))

    if header_field_missing:
        log.warning("Values missing " + ",".join(header_field_missing))
    bcftools_annotate(
        args.input,
        args.output,
        header_field,
        log,
        args.hgmd,
        args.bcftools,
        args.bgzip,
        args.tabix,
    )


if __name__ == "__main__":
    main()
