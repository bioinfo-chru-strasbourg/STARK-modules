# !/usr/bin/python
# !/bin/sh
# !/usr/tcl/bin
# -*- coding: utf-8 -*-
"""
##################################
##							    ##
##  VARIANTANNOTATION ANALYSIS  ##
##    Author : Mateusz RAUCH    ##
##							    ##
##################################
"""

import argparse
from locale import strcoll
import sys
import os
import logging as log
import commons
import checker
import howard_run
import howard_folder


def main(args):
    original_umask = os.umask(0o000)

    # print(mylog.handlers[0].baseFilename)

    # if "varank" is True:
    #     if "run" in args:
    #         print("run")
    #     if "folder" in args:
    #         print("folder")

    if "run" in args:
        commons.set_log_level_run(args)
        commons.set_logger_info_run(args)

        howard_run.launch_run(args)

        commons.set_logger_info_run(args)

    elif "folder" in args:
        commons.set_log_level_default(args)
        howard_folder.launch_folder(args)

    # elif args.allsync18:
    # 	pattern = pattern_generator(args)
    # 	launch_universal_sync18(pattern)
    os.umask(original_umask)


def parse_args():

    # Main parser
    main_parser = argparse.ArgumentParser(prog="variantannotation")

    # Secondary parsers with global arguments
    verbosity_parser = argparse.ArgumentParser(add_help=False)
    verbosity_parser.add_argument(
        "-v",
        "--verbosity",
        type=str,
        default="info",
        help="Verbosity level, DEBUG, INFO, WARNING, ERROR, CRITICAL. Default value is INFO",
    )
    pattern_parser = argparse.ArgumentParser(add_help=False)
    pattern_parser.add_argument(
        "-p",
        "--pattern",
        type=str,
        default=commons.default_pattern,
        help="pattern describing which vcf files to synchronize in STARK folders, actual patterns : */STARK/*.reports/*.final.vcf.gz (default), */POOL/*.final.vcf.gz. You can use two patterns with space as separator",
    )
    mode_parser = argparse.ArgumentParser(add_help=False)
    mode_parser.add_argument(
        "-m",
        "--launchmode",
        type=str,
        default="manual",
        help="hidden argument to know if the varank analysis is manual or listener launched",
    )
    varank_parser = argparse.ArgumentParser(add_help=False)
    varank_parser.add_argument(
        "-va",
        "--varank",
        action="store_true",
        help="additionnal argument to launch a preventive analysis with VaRank software if you want some annotation not generated with howard",
    )
    assembly_parser = argparse.ArgumentParser(add_help=False)
    assembly_parser.add_argument(
        "-a",
        "--assembly",
        required=True,
        type=str,
        help="argument to choose with which assembly you want to do the analysis",
    )
    output_format_parser = argparse.ArgumentParser(add_help=False)
    output_format_parser.add_argument(
        "-of",
        "--output_format",
        action="store_true",
        help="Choose output format, available VCF, Parquet, TSV, CSV, PSV or duckDB",
    )
    param_parser = argparse.ArgumentParser(add_help=False)
    param_parser.add_argument(
        "-pa",
        "--param",
        action="store_true",
        help="Parameters JSON file (or string) defines parameters to process annotations, calculations, prioritizations, convertions and queries.",
    )
    type_parser = argparse.ArgumentParser(add_help=False)
    type_parser.add_argument(
        "-t",
        "--type",
        action="store_true",
        help="Choose howard analysis type, annotation or calculation",
    )

    # Subparser definition
    subparsers = main_parser.add_subparsers(help="sub-command help")

    parser_folder = subparsers.add_parser(
        "folder",
        help="run variantannotation on any folder you want, must containing VCFs and optional custom configfile",
        parents=[
            verbosity_parser,
            varank_parser,
            assembly_parser,
            output_format_parser,
            param_parser,
            type_parser,
        ],
    )
    parser_folder.add_argument(
        "-f",
        "--folder",
        required=True,
        type=checker.absolute_folder_path,
        help="absolute path to the folder where you want to launch variantannotation analysis",
    )

    parser_run = subparsers.add_parser(
        "run",
        help="run variantannotation on STARK run repository",
        parents=[
            verbosity_parser,
            pattern_parser,
            mode_parser,
            varank_parser,
            assembly_parser,
            output_format_parser,
            param_parser,
            type_parser,
        ],
    )
    parser_run.add_argument(
        "-r",
        "--run",
        required=True,
        type=checker.absolute_run_path,
        help="absolute path to STARK run for which you want to launch variantannotation analysis",
    )

    parser_sync = subparsers.add_parser(
        "sync",
        help="used to rsync all VCF files in the described STARK repository/archives and defined pattern",
        parents=[verbosity_parser, pattern_parser, mode_parser],
    )
    sync_group = parser_sync.add_mutually_exclusive_group(required=True)
    sync_group.add_argument(
        "-sr",
        "--syncrepo",
        action="store_true",
    )
    sync_group.add_argument(
        "-sa",
        "--syncarchives",
        action="store_true",
    )

    args = main_parser.parse_args()

    if hasattr(args, "pattern"):
        if commons.default_pattern not in args.pattern:
            setattr(args, "pattern", [commons.default_pattern, args.pattern])
        elif args.pattern == commons.default_pattern:
            setattr(args, "pattern", [commons.default_pattern])

    if len(sys.argv) == 1:
        main_parser.print_help(sys.stderr)
        sys.exit(1)

    main(args)


if __name__ == "__main__":
    parse_args()