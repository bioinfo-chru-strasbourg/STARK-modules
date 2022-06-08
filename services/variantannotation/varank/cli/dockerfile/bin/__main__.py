# !/usr/bin/python
# !/bin/sh
# !/usr/tcl/bin
# -*- coding: utf-8 -*-
"""
###############################
##							 ##
##		VARANK ANALYSIS		 ##
## need alamut-batch license ##
##	Author : Mateusz RAUCH	 ##
##							 ##
###############################
"""

import argparse
from locale import strcoll
import sys
import os
import logging as log
import commons
import checker
import run_processing


def main(args):
    original_umask = os.umask(0o000)
    commons.set_log_level(args)
    checker.varank_config_json_checker()

    # mylog = log.getLogger()
    # print(mylog.handlers[0].baseFilename)

    if args.run:
        run_processing.launch_run(args)

    # elif args.varank_folder:
    # 	launch_folder_analysis(args)
    # elif args.allsync18:
    # 	pattern = pattern_generator(args)
    # 	launch_universal_sync18(pattern)
    os.umask(original_umask)


def parse_args():

    # Main parser
    main_parser = argparse.ArgumentParser(prog="varank")

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
        "--mode",
        type=str,
        default="manual",
        help="hidden argument to know if the varank analysis is manual or listener launched",
    )

    # Subparser definition
    subparsers = main_parser.add_subparsers(help="sub-command help")

    parser_folder = subparsers.add_parser(
        "folder",
        help="run VaRank on any folder you want, must containing VCFs and optional custom configfile",
        parents=[verbosity_parser],
    )
    parser_folder.add_argument(
        "-f",
        "--folder",
        required=True,
        type=checker.absolute_folder_path,
        help="absolute path to the folder where you want to launch VaRank analysis",
    )

    parser_run = subparsers.add_parser(
        "run",
        help="run VaRank on STARK run repository",
        parents=[verbosity_parser, pattern_parser, mode_parser],
    )
    parser_run.add_argument(
        "-r",
        "--run",
        required=True,
        type=checker.absolute_run_path,
        help="absolute path to STARK run for which you want to launch VaRank analysis",
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
