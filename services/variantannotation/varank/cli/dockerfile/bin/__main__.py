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
from commons import set_log_level


def main(args):
    set_log_level(args.verbosity)
    if args.inputFormat.lower() == "decon":
        raise ValueError(
            "DECON is handled as a TSV conversion. Use 'tsv' as input format"
        )

    factory = ConverterFactory()
    converter = factory.get_converter(
        args.inputFormat.lower(), args.outputFormat.lower(), args.configFile
    )

    if args.inputFormat == "varank":
        if args.coordConversionFile == "":
            raise ValueError(
                "Converting from a Varank file requires setting the --coordConversionFile argument to an existing file"
            )
        if not os.path.exists(args.coordConversionFile):
            raise ValueError(
                "coordConversionFile does not exist:" + args.coordConversionFile
            )
        converter.set_coord_conversion_file(args.coordConversionFile)

    converter.convert(args.inputFile, args.outputFile)


def parse_args():

    parser = argparse.ArgumentParser(prog="VaRank")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(
        "-f",
        "--varank_folder",
        type=str,
        help="absolute path to the folder where you want to launch VaRank analysis",
    )
    group.add_argument(
        "-r",
        "--run",
        type=str,
        help="path to STARK run for which you want to launch VaRank analysis",
    )
    group.add_argument(
        "-s",
        "--vcfsync",
        help="Used to rsync all VCF files in the described STARK repository and defined pattern",
        action="store_true",
    )
    parser.add_argument(
        "-p",
        "--pattern",
        type=str,
        nargs="+",
        help="pattern describing which vcf files to synchronize in STARK folders, actual patterns : */STARK/*.reports/*.final.vcf.gz, */POOL/*.final.vcf.gz. You can use two patterns with space as separator",
    )
    parser.add_argument(
        "-v", "--verbosity", type=str, default="info", help="Verbosity level"
    )

    if "verbosity" not in args:
        parser.print_help()

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args
    main(args)
