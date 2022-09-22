"""
Wrapper for DPNI module associated with STARK software https://github.com/bioinfo-chru-strasbourg/STARK
#TODO
"""

# sys.path.append("..")
import subprocess
from configure import config
import argparse
import os
from os.path import join as osj
import json


def check_exists(item):
    assert os.path.exists(item), f"ERROR {item} does not exists"


def parseargs():
    """
    parse arguments
    """
    parser = argparse.ArgumentParser(
        # formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-r", "--run", type=str, help="Absolute path of run folder repository format"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="/app/res",
        help="Absolute path of output folder",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default="/app/config/analysis.json",
        help="Absolute path of <config>.json created",
    )
    parser.add_argument(
        "-g",
        "--genome",
        type=str,
        default="/STARK/databases/genomes/current/hg19.fa",
        help="Absolute path of genome assemby fasta file",
    )
    parser.add_argument(
        "--tools",
        type=str,
        default="/app/config/default.json",
        help="json containing path of tools executable",
    )
    parser.add_argument(
        "-t",
        "--trio",
        type=str,
        default="",
        help="List of sample representing example s1,s2,s3",
    )
    parser.add_argument(
        "-ry",
        "--repository",
        type=str,
        default="",
        help="Absolute path of STARK repository",
    )
    parser.add_argument(
        "-dy",
        "--depository",
        type=str,
        default="",
        help="Absolute path of STARK depository",
    )
    parser.add_argument(
        "--launcher", action="store_true", help="if set process in STARK auto analysis"
    )
    parser.add_argument(
        "--jobs", type=int, default=4, help="Jobs to run in parallel in snakemake"
    )
    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    print("#[INFO] DPNI wrapper")
    print(args)
    for folder in [args.run, args.repository, args.depository, args.genome]:
        if folder:
            check_exists(folder)
    with open(args.tools) as too:
        tools = json.load(too)
    # if args.launcher:
    #    cf = config(
    #        args.run,
    #        args.trio,
    #        args.genome,
    #        tools,
    #        args.output,
    #        args.repository,
    #        args.depository,
    #        args.config,
    #    )
    #    # TODO
    #    subprocess.call(
    #        "python " + osj(os.getenv("MICROSERVICE_CONFIG"), "launcher.py")
    #    )
    # else:
    cf = config(
        args.run,
        args.trio,
        args.genome,
        tools,
        args.output,
        args.repository,
        args.depository,
        args.config,
    )
    cf.configure()
    subprocess.call(
        "snakemake --snakefile /app/lib/snakemake/Snakefile --configfile "
        + args.config
        + " -j "
        + args.jobs
    )


if __name__ == "__main__":
    main()
