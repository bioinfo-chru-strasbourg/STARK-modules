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

# tools = {
#    "java": "java",
#    "gatk": "/STARK/tools/gatk/3.8-1-0/bin/GenomeAnalysisTK.jar",
#    "VaRank": "/home1/TOOLS/tools/varank/VaRank_1.4.3",
#    "Alamut-batch": "/home1/TOOLS/tools/alamut_batch/alamut-batch-standalone-1.11",
#    "howard": "/STARK/tools/howard/current/bin/HOWARD",
#    "bcftools": "/STARK/tools/bcftools/current/bin/bcftools",
#    "bgzip": "/STARK/tools/htslib/current/bin/bgzip",
#    "tabix": "/STARK/tools/htslib/current/bin/tabix",
#    "threads": 6,
# }


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
        "-c", "--config", type=str, help="Absolute path of <config>.json created"
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
        type=dict,
        default="/app/config/default.json",
        help="Dict containing path of tools executable",
    )
    parser.add_argument(
        "-t", "--trio", type=str, help="List of sample representing example s1,s2,s3"
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
        action="store_true", help="if set process in STARK auto analysis"
    )
    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    if args.auto:
        # args.config = "/app/config/config.json"
        cf = config(
            args.run,
            args.trio,
            args.genome,
            args.tools,
            args.output,
            args.repository,
            args.depository,
            args.config,
        )
        subprocess.call(
            "python " + osj(os.getenv("MICROSERVICE_CONFIG"), "launcher.py")
        )


if __name__ == "__main__":
    args = parseargs()
