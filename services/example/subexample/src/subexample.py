"""
@Goal: Provide a minimal example module

@Author: Samuel Nicaise (2026)
"""

import argparse
import logging as log
import os
from os.path import join as osj
import time

def set_log_level(verbosity):
    verbosity = verbosity.lower()
    configs = {
        "debug": log.DEBUG,
        "info": log.INFO,
        "warning": log.WARNING,
        "error": log.ERROR,
        "critical": log.CRITICAL,
    }
    if verbosity not in configs.keys():
        raise ValueError(
            f"Unknown verbosity level: {verbosity}\nPlease use any in: {configs.keys()}"
        )
    log.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )


def main(run_dir: str, genome: str, verbosity: str = "info") -> None:
    set_log_level(verbosity)
    log.info(f"run_dir: {run_dir}")
    log.info(f"genome: {genome}")

    ###############################
    # long duration task example so you can test canceling it
    def fib(n):
        a, b = 0, 1
        out = []
        for _ in range(n):
            out.append(a)
            a, b = b, a + b
        return out
    res = fib(10000000)
    log.info(f"fib done {res[0]}")
    ###############################

    running_file = osj(run_dir,"SUBEXAMPLERunning.txt")
    complete_file = osj(run_dir, "SUBEXAMPLEComplete.txt")
    with open(complete_file, "w") as f:
        f.write(time.ctime())
    if os.path.exists(running_file):
        os.remove(running_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="subexample to launch routine analysis")
    parser.set_defaults(mode=main)
    parser.add_argument(
        "-i",
        "--runDir",
        type=str,
        help="path to run in a STARK 0.9.18 repository",
        required=True,
    )
    parser.add_argument(
        "-g", "--genome", help="genome file", type=str, dest="genome", required=True
    )
    parser.add_argument(
        "-v",
        "--verbosity",
        help="set the logging level (debug, info, warning, error, critical)",
        type=str,
        default="info",
    )
    args = parser.parse_args()
    if not hasattr(args, "mode"):
        parser.print_help()
    else:
        args.mode(args.runDir, args.genome, verbosity=args.verbosity)