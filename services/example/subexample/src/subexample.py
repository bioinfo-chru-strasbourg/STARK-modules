"""
@Goal: Provide a minimal example module

@Author: Samuel Nicaise (2026)
"""

import argparse
import logging as log

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


def main(run_dir: str, genome: str) -> None:
    print("run_dir:", run_dir)
    print("genome:", genome)

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
    args = parser.parse_args()
    if not hasattr(args, "mode"):
        parser.print_help()
    else:
        args.mode(args.runDir, args.genome)