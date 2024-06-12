# -*- coding: utf-8 -*-
import os
import logging as log
from os.path import join as osj
from threading import local
import sys
import time

default_pattern = "*/STARK/*.reports/*.final.vcf.gz"


def set_logger_info_run(args):
    mylog = log.getLogger()
    log_file = mylog.handlers[0].baseFilename
    info_file = log_file.replace(".log", ".info")

    if not os.path.isfile(info_file):
        python_command = " ".join(sys.argv)
        python_version = f"{sys.version.split(' ')[0].split('.')[0]}.{sys.version.split(' ')[0].split('.')[1]}"
        command = f"python{python_version} {python_command}"
        global time_seconds_start
        time_seconds_start = time.time() + 7200
        local_time = time.localtime(time_seconds_start)
        actual_time = time.strftime("%a %b %d %H:%M:%S %Y", local_time)
        start = actual_time

        logging = f"""Command: {command}
Start time: {start}\n"""

        with open(info_file, "a") as write_file:
            write_file.write(logging)

    else:
        time_seconds_end = time.time() + 7200
        local_time = time.localtime(time_seconds_end)
        actual_time = time.strftime("%a %b %d %H:%M:%S %Y", local_time)
        end = actual_time
        delta = time_seconds_end - time_seconds_start

        logging = f"""End time: {end}
Time run: {delta}"""

        with open(info_file, "a") as write_file:
            write_file.write(logging)


def logger_header_run(log_file):
    logging = f"""#########################
#   Howard annotations  #
# Author: Mateusz Rauch #
#########################

####################
#   Release: {os.environ["DOCKER_SERVICE_CLI_RELEASE"]} #
####################

"""
    with open(log_file, "a") as write_file:
        write_file.write(logging)


def set_log_level_run(args):
    verbosity = args.verbosity
    run = args.run
    mode = args.launchmode
    if run.endswith("/"):
        run = run[:-1]

    run_name = run.split("/")[-1]
    run_application = run.split("/")[-2]
    run_platform = run.split("/")[-3]

    time_seconds = time.time() + 7200
    local_time = time.localtime(time_seconds)
    actual_time = time.strftime("%Y%m%d_%H%M%S", local_time)

    if mode == "manual":
        log_file = (
            f"{actual_time}_{run_platform}_{run_application}_{run_name}_{mode}.log"
        )
    elif mode == "listener":
        log_file = (
            f"{actual_time}_{run_platform}_{run_application}_{run_name}_{mode}.log"
        )

    log_file = osj(
        os.environ["DOCKER_SERVICES"],
        "logs",
        log_file,
    )

    logger_header_run(log_file)

    configs = {
        "debug": log.DEBUG,
        "info": log.INFO,
        "warning": log.WARNING,
        "error": log.ERROR,
        "critical": log.CRITICAL,
    }
    if verbosity not in configs.keys():
        raise ValueError(
            "Unknown verbosity level:"
            + verbosity
            + "\nPlease use any in:"
            + configs.keys()
        )

    log.basicConfig(
        filename=log_file,
        force=True,
        filemode="a",
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )


def set_log_level_default(args):
    verbosity = args.verbosity

    configs = {
        "debug": log.DEBUG,
        "info": log.INFO,
        "warning": log.WARNING,
        "error": log.ERROR,
        "critical": log.CRITICAL,
    }
    if verbosity not in configs.keys():
        raise ValueError(
            "Unknown verbosity level:"
            + verbosity
            + "\nPlease use any in:"
            + configs.keys()
        )
    log.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
        level=configs[verbosity],
    )


if __name__ == "__main__":
    pass
