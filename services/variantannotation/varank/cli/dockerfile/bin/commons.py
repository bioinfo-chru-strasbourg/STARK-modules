# -*- coding: utf-8 -*-
import os
import logging as log
from os.path import join as osj
from threading import local
import time

default_pattern = "*/STARK/*.reports/*.final.vcf.gz"


# def logger_infos(args, log_file, seconds_start):
#     info_file = log_file.replace(".log", ".info")
#     if not os.path.isfile(info_file):
#         # command =
#         time_seconds_start = time.time() + 7200
#         local_time = time.localtime(time_seconds_start)
#         actual_time = time.strftime("%a %b %d %H:%M:%S %Y", local_time)
#         start = actual_time
#         print(start)

#         log = f"""
#         """
#         with open(info_file, "a") as write_file:
#             write_file.write(log)

#         return info_file, time_seconds_start

#     else:
#         time_seconds_end = time.time() + 7200
#         local_time = time.localtime(time_seconds_end)
#         actual_time = time.strftime("%a %b %d %H:%M:%S %Y", local_time)
#         end = actual_time
#         delta = time_seconds_end - seconds_start
#         print(delta)


def logger_header(log_file):
    log = f"""#########################
#     VaRank Wrapper    #
#    Variants Ranking   #
# Author: Mateusz Rauch #
#########################

####################
#   Release: {os.environ["DOCKER_STARK_MODULE_SUBMODULE_SERVICE_CLI_RELEASE"]}   #
####################

"""
    with open(log_file, "a") as write_file:
        write_file.write(log)


def set_log_level(args):
    verbosity = args.verbosity
    run = args.run
    mode = args.mode

    run_name = run.split("/")[-2]
    run_application = run.split("/")[-3]
    run_platform = run.split("/")[-4]

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
        os.environ["DOCKER_MODULE_VARIANTANNOTATION_SUBMODULE_VARANK_FOLDER_SERVICES"],
        "logs",
        log_file,
    )

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
        filemode="a",
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )

    # mylog = log.getLogger()
    # print(mylog.handlers[0].baseFilename)


if __name__ == "__main__":
    pass
