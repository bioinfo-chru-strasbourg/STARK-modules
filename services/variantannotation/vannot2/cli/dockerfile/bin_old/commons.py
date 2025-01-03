import os
import logging as log
from os.path import join as osj
from threading import local
import sys
import time
import json

def get_threads(threads_type):
    module_config = osj(os.environ["HOST_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    if not os.path.isfile(module_config):
        log.error(f"{module_config} do not exist, primordial file, check its existence")
        raise ValueError(module_config)
      
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        threads = data[threads_type]

    return threads

def get_memory(memory_type):
    module_config = osj(os.environ["HOST_MODULE_CONFIG"], f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json")
    if not os.path.isfile(module_config):
        log.error(f"{module_config} do not exist, primordial file, check its existence")
        raise ValueError(module_config)
      
    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        memory = data[memory_type]

    return memory

def get_default_pattern():
    default_pattern = "*/STARK/*.reports/*.final.vcf.gz"
    return default_pattern


def set_logger_info():
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


def logger_header(log_file):
    logging = f"""#########################
#         vAnnot        #
# Author: Mateusz Rauch #
#########################

####################
#   Release: {os.environ["DOCKER_SERVICE_CLI_RELEASE"]} #
####################

"""
    with open(log_file, "a") as write_file:
        write_file.write(logging)


def set_log_level(args):
    verbosity = args.verbosity
    time_seconds = time.time() + 7200
    local_time = time.localtime(time_seconds)
    actual_time = time.strftime("%Y%m%d_%H%M%S", local_time)

    if "run" in args:
        run = args.run
        if run.endswith("/"):
            run = run[:-1]
        mode = args.launchmode
        run_name = run.split("/")[-1]
        run_application = run.split("/")[-2]
        run_platform = run.split("/")[-3]
        log_file = f"{actual_time}_{run_platform}_{run_application}_{run_name}_{mode}.log"

    elif "folder" in args:
        folder = args.folder
        if folder.endswith("/"):
            folder = folder[:-1]
        folder_name = folder.split("/")[-1]
        log_file = f"{actual_time}_{folder_name}.log"
        
    elif "dejavu" in args:
        dejavu = args.dejavu
        if dejavu.endswith("/"):
            dejavu = dejavu[:-1]
        dejavu_name = dejavu.split("/")[-1]
        log_file = f"{actual_time}_dejavuonly_{dejavu_name}.log"

    log_file = osj(
        os.environ["HOST_SERVICES"],
        "logs",
        log_file,
    )

    logger_header(log_file)

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
        format="vAnnot %(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )


if __name__ == "__main__":
    pass
