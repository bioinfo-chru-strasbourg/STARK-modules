from os.path import join as osj
import os
import subprocess
import logging as log
import json


def launch(container_name, launch_arguments):
    module_config = osj(
        os.environ["HOST_MODULE_CONFIG"],
        f"{os.environ["DOCKER_SUBMODULE_NAME"]}_config.json",
    )

    if not os.path.isfile(module_config):
        log.error(f"{module_config} do not exist, primordial file, check its existence")
        raise ValueError(module_config)

    with open(module_config, "r") as read_file:
        data = json.load(read_file)
        howard_image = data["howard_image"]
    log.info(f"Using {howard_image}")

    command_list = [
        "docker",
        "run",
        "--rm",
        "--name",
        container_name,
        "--env",
        f"http_proxy={os.environ["http_proxy"]}",
        "--env",
        f"ftp_proxy={os.environ["ftp_proxy"]}",
        "--env",
        f"https_proxy={os.environ["https_proxy"]}",
        "-v",
        f"{os.environ["HOST_TMP"]}:{os.environ["HOST_TMP"]}",
        "-v",
        f"{os.environ["HOST_DATABASES"]}:{os.environ["HOST_DATABASES"]}",
        "-v",
        f"{os.environ["HOST_SERVICES"]}:{os.environ["HOST_SERVICES"]}",
        "-v",
        f"{os.environ["HOST_CONFIG"]}:{os.environ["HOST_CONFIG"]}",
        "-v",
        "/var/run/docker.sock:/var/run/docker.sock",
        howard_image,
    ]

    command_list = command_list + launch_arguments
    log.debug(" ".join(command_list))
    print(" ".join(command_list))
    subprocess.call(command_list, universal_newlines=True)
