# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import shutil
import re
import pandas
import glob
import logging as log
import time

default_pattern = "*/STARK/*.reports/*.final.vcf.gz"


def set_log_level(verbosity):
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
        filename="",
        filemode="a",
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=configs[verbosity],
    )


if __name__ == "__main__":
    pass
