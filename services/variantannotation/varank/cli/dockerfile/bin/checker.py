# -*- coding: utf-8 -*-
import os
import logging as log
from commons import set_log_level


def absolute_dir_path(path):
    if os.path.isabs(path) and os.path.isdir(path):
        return path
    elif not os.path.isabs(path) and os.path.isdir(path):
        return os.path.abspath(path)
    else:
        raise ValueError(path)
