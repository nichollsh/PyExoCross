#!/usr/bin/env python

import os
import sys
import multiprocessing as mp
import concurrent.futures as cf

project_root = os.path.dirname(os.path.abspath(__file__))
src_root = os.path.join(project_root, 'src')
if src_root not in sys.path:
    sys.path.insert(0, src_root)

try:
    mp.set_start_method('fork')
except (RuntimeError, ValueError):
    pass

# macOS stability: avoid nested process-pool deadlocks in long workflows.
if sys.platform == 'darwin':
    cf.ProcessPoolExecutor = cf.ThreadPoolExecutor

from pyexocross.config import Config
from pyexocross.core import get_results
from pyexocross.base.log import (
    output_context,
    parse_logging_info,
    parse_verbose_info,
    setup_logging,
)
from pyexocross.base.input import parse_args

if __name__ == '__main__':
    # get input file path and resume flag from command line
    args = parse_args()
    inp_path = args.path

    # load config from input file
    cfg = Config(inp_filepath=inp_path, force_reload=True)
    cfg.resume = args.resume

    # setup logging
    verbose = parse_verbose_info(inp_path)
    logpath = parse_logging_info(inp_path)
    if logpath is not None:
        setup_logging(logpath, announce=verbose)

    # run the main function with output context
    with output_context(verbose):
        get_results(cfg)
