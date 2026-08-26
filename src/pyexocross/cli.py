"""
Command-line interface for PyExoCross.
"""
import argparse
import multiprocessing as mp
import concurrent.futures as cf
import sys
from pyexocross.config import Config
from pyexocross.core import get_results
from pyexocross.base.log import (
    output_context,
    parse_logging_info,
    parse_verbose_info,
    setup_logging,
)


def main():
    """
    Main CLI entry point.
    """
    parser = argparse.ArgumentParser(
        description='PyExoCross: Calculate molecular cross sections and spectra',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using .inp file
  pyexocross -p config.inp
  
  # Using Python API
  python -c "import pyexocross as pe; pe.stick_spectra(database='ExoMol', ...)"
        """
    )
    parser.add_argument(
        '-p', '--path',
        type=str,
        metavar='PATH',
        required=True,
        help='Path to .inp configuration file'
    )
    parser.add_argument(
        '-r', '--resume',
        action='store_true',
        default=False,
        help=(
            'Resume an interrupted run; skip existing T-P grid points.'
        ),
    )

    args = parser.parse_args()
    try:
        mp.set_start_method('fork')
    except (RuntimeError, ValueError):
        pass
    if sys.platform == 'darwin':
        cf.ProcessPoolExecutor = cf.ThreadPoolExecutor
    
    try:
        verbose = parse_verbose_info(args.path)
        logpath = parse_logging_info(args.path)
        if logpath is not None:
            setup_logging(logpath, announce=verbose)
        with output_context(verbose):
            config = Config(inp_filepath=args.path)
            config.resume = args.resume
            get_results(config)
            print('Done!')
    except Exception as e:
        print(f'Error: {e}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
