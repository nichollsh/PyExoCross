"""
Utilities for handling large transition files.

This module provides functions for efficiently processing large transition files
that may not fit in memory, using chunked reading and parallel processing.
"""
import os
from pathlib import Path
import bz2
import subprocess
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED

TMPDIR = os.getenv('TMPDIR', '/tmp')

from .constants import (
    DEFAULT_CHUNK_SIZE,
    DEFAULT_CHUNK_CACHE_GB,
    DEFAULT_MAX_MEMORY_GB,
    DEFAULT_LARGE_FILE_GB,
    LARGE_TRANS_FILE_BYTES,
    MAX_LARGE_FILE_WORKERS,
    MAX_INFLIGHT_MULTIPLIER,
    LARGE_WRITE_CHUNK_ROWS,
)
from .utils import ensure_dir
from .log import log_tqdm


@dataclass(frozen=True)
class TransSource:
    """Reusable transition source backed by a file or an in-memory DataFrame."""

    path: str
    origin: str
    frame: object = field(default=None, compare=False, hash=False, repr=False)
    size: int = 0
    minwn: object = None
    maxwn: object = None


# Run-time memory settings (set per run from LargeFileSize / ChunkCacheMemory by set_memory_limits).
_large_file_bytes = LARGE_TRANS_FILE_BYTES
_chunk_cache_bytes = int(DEFAULT_CHUNK_CACHE_GB * 1e9)    # None means unlimited


def normalize_large_file_size(value):
    """Return the large-file threshold in GB (decimal) as a positive float."""
    if value is None or (isinstance(value, str) and value.strip().upper() in ('', 'AUTO', 'NONE')):
        return float(DEFAULT_LARGE_FILE_GB)
    if isinstance(value, float) and np.isnan(value):
        return float(DEFAULT_LARGE_FILE_GB)
    try:
        size = float(value)
    except (TypeError, ValueError):
        raise ValueError('LargeFileSize / large_file_size must be a positive number of GB.')
    if not size > 0:
        raise ValueError('LargeFileSize / large_file_size must be a positive number of GB.')
    return size


def normalize_chunk_cache_memory(value):
    """Return the chunk-cache budget in GB (decimal) as a float >= 0, or None for unlimited."""
    if value is None:
        return None
    if isinstance(value, float) and np.isnan(value):
        return None
    if isinstance(value, str):
        text = value.strip().upper()
        if text in ('NONE', 'UNLIMITED', 'INF'):
            return None
        if text in ('', 'AUTO'):
            return float(DEFAULT_CHUNK_CACHE_GB)
    try:
        size = float(value)
    except (TypeError, ValueError):
        raise ValueError('ChunkCacheMemory / chunk_cache_memory must be a non-negative number of GB or UNLIMITED.')
    if not size >= 0:
        raise ValueError('ChunkCacheMemory / chunk_cache_memory must be a non-negative number of GB or UNLIMITED.')
    return size


def normalize_max_memory(value):
    """Return the Cache=Auto in-memory limit in GB (decimal) as a float >= 0."""
    if value is None or (isinstance(value, str) and value.strip().upper() in ('', 'AUTO')):
        return float(DEFAULT_MAX_MEMORY_GB)
    if isinstance(value, float) and np.isnan(value):
        return float(DEFAULT_MAX_MEMORY_GB)
    try:
        size = float(value)
    except (TypeError, ValueError):
        raise ValueError('MaxMemory / max_memory must be zero or a positive number of GB.')
    if not size >= 0:
        raise ValueError('MaxMemory / max_memory must be zero or a positive number of GB.')
    return size


def set_memory_limits(large_file_gb=None, chunk_cache_gb=DEFAULT_CHUNK_CACHE_GB):
    """Apply the large-file threshold and chunk-cache budget, both in GB (None cache = unlimited)."""
    global _large_file_bytes, _chunk_cache_bytes
    _large_file_bytes = int(normalize_large_file_size(large_file_gb) * 1e9)
    gb = normalize_chunk_cache_memory(chunk_cache_gb)
    _chunk_cache_bytes = None if gb is None else int(gb * 1e9)


def large_file_bytes():
    """Return the size (bytes) at or above which a transition file is always streamed."""
    return _large_file_bytes


def chunk_cache_bytes():
    """Return the total budget (bytes) for retained parsed chunks, or None if unlimited."""
    return _chunk_cache_bytes


def sourcepath(source):
    """Return the physical or identifying path for a transition source."""
    return source.path if isinstance(source, TransSource) else str(source)


def sourcename(source):
    """Return a display name for a transition source."""
    return os.path.basename(sourcepath(source))


def is_large_trans_file(trans_filepath):
    """
    Check if a transition file exceeds the large file size threshold.

    Parameters
    ----------
    trans_filepath : str
        Path to the transition file to check.

    Returns
    -------
    bool
        True if the file size is greater than or equal to the large-file threshold (LargeFileSize),
        False otherwise or if an error occurs accessing the file.
    """
    if isinstance(trans_filepath, TransSource):
        return trans_filepath.size >= _large_file_bytes
    if not os.path.exists(trans_filepath):
        return False
    return os.path.getsize(trans_filepath) >= _large_file_bytes

# Decompress Large .trans.bz2 Files
def command_decompress(trans_filename, to_tmpdir=True):
    """
    Decompress a bzip2-compressed transition file.

    Decompresses the file to a dedicated decompressed directory if it doesn't
    already exist. Uses system bunzip2 command for decompression.

    Parameters
    ----------
    trans_filename : str
        Path to the compressed .trans.bz2 file to decompress.

    Returns
    -------
    tuple of (str, int)
        A tuple containing:
        - Path to the decompressed file
        - Flag indicating if decompression was performed (1) or file already existed (0)
    """
    # Directory where the decompressed .trans files will be saved
    if to_tmpdir:
        trans_dir = Path(TMPDIR) / Path(trans_filename).name / 'decompressed'
    else:
        trans_dir = Path(trans_filename).parent / 'decompressed'
    ensure_dir(str(trans_dir))
    trans_file = os.path.join(trans_dir, trans_filename.split('/')[-1].replace('.bz2', ''))
    if os.path.exists(trans_file):
        num = 0
    else:
        command = f'bunzip2 < {trans_filename} > {trans_file}'
        print('Decompressing', trans_filename.split('/')[-1], '>', trans_file)
        subprocess.run(command, shell=True)
        num = 1
    return(trans_file, num)

def read_trans_chunks(trans_filepath, usecols, names, chunk_sz=None, extracols=None):
    """
    Read a transition file in chunks for memory-efficient processing.

    Reads the file as an iterator of DataFrame chunks, automatically handling
    bzip2 compression if the file extension indicates it.

    Parameters
    ----------
    trans_filepath : str
        Path to the transition file to read.
    usecols : list of int
        Column indices to read from the file.
    names : list of str
        Column names to assign to the selected columns.
    chunk_sz : int, optional
        Number of rows per chunk. If None, uses the global chunk_size variable
        or DEFAULT_CHUNK_SIZE as fallback.
    extracols : list of str, optional
        Named state-derived columns to project from a Parquet transition cache.

    Returns
    -------
    TextFileReader
        A pandas TextFileReader iterator that yields DataFrame chunks.
    """
    if chunk_sz is None:
        import pyexocross.core as core_module
        chunk_sz = getattr(core_module, 'chunk_size', DEFAULT_CHUNK_SIZE)
    if isinstance(trans_filepath, TransSource) and trans_filepath.frame is not None:
        frame = trans_filepath.frame

        def memorychunks():
            for start in range(0, len(frame), chunk_sz):
                chunk = frame.iloc[start:start + chunk_sz, usecols].copy()
                chunk.columns = names
                yield chunk

        return memorychunks()

    path = sourcepath(trans_filepath)
    if path.endswith('.parquet'):
        import pyarrow.parquet as pq

        columns = [f'c{col}' for col in usecols]

        def parquetchunks():
            parquet_file = pq.ParquetFile(path)
            available = set(parquet_file.schema_arrow.names)
            extras = [
                name for name in (extracols or ())
                if name in available and name not in columns
            ]
            readcolumns = columns + extras
            for batch in parquet_file.iter_batches(batch_size=chunk_sz, columns=readcolumns):
                chunk = batch.to_pandas()
                chunk = chunk.rename(columns=dict(zip(columns, names)))
                yield chunk

        return parquetchunks()

    read_kwargs = dict(sep=r'\s+', header=None, usecols=usecols, names=names,
                       chunksize=chunk_sz, iterator=True, low_memory=False, encoding='utf-8')
    if path.endswith('.bz2'):
        read_kwargs['compression'] = 'bz2'
    return pd.read_csv(path, **read_kwargs)

def chunks_nbytes(chunk):
    """Return the in-memory size of one parsed chunk in bytes."""
    return int(chunk.memory_usage(index=True, deep=True).sum())


def collect_chunks_within_budget(reader, budget_bytes):
    """
    Materialize chunks from a reader, giving up as soon as the budget is exceeded.

    Returns (chunks, nbytes), or (None, nbytes) if the budget was exceeded, in which
    case the partial chunks are dropped so peak memory stays near the budget.
    budget_bytes=None means unlimited.
    """
    chunks = []
    nbytes = 0
    for chunk in reader:
        nbytes += chunks_nbytes(chunk)
        if budget_bytes is not None and nbytes > budget_bytes:
            close = getattr(reader, 'close', None)
            if close is not None:
                close()
            return None, nbytes
        chunks.append(chunk)
    return chunks, nbytes


def cache_small_file_chunks(trans_filepaths, use_cols, use_names, extracols=None, budget_bytes='default'):
    """
    Materialize and cache parsed chunks once for small files, within a total memory budget.

    Reads each non-large transition file via read_trans_chunks and stores the resulting
    chunk DataFrames, so callers can reuse the same parsed data across a whole
    temperature/pressure grid instead of re-reading the file per grid point. Files at or
    above the large-file threshold (LargeFileSize), and files that would push the total
    retained size over the budget (ChunkCacheMemory), are left out of the returned dict;
    callers must stream those lazily so memory usage stays bounded.

    Parameters
    ----------
    trans_filepaths : iterable
        Transition file paths or TransSource objects to consider caching.
    use_cols : list of int
        Column indices to read, passed through to read_trans_chunks.
    use_names : list of str
        Column names to assign, passed through to read_trans_chunks.
    extracols : list of str, optional
        Extra state-derived columns to project from a Parquet transition cache.
    budget_bytes : int or None, optional
        Total bytes of parsed chunks to retain across all files; None is unlimited.
        By default the run's ChunkCacheMemory setting is used.

    Returns
    -------
    dict
        Mapping from trans_filepath to a list of already-parsed chunk DataFrames,
        containing only entries for files that were cached. Files absent from the
        mapping must be streamed.
    """
    if isinstance(budget_bytes, str):
        budget_bytes = _chunk_cache_bytes
    cache = {}
    used = 0
    streamed = []
    for trans_filepath in trans_filepaths:
        if is_large_trans_file(trans_filepath):
            continue
        remaining = None if budget_bytes is None else budget_bytes - used
        if remaining is not None and remaining <= 0:
            streamed.append(sourcename(trans_filepath))
            continue
        chunks, nbytes = collect_chunks_within_budget(
            read_trans_chunks(trans_filepath, use_cols, use_names, extracols=extracols),
            remaining,
        )
        if chunks is None:
            streamed.append(sourcename(trans_filepath))
        else:
            cache[trans_filepath] = chunks
            used += nbytes
    if budget_bytes is not None:
        print(f'Chunk cache: holding {used / 1e9:.2f} GB of {budget_bytes / 1e9:.2f} GB budget '
              f'across {len(cache)} file(s); streaming {len(streamed)} file(s) that do not fit.')
    return cache


def process_large_chunks(trans_reader, handler, combine_fn, zero_factory, desc,
                         max_workers=MAX_LARGE_FILE_WORKERS, max_inflight=None, reducer=None,
                         executor_class=ProcessPoolExecutor):
    """
    Process large transition file chunks in parallel with bounded memory usage.

    Processes chunks from a file reader using a process pool executor, with
    controlled parallelism to limit memory consumption. Supports both accumulation
    and collection modes for combining results.

    Parameters
    ----------
    trans_reader : iterable
        Iterator that yields DataFrame chunks to process.
    handler : callable
        Function to process each chunk. Should accept a DataFrame chunk and return a result.
    combine_fn : callable
        Function to combine multiple results into a single result (used when reducer is None).
    zero_factory : callable
        Function that returns an empty result of the appropriate type.
    desc : str
        Description string for progress bar display.
    max_workers : int, optional
        Maximum number of worker processes. Default is MAX_LARGE_FILE_WORKERS.
    max_inflight : int, optional
        Maximum number of concurrent futures. If None, calculated as
        max_workers * MAX_INFLIGHT_MULTIPLIER.
    reducer : callable, optional
        Function to reduce two results into one. If provided, uses accumulation
        mode instead of collection mode. Should accept (accumulator, new_result)
        and return combined result.

    Returns
    -------
    object
        Combined result from all chunks, or empty result from zero_factory if no chunks processed.
    """
    if max_inflight is None:
        max_inflight = max_workers * MAX_INFLIGHT_MULTIPLIER
    results = [] if reducer is None else None
    accumulator = None
    with executor_class(max_workers=max_workers) as executor:
        futures = set()
        for trans_df_chunk in log_tqdm(trans_reader, desc=desc):
            futures.add(executor.submit(handler, trans_df_chunk))
            if len(futures) >= max_inflight:
                done, futures = wait(futures, return_when=FIRST_COMPLETED)
                for done_future in done:
                    result = done_future.result()
                    if reducer is None:
                        results.append(result)
                    else:
                        accumulator = result if accumulator is None else reducer(accumulator, result)
        for future in as_completed(futures):
            result = future.result()
            if reducer is None:
                results.append(result)
            else:
                accumulator = result if accumulator is None else reducer(accumulator, result)
    if reducer is None:
        return combine_fn(results) if results else zero_factory()
    return accumulator if accumulator is not None else zero_factory()

def _prepare_array_writer(data):
    """
    Prepare data for chunked writing by creating a row getter function.

    Converts input data (DataFrame or array-like) into a format suitable for
    chunked writing, returning the total number of rows and a getter function.

    Parameters
    ----------
    data : pd.DataFrame or array-like
        Data to prepare for writing. Can be a DataFrame or any array-like structure.

    Returns
    -------
    tuple of (int, callable)
        A tuple containing:
        - Total number of rows in the data
        - Getter function that accepts (start, end) indices and returns a numpy array slice
    """
    if isinstance(data, pd.DataFrame):
        total_rows = len(data)
        def getter(start, end):
            return data.iloc[start:end].to_numpy()
    else:
        arr = np.asarray(data)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        total_rows = arr.shape[0]
        def getter(start, end):
            return arr[start:end]
    return total_rows, getter

def save_large_txt(filepath_or_buffer, data, fmt, header=None, comments='', chunk_rows=LARGE_WRITE_CHUNK_ROWS, mode='w'):
    """
    Save large arrays or DataFrames to text file in chunks to manage memory.

    Writes data to a file in chunks to avoid loading the entire dataset into memory.
    Supports both regular files and bzip2-compressed files. Can write to a file path
    or an open file-like object.

    Parameters
    ----------
    filepath_or_buffer : str or file-like object
        Path to output file or an open file handle. If string ends with '.bz2',
        file will be compressed with bzip2.
    data : pd.DataFrame or array-like
        Data to write. Can be a DataFrame or any array-like structure.
    fmt : str or sequence of str
        Format string(s) for numpy.savetxt formatting.
    header : str, optional
        Header string to write at the beginning of the file. Default is None.
    comments : str, optional
        Comment prefix for header line. Default is empty string.
    chunk_rows : int, optional
        Number of rows to write per chunk. Default is LARGE_WRITE_CHUNK_ROWS.
    mode : str, optional
        File open mode ('w' for write, 'a' for append). Default is 'w'.
    """
    total_rows, getter = _prepare_array_writer(data)
    if hasattr(filepath_or_buffer, 'write'):
        f = filepath_or_buffer
        close_after = False
    else:
        # Handle bz2 compressed files
        if isinstance(filepath_or_buffer, str) and filepath_or_buffer.endswith('.bz2'):
            f = bz2.open(filepath_or_buffer, mode + 't')  # 't' for text mode
        else:
            f = open(filepath_or_buffer, mode)
        close_after = True
    try:
        if header is not None and header != '':
            f.write((comments or '') + header + '\n')
        if total_rows == 0:
            np.savetxt(f, np.empty((0, 1)), fmt=fmt)
            return
        for start in range(0, total_rows, chunk_rows):
            end = min(start + chunk_rows, total_rows)
            chunk = getter(start, end)
            if chunk.size == 0:
                continue
            np.savetxt(f, chunk, fmt=fmt)
    finally:
        if close_after:
            f.close()
