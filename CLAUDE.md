# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

PyExoCross is a Python framework for LTE/non-LTE spectroscopic post-processing of molecular and
atomic line lists (ExoMol, ExoMolHR, ExoAtom, HITRAN, HITEMP). It computes cross sections, stick
spectra, partition functions, specific heats, cooling functions, lifetimes, and oscillator
strengths, with optional CUDA/MPS GPU acceleration, and converts line lists between the ExoMol and
HITRAN formats. Fortran predecessor: ExoCross. Published as `pyexocross` on PyPI; source package
lives under `src/pyexocross/`.

Two calling conventions produce identical results:
- `.inp` file driven (legacy/CLI): `python run.py -p input/H2O_ExoMol_xsec.inp` or `pyexocross -p ...`
- Python API: `import pyexocross as px; px.cross_sections(database=..., molecule=..., ...)`

`docs/source/python_api/inp_mapping.md` is the authoritative table mapping every `.inp` keyword to
its Python kwarg — consult it rather than guessing kwarg names.

## Commands

```bash
# Install (editable, with dev deps for pytest/black/flake8)
pip install -e ".[dev]"

# Run via .inp file
python run.py -p ./input/H2O_ExoMol_xsec.inp
# or, after install:
pyexocross -p ./input/H2O_ExoMol_xsec.inp

# Background run (see docs/source/command.md for details)
nohup python -u run.py -p ./input/H2O_ExoMol_xsec.inp > /dev/null 2>&1 &
```

There is no pytest suite in the CI sense. Files under `test/` (`test_api_*.py`, `test_download.py`)
are example/smoke scripts: their `test_*` functions call the public API and print/save output but
contain no `assert`s, and `COMMON` dicts at the top of each file hardcode absolute
`read_path`/`save_path`/`logs_path` values from the original author's machine. To use one, edit
those paths first, then run directly (`python test/test_api_exomol.py`), not `pytest` — invoking
them under `pytest` unmodified will attempt real filesystem I/O against paths that don't exist here.
`slurm.sh` is an HPC (Slurm) submission template for large runs; `AGNI` and `SOCRATES` at the repo
root are symlinks to sibling projects, not part of this package.

## Architecture

**Config resolution (`config.py`, `base/config_manager.py`, `base/input.py`).** A `Config` object is
the single source of truth for a run, built either from an `.inp` file (`base/input.py:inp_para`)
or from API kwargs (`Config._load_from_kwargs` → `_set_defaults` → `_resolve_metadata` →
`_set_attributes`). `ConfigManager` caches parsed `.inp` files keyed by path so repeated `px.run()`
calls on the same file don't re-parse it (`force_reload=True` bypasses the cache — `run.py` always
does this). `_resolve_metadata` delegates to `api/_resolve.py:resolve_database_metadata`, which
auto-detects states-file columns/formats and database-specific physical constants (species IDs,
abundance, mass, capability flags like uncertainty/lifetime/g-factor/predissociation) by reading the
actual ExoMol `.def` file or the bundled `database/meta/hitran_molparam.txt`, so the Python API
doesn't need an `.inp` file's worth of explicit metadata.

**Execution (`core.py:get_results`)** is the orchestrator every entry point funnels through
(`run.py`, `cli.py`, and every `api/` function via `px.run`/individual calls). It: prints device
(CPU/GPU) and database info tables; lazily imports all `calculation/`, `save/`, `convert/`,
`process/`, `database/` submodules (kept lazy to avoid circular imports and to skip loading unused
database-format code); recomputes bin-size- and mass-dependent physics constants into
`base/constants.py` and `calculation/calcualte_line_profile.py` globals for the current run; then
branches on `config.database` to decide which files/columns are actually needed
(`NeedAllStates`/`NeedPartStates`/`NeedAllTransitions`) before touching any transition data, since
transition files can be very large. Every function toggle (`PartitionFunctions`, `StickSpectra`,
`CrossSections`, etc.) is an independent boolean read straight off `Config`.

**Data loading is database-pluggable.** `database/load_exomol.py`, `load_exomolhr.py`,
`load_hitran.py` each know how to read that database's native states/transitions format;
`database/data.py` (`loaddata`/`LoadedData`) adds an optional on-disk cache layer (Parquet, keyed by
wavenumber range and states-column set) so repeated calculations over the same line list skip
re-parsing raw ExoMol `.trans`/HITRAN `.par` files — `px.load(...)` returns a `LoadedData` that can
be reused across multiple `get_results` calls (see the `NeedAllTransitions`/`alltrans` check in
`core.py`, which raises if a calculation needs transitions the loaded data wasn't prepared for).
`base/large_file.py` handles chunked reading of oversized transition files
(`read_trans_chunks`/`is_large_trans_file`) so full line lists don't have to fit in memory at once.

**Calculation vs. save vs. convert are separate layers**, mirroring the database split:
- `calculation/` — pure numerics (line profiles, cross sections, partition functions, cooling
  functions, oscillator strengths, lifetimes, specific heats), CPU-only.
- `gpu/` — GPU-accelerated counterparts for cooling functions, stick spectra, cross sections, and
  combined stick-spectra+cross-section; `gpu/base_gpu.py` handles backend selection/fallback
  (`PyTorch-CUDA` → `CuPy-CUDA` → `MPS` → CPU) and array marshalling (`to_numpy`, `_backend_arrays`)
  so calculation code stays backend-agnostic. GPU mode never hard-fails: it auto-falls back to CPU
  formulas when CUDA/MPS runtimes or devices are unavailable. Note: MPS (macOS) only supports
  float32 and is documented as lower precision — CUDA or CPU is preferred for high-precision runs.
- `save/exomol/`, `save/exomolhr/`, `save/hitran/` — per-database, per-quantity output writers
  (mirror the `calculation/` quantity names). Adding support for a new output format/quantity
  combination usually means adding one file here, not touching `calculation/`.
- `convert/` — cross-database line-list format conversion (ExoMol⇄HITRAN, ExoMolHR→HITRAN),
  independent of the cross-section/spectra calculation path.
- `process/` — cross-cutting transforms used by multiple calculations: quantum-number filtering
  (`filter_qn.py`), HITRAN quantum-number/state handling (`hitran_qn.py`, `hitran_states.py`),
  multi-temperature partition function tables (`Q_multi_T.py`), and shared intensity math
  (`S_for_LTE_NLTE_Ab_Em.py`) used by both LTE and non-LTE (2-temperature) calculations.
- `download/` — fetches line-list/database files from ExoMol, ExoMolHR, ExoAtom, and HITRAN
  (`common.py` has shared HTTP helpers; one `download_<database>.py` per source).

**Concurrency.** `run.py` sets the multiprocessing start method to `'fork'` (POSIX only) and, on
macOS, monkeypatches `ProcessPoolExecutor` to `ThreadPoolExecutor` to avoid nested process-pool
deadlocks in long workflows — keep that platform guard in mind before assuming `ProcessPoolExecutor`
means real process parallelism. `ncpufiles`/`ncputrans` (transition-file-level vs. within-file
parallelism) and `chunk_size` control CPU parallelism/batching; `gpu_batch_lines`/`gpu_batch_grid`
control GPU memory usage the equivalent way.

**Logging (`base/log.py`).** `setup_logging(logpath, announce=verbose)` tees stdout/stderr to a log
file (`TeeStream`) when `LogFilePath`/`logs_path` is set; `LogFilePath None` (or API `log='none'`)
disables file logging independent of terminal `Verbose`/`verbose` output. `output_context(verbose)`
in `run.py`/`cli.py` scopes console verbosity for a single run.

**Result capture (`base/result.py`).** When `config.output` is `'memory'` or `'both'`,
`get_results` builds a `Result`/`Parameters` object (`active_result` module global in `core.py`) so
API callers can get results back as Python objects instead of only reading files from `save_path`.

## Working in this repo

- The package uses a `src/` layout (`src/pyexocross/`); `run.py` and every `test/*.py` script
  manually prepend `src/` to `sys.path` rather than relying on an editable install being present —
  preserve that pattern in any new top-level script.
- Heavy imports of calculation/save/convert submodules inside `get_results` are intentionally lazy
  (comment: "to avoid circular import issues") — don't hoist them to module level.
- `.inp` files are positional/keyword text config (see `input/templates/` and `input/examples/`);
  the first column of each row is a fixed keyword string that must not be renamed (per
  `docs/source/command.md`).
- Citation/versioning: package version lives in `pyproject.toml` and `src/pyexocross/__init__.py`
  (`__version__`) — keep them in sync when bumping.
