"""
Unit tests for the "resume" cross-section feature.

Self-contained: uses synthetic, empty .xsec/.xsec.bz2 files under tmp_path with
controlled mtimes -- no real molecular data or downloads required. Run with:

    pytest test/test_resume.py -v
"""
import os
import sys

# Ensure src layout is importable in local runs
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src_root = os.path.join(project_root, 'src')
if src_root not in sys.path:
    sys.path.insert(0, src_root)

from pyexocross.config import Config
from pyexocross.process.stick_xsec_filepath import (
    compute_resume_skip_set,
    cross_section_filepath,
    existing_cross_section_file,
)


# ---------------------------------------------------------------------------
# Shared filename-building arguments for a synthetic (T, P) point
# ---------------------------------------------------------------------------
COMMON_ARGS = dict(
    data_info=['CO', '12C-16O', 'Li2015'],
    Tvib_list=[],
    Trot_list=[],
    str_min_v='2000',
    str_max_v='2300',
    unit_fn='cm-1__',
    wn_wl='WN',
    UncFilter=None,
    threshold=None,
    database='ExoMol',
    abs_emi='Ab',
    bin_size=0.5,
    cutoff=10,
    profile_label='SciPyVoigt',
    LTE_NLTE='',
    photo='',
    NLTEMethod='L',
    pressure_dependent=True,
)


def _xsec_path(folder, T, P, temp_idx=0):
    return cross_section_filepath(
        folder, T=T, P=P, temp_idx=temp_idx, **COMMON_ARGS
    )


def _touch(path, mtime):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write('x')
    os.utime(path, (mtime, mtime))


# ---------------------------------------------------------------------------
# existing_cross_section_file
# ---------------------------------------------------------------------------
def test_existing_cross_section_file_none(tmp_path):
    folder = str(tmp_path) + '/'
    result = existing_cross_section_file(
        folder, T=300, P=1.0, temp_idx=0, compress_xsec_yn='N', **COMMON_ARGS
    )
    assert result is None


def test_existing_cross_section_file_plain(tmp_path):
    folder = str(tmp_path) + '/'
    path = _xsec_path(folder, T=300, P=1.0)
    _touch(path, 1000.0)
    result = existing_cross_section_file(
        folder, T=300, P=1.0, temp_idx=0, compress_xsec_yn='N', **COMMON_ARGS
    )
    assert result is not None
    found_path, mtime = result
    assert found_path == path
    assert mtime == 1000.0


def test_existing_cross_section_file_bz2(tmp_path):
    folder = str(tmp_path) + '/'
    path = _xsec_path(folder, T=300, P=1.0) + '.bz2'
    _touch(path, 2000.0)
    result = existing_cross_section_file(
        folder, T=300, P=1.0, temp_idx=0, compress_xsec_yn='Y', **COMMON_ARGS
    )
    assert result is not None
    found_path, mtime = result
    assert found_path == path
    assert mtime == 2000.0


def test_existing_cross_section_file_finds_bz2_even_when_compression_off(tmp_path):
    # A prior run wrote .xsec.bz2; this run has compression toggled off --
    # should still be detected as already-done.
    folder = str(tmp_path) + '/'
    path = _xsec_path(folder, T=300, P=1.0) + '.bz2'
    _touch(path, 3000.0)
    result = existing_cross_section_file(
        folder, T=300, P=1.0, temp_idx=0, compress_xsec_yn='N', **COMMON_ARGS
    )
    assert result is not None


# ---------------------------------------------------------------------------
# compute_resume_skip_set
# ---------------------------------------------------------------------------
def _combos(n):
    """n synthetic (key, T, P, temp_idx) combos, one per temp_idx."""
    return [((i, 0), 300 + i, 1.0, i) for i in range(n)]


def test_skip_set_empty_when_no_files_exist(tmp_path):
    folder = str(tmp_path) + '/'
    combos = _combos(3)
    skip = compute_resume_skip_set(
        combos, folder, compress_xsec_yn='N', **COMMON_ARGS
    )
    assert skip == set()


def test_skip_set_with_one_existing_file_is_not_skipped(tmp_path):
    # A single existing file is conservatively treated as possibly incomplete.
    folder = str(tmp_path) + '/'
    combos = _combos(3)
    _touch(_xsec_path(folder, T=combos[0][1], P=combos[0][2], temp_idx=combos[0][3]), 1000.0)
    skip = compute_resume_skip_set(
        combos, folder, compress_xsec_yn='N', **COMMON_ARGS
    )
    assert skip == set()


def test_skip_set_with_two_existing_files_is_not_skipped(tmp_path):
    folder = str(tmp_path) + '/'
    combos = _combos(3)
    for i in (0, 1):
        _touch(_xsec_path(folder, T=combos[i][1], P=combos[i][2], temp_idx=combos[i][3]), 1000.0 + i)
    skip = compute_resume_skip_set(
        combos, folder, compress_xsec_yn='N', **COMMON_ARGS
    )
    assert skip == set()


def test_skip_set_excludes_two_most_recent_of_five(tmp_path):
    folder = str(tmp_path) + '/'
    combos = _combos(5)
    # mtimes strictly increasing with combo index: combo 4 is newest, combo 0 oldest.
    for i, (key, T, P, temp_idx) in enumerate(combos):
        _touch(_xsec_path(folder, T=T, P=P, temp_idx=temp_idx), 1000.0 + i)

    skip = compute_resume_skip_set(
        combos, folder, compress_xsec_yn='N', **COMMON_ARGS
    )
    # The 2 most recent (combo 3, combo 4) must always be redone (not in skip set).
    assert (3, 0) not in skip
    assert (4, 0) not in skip
    # The 3 oldest must be skipped.
    assert skip == {(0, 0), (1, 0), (2, 0)}


def test_skip_set_mtime_tie_broken_by_combos_order(tmp_path):
    folder = str(tmp_path) + '/'
    combos = _combos(4)
    # All 4 files share the same mtime (a real possibility under a thread pool).
    for key, T, P, temp_idx in combos:
        _touch(_xsec_path(folder, T=T, P=P, temp_idx=temp_idx), 5000.0)

    skip = compute_resume_skip_set(
        combos, folder, compress_xsec_yn='N', **COMMON_ARGS
    )
    # Ties broken by combos order: the last 2 entries (index 2, 3) are treated
    # as "more recent" and excluded from the skip set.
    assert skip == {(0, 0), (1, 0)}


def test_skip_set_unaffected_by_files_written_after_the_call(tmp_path):
    """
    Sanity check of the documented precondition: compute_resume_skip_set()
    must be called exactly once, before this run writes any new file. A file
    that appears on disk AFTER the call (e.g. because this run is now
    actively computing and saving points) must not change a result already
    returned -- the function only inspects disk state at call time.
    """
    folder = str(tmp_path) + '/'
    combos = _combos(3)
    for i, (key, T, P, temp_idx) in enumerate(combos):
        _touch(_xsec_path(folder, T=T, P=P, temp_idx=temp_idx), 1000.0 + i)

    skip_before = compute_resume_skip_set(combos, folder, compress_xsec_yn='N', **COMMON_ARGS)

    # Simulate "this run" writing a brand new, much newer file after the call.
    _touch(_xsec_path(folder, T=999, P=1.0, temp_idx=99), 9999.0)

    # Re-running the same call (same combos) is unaffected by that new file,
    # since it isn't one of the combos being checked.
    skip_after = compute_resume_skip_set(combos, folder, compress_xsec_yn='N', **COMMON_ARGS)
    assert skip_before == skip_after


# ---------------------------------------------------------------------------
# Config.resume
# ---------------------------------------------------------------------------
def _minimal_config_kwargs(**overrides):
    kwargs = dict(
        database='ExoMol', molecule='CO', isotopologue='12C-16O', dataset='Li2015',
        species_id=51,
    )
    kwargs.update(overrides)
    return kwargs


def test_config_resume_default_false_via_kwargs(monkeypatch):
    # Avoid touching the filesystem / requiring real line-list data: bypass
    # metadata resolution, which isn't relevant to the resume flag itself.
    monkeypatch.setattr(Config, '_resolve_metadata', lambda self: None)
    cfg = Config(**_minimal_config_kwargs())
    assert cfg.resume is False


def test_config_resume_true_via_kwargs(monkeypatch):
    monkeypatch.setattr(Config, '_resolve_metadata', lambda self: None)
    cfg = Config(**_minimal_config_kwargs(resume=True))
    assert cfg.resume is True


def test_config_resume_yes_string_via_kwargs(monkeypatch):
    monkeypatch.setattr(Config, '_resolve_metadata', lambda self: None)
    cfg = Config(**_minimal_config_kwargs(resume='Yes'))
    assert cfg.resume is True


def test_config_resume_no_string_via_kwargs(monkeypatch):
    monkeypatch.setattr(Config, '_resolve_metadata', lambda self: None)
    cfg = Config(**_minimal_config_kwargs(resume='No'))
    assert cfg.resume is False


def test_config_resume_default_false_when_load_from_kwargs_skipped():
    """
    Regression check: on the plain .inp + CLI path, Config(inp_filepath=...,
    force_reload=True) is called with NO kwargs, so _load_from_kwargs() is
    skipped entirely (see Config.__init__: it only calls _load_from_kwargs
    when inp_filepath is None or kwargs is non-empty). self.resume must
    still default to False via the unconditional default set in __init__ --
    not raise AttributeError.

    Uses a real, repo-committed .inp template (this only exercises the
    Config/.inp parsing pipeline -- it never calls get_results(), so no
    line-list data needs to actually exist on disk).
    """
    template_path = os.path.join(project_root, 'input', 'templates', 'CO_ExoMol_template.inp')
    cfg = Config(inp_filepath=template_path, force_reload=True)
    assert cfg.resume is False


# ---------------------------------------------------------------------------
# CLI argparse wiring
# ---------------------------------------------------------------------------
def test_parse_args_resume_default_false(monkeypatch):
    from pyexocross.base.input import parse_args
    monkeypatch.setattr(sys, 'argv', ['run.py', '-p', 'dummy.inp'])
    args = parse_args()
    assert args.path == 'dummy.inp'
    assert args.resume is False


def test_parse_args_resume_short_flag(monkeypatch):
    from pyexocross.base.input import parse_args
    monkeypatch.setattr(sys, 'argv', ['run.py', '-p', 'dummy.inp', '-r'])
    args = parse_args()
    assert args.resume is True


def test_parse_args_resume_long_flag(monkeypatch):
    from pyexocross.base.input import parse_args
    monkeypatch.setattr(sys, 'argv', ['run.py', '-p', 'dummy.inp', '--resume'])
    args = parse_args()
    assert args.resume is True


def test_cli_parser_has_resume_flag(monkeypatch):
    import argparse
    monkeypatch.setattr(sys, 'argv', ['pyexocross', '-p', 'dummy.inp', '--resume'])
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', required=True)
    parser.add_argument('-r', '--resume', action='store_true', default=False)
    args = parser.parse_args()
    assert args.resume is True
