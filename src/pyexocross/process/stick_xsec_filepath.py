import os

import numpy as np


def temperature_string_base(T, Tvib, Trot, NLTEMethod):
    """
    Base temperature string (no pressure/range), shared between stick spectra and cross sections.
    """
    if NLTEMethod == 'T' and Tvib is not None and Trot is not None:
        return f'Tvib{Tvib}K__Trot{Trot}K'
    elif NLTEMethod == 'D' and T is not None and Trot is not None:
        return f'T{T}K__Trot{Trot}K'
    else:
        # L or P (or fallback)
        T_val = T if T is not None else 0
        return f'T{T_val}K'


def temperature_pressure_string(T, P, temp_idx, NLTEMethod,
                                Tvib_list=None, Trot_list=None,
                                pressure_dependent=False):
    """
    Build temperature + pressure string for filenames.

    - Uses the same temperature form as stick spectra:
      * LTE / P:      T{T}K
      * Non-LTE T:    Tvib{Tvib}K__Trot{Trot}K
      * Non-LTE D:    T{T}K__Trot{Trot}K
    - Always uses the actual T / P for this file (no ranges).
    """
    # Temperature part
    if NLTEMethod == 'T' and temp_idx is not None and Tvib_list is not None and Trot_list is not None:
        Tvib_val = Tvib_list[temp_idx]
        Trot_val = Trot_list[temp_idx]
        T_str = temperature_string_base(T=None, Tvib=Tvib_val, Trot=Trot_val, NLTEMethod='T')
    else:
        T_str = temperature_string_base(T=T, Tvib=None, Trot=None, NLTEMethod=NLTEMethod)

    # Optional pressure part (single value, no range)
    if pressure_dependent and P is not None:
        # if P < 0.001 or P >= 1000:
        P_str = f'__P{P:.2e}bar'
        # else:
        #     P_str = f'__P{P}bar'
        T_str = T_str + P_str

    return T_str


def binsizestring(bin_size):
    """Format bin size consistently for output filenames."""
    return f'{float(bin_size):.4f}'


def crosssectiondetails(bin_size, unit_fn, cutoff, profile_label):
    """Build the shared bin-size, cutoff, and profile filename segment."""
    return (
        '__BinSize' + binsizestring(bin_size) + unit_fn
        + 'Cutoff' + str(cutoff) + '__' + profile_label.replace(' ', '')
    )


def stick_spectra_filepath(ss_folder, T, Tvib, Trot, str_min_wnl, str_max_wnl, unit_fn,
                           data_info, wn_wl, UncFilter, threshold, database, abs_emi, LTE_NLTE, photo,
                           NLTEMethod):
    """
    Build stick spectra output file path (shared naming for ExoMol and HITRAN).
    """
    temp_part = temperature_string_base(T, Tvib, Trot, NLTEMethod)
    prefix = '__'.join(data_info) + '__' + temp_part + '__'
    return (ss_folder + prefix + wn_wl.lower() + str_min_wnl + '-' + str_max_wnl + unit_fn
            + 'unc' + str(UncFilter) + '__thres' + str(threshold) + '__' + database + '__'
            + abs_emi + photo + LTE_NLTE + '.stick')


def cross_section_filepath(xsecs_folder, data_info,
                           T, P, temp_idx,
                           Tvib_list, Trot_list,
                           str_min_v, str_max_v, unit_fn, wn_wl,
                           UncFilter, threshold, database, abs_emi,
                           bin_size, cutoff, profile_label, LTE_NLTE, photo,
                           NLTEMethod, pressure_dependent):
    """
    Build cross-section (.xsec) output file path (shared naming for ExoMol and HITRAN).

    The temperature/pressure part follows the same rules as stick spectra and
    always uses the actual T / P for this file (no ranges).
    """
    temp_part = temperature_pressure_string(
        T=T,
        P=P,
        temp_idx=temp_idx,
        NLTEMethod=NLTEMethod,
        Tvib_list=Tvib_list,
        Trot_list=Trot_list,
        pressure_dependent=pressure_dependent,
    )
    prefix = '__'.join(data_info) + '__' + temp_part + '__'
    return (xsecs_folder + prefix + wn_wl.lower() + str_min_v + '-' + str_max_v + unit_fn
            + 'unc' + str(UncFilter) + '__thres' + str(threshold) 
            + crosssectiondetails(bin_size, unit_fn, cutoff, profile_label)
            + '__' + database + '__' + abs_emi + photo + LTE_NLTE + '.xsec'
    )


def xsecs_files_foldername(save_path, data_info, database):
    """
    Build the .xsec output folder path.

    Matches the construction in plot/plot_cross_section.py's save_xsec_file_plot(),
    so resume-related existence checks and the real writer can never drift apart.
    """
    return save_path + 'xsecs/files/' + data_info[0] + '/' + database + '/'


def cross_section_wn_range_strings(wn_grid, wn_wl, wn_wl_unit):
    """
    Predict the (str_min_v, str_max_v, unit_fn) triple that save_xsec_file_plot()
    will compute for this run's fixed wn_grid, for use by the resume pre-pass
    (called once, before any (T, P) point has been computed this run).

    Mirrors plot/plot_cross_section.py's save_xsec_file_plot() exactly:
    - wn_wl == 'WN': str_min_v/str_max_v come from floor(min(wn))/ceil(max(wn)),
      unit_fn = 'cm-1__'. Exact -- depends only on the raw wn_grid array, never on
      the computed cross-section values.
    - wn_wl containing 'L' (wavelength output): the real code converts wn to
      wavelength and masks to wn > 0 AND isfinite(wavelength) AND isfinite(xsec),
      via _axis_values_from_wavenumber(). This function reuses that same helper
      with an all-zero placeholder standing in for the not-yet-computed xsec array
      (zeros are finite), so the mask correctly reduces to wn > 0 -- matching the
      real run UNLESS the real cross-section values are non-finite somewhere
      inside wn > 0, which can't be predicted ahead of time. Worst case if that
      happens: a stale filename prediction causes a point to be recomputed rather
      than skipped -- not data loss.
    """
    wn_arr = np.asarray(wn_grid)
    if wn_wl == 'WN':
        min_v = np.min(wn_arr)
        max_v = np.max(wn_arr)
        return str(int(np.floor(min_v))), str(int(np.ceil(max_v))), 'cm-1__'
    elif 'L' in wn_wl:
        # Local import: plot_cross_section.py imports from this module at its own
        # top level, so a module-level import here would be circular.
        from pyexocross.plot.plot_cross_section import _axis_values_from_wavenumber
        placeholder_xsec = np.zeros_like(wn_arr, dtype=float)
        v_value, _, unit_fn, _ = _axis_values_from_wavenumber(wn_arr, placeholder_xsec, 'WL', wn_wl_unit)
        if len(v_value) == 0:
            return None, None, unit_fn
        return str(int(np.floor(np.min(v_value)))), str(int(np.ceil(np.max(v_value)))), unit_fn
    else:
        raise ValueError('Please choose wavenumber or wavelength and type in correct format: wn or wl.')


def existing_cross_section_file(xsecs_folder, data_info, T, P, temp_idx, Tvib_list, Trot_list,
                                str_min_v, str_max_v, unit_fn, wn_wl,
                                UncFilter, threshold, database, abs_emi,
                                bin_size, cutoff, profile_label, LTE_NLTE, photo,
                                NLTEMethod, pressure_dependent, compress_xsec_yn):
    """
    Return (actual_path, mtime) for the .xsec output this (T, P) point would
    produce this run, if it (or its bz2-compressed form) already exists on disk --
    else return None.

    Checks both the plain path and '<path>.bz2', matching _save_xsec_array's
    behaviour (plot/plot_cross_section.py), regardless of this run's
    compress_xsec_yn -- so a prior compressed-output run is still detected as
    "done" even if this run has compression toggled off, and vice versa.
    """
    base_path = cross_section_filepath(
        xsecs_folder, data_info, T, P, temp_idx, Tvib_list, Trot_list,
        str_min_v, str_max_v, unit_fn, wn_wl, UncFilter, threshold,
        database, abs_emi, bin_size, cutoff, profile_label, LTE_NLTE, photo,
        NLTEMethod, pressure_dependent,
    )
    for candidate in (base_path, base_path + '.bz2'):
        if os.path.isfile(candidate):
            return candidate, os.path.getmtime(candidate)
    return None


def compute_resume_skip_set(combos, xsecs_folder, data_info, Tvib_list, Trot_list,
                            str_min_v, str_max_v, unit_fn, wn_wl,
                            UncFilter, threshold, database, abs_emi,
                            bin_size, cutoff, profile_label, LTE_NLTE, photo,
                            NLTEMethod, pressure_dependent, compress_xsec_yn,
                            keep_most_recent=2):
    """
    Given `combos` -- an iterable of (key, T, P, temp_idx) describing every (T, P)
    grid point a run intends to compute -- return the set of `key` values to SKIP
    because their output .xsec (or .xsec.bz2) file already exists on disk.

    Among combos with an existing output file, the `keep_most_recent` (default 2)
    with the latest mtime are always EXCLUDED from the skip set (i.e. always
    redone), since a killed/timed-out prior run may have left them half-written.
    Ties in mtime are broken by `combos` order, later entries treated as more
    recent.

    IMPORTANT: must be called exactly once, before this run writes any new
    output file -- correctness of the mtime-based "most recent" detection
    depends on only ever comparing files that existed before this run started.
    This is what makes resume safe even when (T, P) points are computed
    concurrently and complete out of grid order (HITRAN, ExoMolHR).

    Returns
    -------
    set
        `key` values (as given in `combos`) to skip.
    """
    existing = []  # (key, mtime, order)
    for order, (key, T, P, temp_idx) in enumerate(combos):
        found = existing_cross_section_file(
            xsecs_folder, data_info, T, P, temp_idx, Tvib_list, Trot_list,
            str_min_v, str_max_v, unit_fn, wn_wl, UncFilter, threshold,
            database, abs_emi, bin_size, cutoff, profile_label, LTE_NLTE,
            photo, NLTEMethod, pressure_dependent, compress_xsec_yn,
        )
        if found is not None:
            _, mtime = found
            existing.append((key, mtime, order))

    if not existing:
        return set()

    existing_sorted = sorted(existing, key=lambda item: (item[1], item[2]), reverse=True)
    return {key for key, _, _ in existing_sorted[keep_most_recent:]}
