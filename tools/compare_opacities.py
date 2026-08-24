#!/usr/bin/env python3
"""
Compare opacities from ExoCross output (data/xsc/) with reference data (data/opacity_ruizhi/).

For a given molecule, temperature, and pressure, this script:
1. Reads the ExoCross .xsec output file
2. Reads the reference opacity from HDF5/NetCDF file
3. Interpolates to common wavelength grid if needed
4. Plots both opacities vs wavelength
5. Shows residual difference

Usage:
    python compare_opacities.py --molecule CO --temperature 500 --pressure 10.0
    python compare_opacities.py --molecule H2O --T 1000 --P 0.1
"""

import argparse
import os
import re
import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

# reference data 
REF_DIR = Path(__file__).parent.parent / "data" / "opacity_ruizhi"

# generated data
XSC_DIR = Path(__file__).parent.parent / "data" / "xsc" / "exocross" / "xsecs" / "files"

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Compare ExoCross opacities with reference data (m=species, T=kelvin, P=bar)."
    )
    parser.add_argument(
        "--molecule", "-m",
        type=str,
        required=True,
        help="Molecule name (e.g., CO, H2O, Na)"
    )
    parser.add_argument(
        "--temperature", "-T",
        type=float,
        required=True,
        help="Temperature in Kelvin"
    )
    parser.add_argument(
        "--pressure", "-P",
        type=float,
        required=True,
        help="Pressure in bar"
    )
    parser.add_argument(
        "--wavelength-range",
        type=str,
        default=None,
        help="Wavelength range for plot (e.g., '1,10' microns). Default: auto-detect"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory for plot. Default: current directory"
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display plot on screen (in addition to saving)"
    )
    return parser.parse_args()


def find_exocross_xsec_files(molecule, temperature, pressure):
    """
    Find all ExoCross .xsec files matching the given T/P.
    """
    if not XSC_DIR.exists():
        raise FileNotFoundError(f"ExoCross directory not found: {XSC_DIR}")

    # Build expected filename pattern
    # Format: MOLECULE__...__T{T}K__P{P}bar__...__.xsec
    temp_str = f"T{int(temperature)}K"

    # Format pressure - use scientific notation for very small/large values
    if pressure < 0.001 or pressure >= 1000:
        pressure_str = f"P{pressure:.2e}bar"
    else:
        pressure_str = f"P{pressure}bar"

    # Search recursively for matching files
    matching_files = list(XSC_DIR.glob("**/*.xsec"))
    matching_files = [f for f in matching_files if re.match(
        f".*{molecule}.*{temp_str}.*{pressure_str}.*",
        f.name
    )]

    if not matching_files:
        raise FileNotFoundError(
            f"No ExoCross file found for {molecule} at T={temperature}K, P={pressure}bar"
        )

    return sorted(matching_files)


def describe_xsec_source(path):
    """
    Build a short legend label (e.g. 'ExoMol (Li2015)') for an .xsec file.
    """
    database = path.parent.name
    tokens = path.name.split('__')
    dataset = None
    for i, tok in enumerate(tokens):
        if re.match(r'^(T\d|Tvib\d)', tok):
            data_info_tokens = tokens[:i]
            if len(data_info_tokens) >= 2:
                dataset = data_info_tokens[-1]
            break
    return f"{database} ({dataset})" if dataset else database


def read_xsec_file(filepath):
    """Read ExoCross .xsec file and return wavenumber and cross-section."""
    data = np.loadtxt(filepath)
    if data.size == 0:
        raise ValueError(f"Empty file: {filepath}")
    
    wn = data[:, 0]  # wavenumber in cm^-1
    xsec = data[:, 1]  # cross-section in cm^2
    
    # Convert wavenumber to wavelength in microns
    wl = 1e4 / wn  # 1 cm^-1 = 10^4 microns
    
    # Sort by wavelength (descending, since larger wn = shorter wl)
    sort_idx = np.argsort(wl)
    wl = wl[sort_idx]
    xsec = xsec[sort_idx]

    # Get bin widths in wavelength space
    wl_wid = np.abs(np.diff(wl))
    
    return wl[:-1], wl_wid, xsec[:-1]


def find_reference_hdf5_file(molecule):
    """Find the reference HDF5 file for the given molecule."""
    
    if not REF_DIR.exists():
        raise None
    
    # Look for HDF5 file with molecule name
    h5_files = list(REF_DIR.glob(f"*{molecule}*.h5"))
    
    if not h5_files:
        return None
    
    if len(h5_files) > 1:
        print(f"Warning: Found {len(h5_files)} HDF5 files. Using first one.")
    
    return h5_files[0]


def read_reference_opacity(filepath, temperature, pressure):
    """Read reference opacity from HDF5 file for given T/P."""
    with h5py.File(filepath, 'r') as f:
        temps = np.array(f['t'][:])  # Temperature grid
        pressures = np.array(f['p'][:])  # Pressure grid
        bin_edges = np.array(f['bin_edges'][:])  # Wavenumber bin edges
        xsec_arr = np.array(f['xsecarr'][:])  # [pressure, temperature, wavenumber]
        
        # Find closest temperature and pressure indices
        temp_idx = np.argmin(np.abs(temps - temperature))
        pres_idx = np.argmin(np.abs(pressures - pressure))
        
        actual_temp = temps[temp_idx]
        actual_pres = pressures[pres_idx]
        
        if abs(actual_temp - temperature) > 50:
            print(f"Warning: Using T={actual_temp}K (requested {temperature}K)")
        if abs(actual_pres - pressure) > 0.1 * pressure:
            print(f"Warning: Using P={actual_pres}bar (requested {pressure}bar)")
        
        # Extract opacity for this T/P
        # Note: dimensions are [pressure, temperature, wavelength]
        opacity = xsec_arr[pres_idx, temp_idx, :]

        # Convert bin_edges to wavelength in microns
        # bin_edges are in cm^-1, convert to microns
        wl_bin_edges = 1e4 / bin_edges  # wavelength in microns

        # Get bin centers
        wl_ref = 0.5 * (wl_bin_edges[:-1] + wl_bin_edges[1:])
        wl_wid = np.abs(np.diff(wl_bin_edges))  # bin widths in microns

        # Sort by wavelength (ascending)
        sort_idx = np.argsort(wl_ref)
        wl_wid = wl_wid[sort_idx]
        wl_ref = wl_ref[sort_idx]
        opacity = opacity[sort_idx]

        return wl_ref, wl_wid, opacity, actual_temp, actual_pres


def interpolate_to_common_grid(wl1, xsec1, wl2, xsec2):
    """Interpolate both datasets to a common wavelength grid."""
    # Use the denser grid as the target
    if len(wl1) >= len(wl2):
        wl_common = wl1
        # Interpolate reference data to ExoCross grid
        f = interp1d(wl2, xsec2, kind='linear', bounds_error=False, fill_value=np.nan)
        xsec2_interp = f(wl_common)
        xsec1_interp = xsec1
    else:
        wl_common = wl2
        # Interpolate ExoCross data to reference grid
        f = interp1d(wl1, xsec1, kind='linear', bounds_error=False, fill_value=np.nan)
        xsec1_interp = f(wl_common)
        xsec2_interp = xsec2
    
    # Remove NaN values
    valid_idx = ~(np.isnan(xsec1_interp) | np.isnan(xsec2_interp))
    wl_common = wl_common[valid_idx]
    xsec1_interp = xsec1_interp[valid_idx]
    xsec2_interp = xsec2_interp[valid_idx]
    
    return wl_common, xsec1_interp, xsec2_interp


def calculate_residual(xsec1, xsec2):
    """Calculate residual and relative error."""
    # Absolute residual
    residual = xsec1 - xsec2
    
    # Relative error (avoid division by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        relative_error = residual / xsec2
        relative_error = np.where(np.isfinite(relative_error), relative_error, 0)
    
    return residual, relative_error


def _masked(wl, values, wl_range):
    if wl_range is None:
        return wl, values
    mask = (wl >= wl_range[0]) & (wl <= wl_range[1])
    return wl[mask], values[mask]


def make_comparison_plot(sources, reference, molecule, temperature, pressure, wl_range=None):
    """Create a comparison plot overlaying every calculated opacity source.

    ``sources`` is a list of dicts, one per ``.xsec`` file found (e.g. one
    from ExoMol, one from HITRAN), each with keys ``label``, ``wl``,
    ``xsec_per_um`` and, when a reference is available, ``wl_common``,
    ``residual``, ``rel_error``. ``reference`` is ``None`` or a dict with
    ``wl``/``xsec_per_um`` for the reference HDF5 opacity. Residual/relative
    error panels are only drawn when a reference is available.
    """
    have_ref = reference is not None
    n_panels = 3 if have_ref else 1
    fig, axes = plt.subplots(n_panels, 1, figsize=(12, 10 if have_ref else 5), sharex=True)
    axes = np.atleast_1d(axes)

    # Plot 1: Opacities comparison (every calculated source, plus reference)
    ax1 = axes[0]
    for src in sources:
        wl_plot, xsec_plot = _masked(src['wl'], src['xsec_per_um'], wl_range)
        ax1.plot(wl_plot, xsec_plot, label=src['label'], linewidth=2, alpha=0.5)
    if have_ref:
        wl_plot, xsec_plot = _masked(reference['wl'], reference['xsec_per_um'], wl_range)
        ax1.plot(wl_plot, xsec_plot, label='Reference', linewidth=2,  color='black', alpha=0.5)
    ax1.set_ylabel('Cross-section (cm²/µm)', fontsize=11)
    ax1.legend(fontsize=10)
    ax1.set_yscale('log')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.set_title(f'{molecule}: T={temperature}K, P={pressure}bar', fontsize=12, fontweight='bold')

    if not have_ref:
        ax1.set_xlabel('Wavelength (µm)', fontsize=11)
        ax1.set_xscale('log')
        plt.tight_layout()
        return fig

    # Plot 2: Absolute residual (each source vs. reference)
    ax2 = axes[1]
    for src in sources:
        wl_plot, residual_plot = _masked(src['wl_common'], src['residual'], wl_range)
        ax2.plot(wl_plot, np.abs(residual_plot), label=src['label'], linewidth=2, alpha=0.5)
    ax2.axhline(y=0, color='black', alpha=0.5)
    ax2.set_ylabel('|Residual| (cm²)', fontsize=11)
    ax2.set_yscale('symlog', linthresh=1e-30)
    ax2.set_ylim(bottom=0)
    ax2.legend(fontsize=10)
    ax2.grid(True, which='both', alpha=0.3)

    # Plot 3: Relative error (each source vs. reference)
    ax3 = axes[2]
    for src in sources:
        wl_plot, rel_error_plot = _masked(src['wl_common'], src['rel_error'], wl_range)
        ax3.plot(wl_plot, np.abs(rel_error_plot), label=src['label'], linewidth=2, alpha=0.5)
    ax3.set_ylabel('|Relative Error|', fontsize=11)
    ax3.set_xlabel('Wavelength (µm)', fontsize=11)
    ax3.set_yscale('symlog', linthresh=1e-6)
    ax3.set_ylim(bottom=0)
    ax3.legend(fontsize=10)
    ax3.grid(True, which='both', alpha=0.3)
    ax3.set_xscale("log")

    plt.tight_layout()
    return fig

def calc_integrated(wl, wl_wid, xsec, wlmin, wlmax):
    """Calculate integrated cross-section over specified wavelength range."""
    mask = (wl >= wlmin) & (wl <= wlmax)
    if not np.any(mask):
        raise ValueError(f"No data points in the range {wlmin}-{wlmax} µm")
    
    xsec_selected = xsec[mask]
    wl_wid_selected = wl_wid[mask]
    
    # Use rectangular integration
    integrated_xsec = np.sum(xsec_selected * wl_wid_selected)
    
    return integrated_xsec

def main():
    args = parse_arguments()
    
    # Parse wavelength range if provided
    wl_range = None
    if args.wavelength_range:
        try:
            wl_min, wl_max = map(float, args.wavelength_range.split(','))
            wl_range = (wl_min, wl_max)
        except ValueError:
            print(f"Invalid wavelength range format: {args.wavelength_range}")
            sys.exit(1)
    
    # Output directory
    output_dir = Path(args.output_dir) if args.output_dir else Path.cwd()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Comparing opacities for {args.molecule} at T={args.temperature}K, P={args.pressure}bar")
    print()
    
    # Find and read every calculated (ExoCross) cross-section available
    print("Reading calculated cross-sections...")
    xsec_files = find_exocross_xsec_files(args.molecule, args.temperature, args.pressure)
    print(f"  Found {len(xsec_files)} source(s):")

    sources = []
    for xsec_file in xsec_files:
        label = describe_xsec_source(xsec_file)
        wl_exo, wl_wid_exo, xsec_exo = read_xsec_file(xsec_file)
        xsec_exo_per_um = xsec_exo / wl_wid_exo  # Convert to per micron for plotting
        print(f"    {label}")
        print(f"        File: {xsec_file}")
        print(f"        Wavelength range: {wl_exo[0]:.4f} - {wl_exo[-1]:.4f} µm, {len(wl_exo)} points")
        sources.append({
            'label': label,
            'path': xsec_file,
            'wl': wl_exo,
            'wl_wid': wl_wid_exo,
            'xsec': xsec_exo,
            'xsec_per_um': xsec_exo_per_um,
        })

    # Find and read reference file, if one exists for this molecule.
    # A reference isn't required to show/compare the calculated sources
    # with each other.
    print("\nLooking for reference opacity...")
    h5_file = find_reference_hdf5_file(args.molecule)
    if h5_file is None:
        print("  No reference opacity file found. Continuing without a reference; only calculated sources will be shown.")
        reference = None
    else:
        print(f"  Found: {h5_file}")
        wl_ref, wl_wid_ref, xsec_ref, actual_temp, actual_pres = read_reference_opacity(
            h5_file, args.temperature, args.pressure
        )
        xsec_ref_per_um = xsec_ref / wl_wid_ref  # Convert to per micron for plotting
        print(f"  Wavelength range: {wl_ref[0]:.4f} - {wl_ref[-1]:.4f} µm")
        print(f"  Data points: {len(wl_ref)}")
        print(f"  Using T={actual_temp}K, P={actual_pres}bar")
        reference = {
            'wl': wl_ref, 'wl_wid': wl_wid_ref, 'xsec_per_um': xsec_ref_per_um,
            'xsec': xsec_ref, 'actual_temp': actual_temp, 'actual_pres': actual_pres,
        }

        # When a reference is available, compare every calculated source
        # against it individually 
        print("\nComparing each calculated source against the reference...")
        for src in sources:
            wl_common, xsec_interp, xsec_ref_interp = interpolate_to_common_grid(
                src['wl'], src['xsec'],
                reference['wl'], reference['xsec']
            )

            wl_overlap_min = np.amin(wl_common) * 1.01
            wl_overlap_max = np.amax(wl_common) / 1.01
            integrated_src = calc_integrated(
                src['wl'], src['wl_wid'], src['xsec_per_um'], wl_overlap_min, wl_overlap_max
            )
            integrated_ref = calc_integrated(
                reference['wl'], reference['wl_wid'], reference['xsec_per_um'],
                wl_overlap_min, wl_overlap_max
            )

            residual, rel_error = calculate_residual(xsec_interp, xsec_ref_interp)
            valid_mask = np.isfinite(rel_error)
            mean_abs_error = np.mean(np.abs(residual[valid_mask]))
            mean_rel_error = np.mean(np.abs(rel_error[valid_mask]))

            print(f"\n  {src['label']} versus reference, overlap {wl_overlap_min:.4f}-{wl_overlap_max:.4f} µm:")
            print(f"      Integrated cross-section -- source: {integrated_src:.3e} cm², reference: {integrated_ref:.3e} cm²")
            print(f"      Relative difference (integrated): {((integrated_src - integrated_ref) / integrated_ref * 100):.3f}%")
            print(f"      Mean absolute error: {mean_abs_error:.3e} cm², mean relative error: {mean_rel_error:.3%}")
            print(f"      Max absolute error: {np.max(np.abs(residual[valid_mask])):.3e} cm²")
            print(f"      Max relative error: {np.max(np.abs(rel_error[valid_mask])):.3%} at wl={wl_common[np.argmax(np.abs(rel_error[valid_mask]))]:.4f} µm")

            src['wl_common'] = wl_common
            src['residual'] = residual
            src['rel_error'] = rel_error

    # Create plot
    print("\nGenerating comparison plot...")
    fig = make_comparison_plot(
        sources, reference,
        args.molecule, args.temperature, args.pressure,
        wl_range=wl_range
    )

    # Save plot
    fmt = 'pdf'
    output_file = output_dir / f"opacity_comparison_{args.molecule}_T{int(args.temperature)}K_P{args.pressure}bar.{fmt}"
    plt.savefig(output_file, dpi=250, bbox_inches='tight')
    print(f"  Saved to: {output_file}")

    if args.show:
        plt.show()

    plt.close(fig)
    print("\nDone!")


if __name__ == "__main__":
    main()
