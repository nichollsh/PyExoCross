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


def find_exocross_xsec_file(molecule, temperature, pressure):
    """Find the ExoCross .xsec file matching the given T/P."""
    
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
    
    pattern = f"{molecule}*{temp_str}*{pressure_str}*.xsec"
    
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
    
    if len(matching_files) > 1:
        print(f"Warning: Found {len(matching_files)} matching files. Using first one.")
    
    return matching_files[0]


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
    
    return wl, xsec


def find_reference_hdf5_file(molecule):
    """Find the reference HDF5 file for the given molecule."""
    
    if not REF_DIR.exists():
        raise FileNotFoundError(f"Reference opacity directory not found: {REF_DIR}")
    
    # Look for HDF5 file with molecule name
    h5_files = list(REF_DIR.glob(f"*{molecule}*.h5"))
    
    if not h5_files:
        raise FileNotFoundError(
            f"No HDF5 reference file found for molecule {molecule} in {REF_DIR}"
        )
    
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
        
        # Sort by wavelength (ascending)
        sort_idx = np.argsort(wl_ref)
        wl_ref = wl_ref[sort_idx]
        opacity = opacity[sort_idx]
        
        return wl_ref, opacity, actual_temp, actual_pres


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


def make_comparison_plot(wl, xsec_exocross, xsec_ref, molecule, temperature, 
                        pressure, residual, rel_error, wl_range=None):
    """Create comparison plot with opacity and residuals."""
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # Apply wavelength range if specified
    if wl_range is not None:
        mask = (wl >= wl_range[0]) & (wl <= wl_range[1])
        wl_plot = wl[mask]
        xsec_exocross_plot = xsec_exocross[mask]
        xsec_ref_plot = xsec_ref[mask]
        residual_plot = residual[mask]
        rel_error_plot = rel_error[mask]
    else:
        wl_plot = wl
        xsec_exocross_plot = xsec_exocross
        xsec_ref_plot = xsec_ref
        residual_plot = residual
        rel_error_plot = rel_error
    
    # Plot 1: Opacities comparison
    ax1 = axes[0]
    ax1.loglog(wl_plot, xsec_exocross_plot, label='ExoCross', linewidth=2, alpha=0.8)
    ax1.loglog(wl_plot, xsec_ref_plot, label='Reference', linewidth=2, alpha=0.8)
    ax1.set_ylabel('Cross-section (cm²)', fontsize=11)
    ax1.legend(fontsize=10)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.set_title(f'{molecule}: T={temperature}K, P={pressure}bar', fontsize=12, fontweight='bold')
    
    # Plot 2: Absolute residual
    ax2 = axes[1]
    ax2.semilogy(wl_plot, np.abs(residual_plot), color='red', linewidth=2, alpha=0.8)
    ax2.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    ax2.set_ylabel('|Residual| (cm²)', fontsize=11)
    ax2.grid(True, which='both', alpha=0.3)
    
    # Plot 3: Relative error
    ax3 = axes[2]
    ax3.semilogy(wl_plot, np.abs(rel_error_plot), color='green', linewidth=2, alpha=0.8)
    ax3.set_ylabel('|Relative Error|', fontsize=11)
    ax3.set_xlabel('Wavelength (µm)', fontsize=11)
    ax3.grid(True, which='both', alpha=0.3)
    
    plt.tight_layout()
    return fig


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
    
    try:
        # Find and read ExoCross file
        print("Reading ExoCross cross-section...")
        xsec_file = find_exocross_xsec_file(args.molecule, args.temperature, args.pressure)
        print(f"  Found: {xsec_file}")
        wl_exo, xsec_exo = read_xsec_file(xsec_file)
        print(f"  Wavelength range: {wl_exo[0]:.4f} - {wl_exo[-1]:.4f} µm")
        print(f"  Data points: {len(wl_exo)}")
        
        # Find and read reference file
        print("\nReading reference opacity...")
        h5_file = find_reference_hdf5_file(args.molecule)
        print(f"  Found: {h5_file}")
        wl_ref, xsec_ref, actual_temp, actual_pres = read_reference_opacity(
            h5_file, args.temperature, args.pressure
        )
        print(f"  Wavelength range: {wl_ref[0]:.4f} - {wl_ref[-1]:.4f} µm")
        print(f"  Data points: {len(wl_ref)}")
        print(f"  Using T={actual_temp}K, P={actual_pres}bar")
        
        # Interpolate to common grid
        print("\nInterpolating to common wavelength grid...")
        wl_common, xsec_exo_interp, xsec_ref_interp = interpolate_to_common_grid(
            wl_exo, xsec_exo, wl_ref, xsec_ref
        )
        print(f"  Common grid points: {len(wl_common)}")
        
        # Calculate residuals
        residual, rel_error = calculate_residual(xsec_exo_interp, xsec_ref_interp)
        
        # Statistics
        valid_mask = np.isfinite(rel_error)
        mean_abs_error = np.mean(np.abs(residual[valid_mask]))
        mean_rel_error = np.mean(np.abs(rel_error[valid_mask]))
        
        print(f"\nComparison statistics:")
        print(f"  Mean absolute error: {mean_abs_error:.3e} cm²")
        print(f"  Mean relative error: {mean_rel_error:.3%}")
        print(f"  Max absolute error: {np.max(np.abs(residual[valid_mask])):.3e} cm²")
        print(f"  Max relative error: {np.max(np.abs(rel_error[valid_mask])):.3%}")
        
        # Create plot
        print("\nGenerating comparison plot...")
        fig = make_comparison_plot(
            wl_common, xsec_exo_interp, xsec_ref_interp,
            args.molecule, args.temperature, args.pressure,
            residual, rel_error,
            wl_range=wl_range
        )
        
        # Save plot
        output_file = output_dir / f"opacity_comparison_{args.molecule}_T{int(args.temperature)}K_P{args.pressure}bar.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"  Saved to: {output_file}")
        
        if args.show:
            plt.show()
        
        plt.close(fig)
        print("\nDone!")
        
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
