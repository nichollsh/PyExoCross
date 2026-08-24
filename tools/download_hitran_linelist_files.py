# Import all what we need
import os
import sys

# Ensure the src-layout package is importable
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src_root = os.path.join(project_root, 'src')
if src_root not in sys.path:
    sys.path.insert(0, src_root)
import pyexocross as px

# File Paths and Molecules
file_path = '/scratch/p321409/opacity_lbl/hitran/'

# Wavenumber range [cm-1] to download from HITRAN.
wn_range = None

# Molecules, isotopologues (ExoMol-style slug, e.g. '1H2-16O') and wavenumber
# ranges (cm-1) to download from HITRAN. wn_range is required for every
# isotopologue, since HITRAN's line-by-line API only serves a bounded range.
molecule_isotopologues = {
    'H2O': {'1H2-16O':  {'wn_range': wn_range}},
    'CO':  {'12C-16O':  {'wn_range': wn_range}},
    'O2':  {'16O2':     {'wn_range': wn_range}},
    'CO2': {'12C-16O2': {'wn_range': wn_range}},
    'NO':  {'14N-16O':  {'wn_range': wn_range}},
    'O3':  {'16O3':     {'wn_range': wn_range}},
}

# Name a specific molecule to download
molec_tgt = 'O2'
download_tgt = {molec_tgt: molecule_isotopologues[molec_tgt]}

# Download line list files.
if __name__ == '__main__':
    px.download_hitran(
        species_info=download_tgt,
        file_path=file_path,
        download=True,
    )
