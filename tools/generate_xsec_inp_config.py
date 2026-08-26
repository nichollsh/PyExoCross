#!/usr/bin/env python
 
"""
Generate ExoMol/ExoAtom/HITRAN cross-section input files.

For ExoMol/ExoAtom, quantum-number labels/formats and dataset metadata come
from the downloaded ``.def.json``/``.adef.json`` files. HITRAN has no such
definition file, so its input files are generated instead from the
``.par`` files downloaded by
``ExoMol_extra_tools/download_hitran_linelist_files.py`` (for which
molecule/isotopologue directories are discovered under ``--data-dir``) plus
the HITRAN molecule/isotopologue classification tables built into
``pyexocross.process.hitran_qn`` and the HITRAN isotopologue metadata page
(``hitran.org/docs/iso-meta/``, via ``pyexocross.download.download_hitran``)
used to resolve each isotopologue's ``SpeciesID``. Resolving HITRAN metadata
therefore requires ``pyexocross`` to be importable and requires internet
access.

Usage:
  python ExoMol_extra_tools/generate_exomolatom_xsec_inp.py --database ExoMol --molecule CO
  python ExoMol_extra_tools/generate_exomolatom_xsec_inp.py --database ExoAtom --molecule Na
  python ExoMol_extra_tools/generate_exomolatom_xsec_inp.py --database HITRAN --molecule H2O \\
      --data-dir /scratch/p321409/opacity_lbl/hitran/
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Iterable, List, Tuple


EXCLUDE_FIELD_NAMES = {
    "id",
    "e",
    "gtot",
    "j",
    "unc",
    "tau",
    "gfactor",
    "gfac",
    "g-factor",
}

LABEL_PREFIXES = (
    "Herzberg:",
    "Polyad:",
    "Hundb:",
    "hundb:",
    "hunda:",
    "TROVE:",
    "AFGL:",
)

UNC_DEFAULT = 0.1  # default uncertainty value to write in xsec input files

def _sanitize_label(label: str) -> str:
    return re.sub(r"\s+", "", label.strip())


def _format_label(raw_label: str, used: set[str]) -> str:
    raw_label = _sanitize_label(raw_label)
    base_label = raw_label
    for prefix in LABEL_PREFIXES:
        if base_label.lower().startswith(prefix.lower()):
            base_label = base_label[len(prefix) :]
            break
    base_label = _sanitize_label(base_label)
    if base_label == "Grve":
        base_label = "Gamma_rve"
    label = base_label if base_label and base_label not in used else raw_label
    if label in used:
        suffix = 2
        while f"{label}_{suffix}" in used:
            suffix += 1
        label = f"{label}_{suffix}"
    return label


def _should_skip_field(name: str) -> bool:
    lower_name = name.lower()
    if lower_name in EXCLUDE_FIELD_NAMES:
        return True
    if lower_name.startswith("auxiliary:"):
        return True
    if "coef" in lower_name:
        return True
    return False


def _extract_qns(fields: Iterable[dict]) -> Tuple[List[str], List[str]]:
    labels: List[str] = []
    formats: List[str] = []
    used: set[str] = set()
    for field in fields:
        name = field.get("name", "")
        if not name or _should_skip_field(name):
            continue
        cfmt = str(field.get("cfmt", "")).strip()
        if not cfmt:
            continue
        label = _format_label(name, used)
        labels.append(label)
        formats.append(cfmt)
        used.add(label)
    return labels, formats


def _read_template(template_path: Path) -> List[str]:
    return template_path.read_text().splitlines(keepends=True)


def _replace_line(lines: List[str], key: str, value: str) -> bool:
    for idx, line in enumerate(lines):
        if line.lstrip().startswith(key):
            lines[idx] = f"{key:<40}{value}\n"
            return True
    return False


def _update_unc_filter(lines: List[str], has_uncertainty: bool, unc_default: float) -> None:
    value = "Y" if has_uncertainty else "N"
    for idx, line in enumerate(lines):
        if line.lstrip().startswith("UncFilter(Y/N)"):
            lines[idx] = f"{'UncFilter(Y/N)':<40}{value}          {unc_default}\n"
            return


def _find_line_value(lines: List[str], key: str) -> str:
    for line in lines:
        if line.lstrip().startswith(key):
            parts = line.split()
            return parts[1] if len(parts) > 1 else ""
    return ""


def _update_log_file_path(lines: List[str], species: str, dataset: str, database: str) -> bool:
    # Point the log at the same directory PyExoCross writes xsec files to
    # (save_path+'xsecs/files/<species>/<database>/', see plot_cross_section.py),
    # so the log ends up alongside the .xsec output instead of in SavePath's root.
    save_path = _find_line_value(lines, "SavePath")
    log_dir = os.path.join(save_path, "xsecs", "files", species, database) if save_path else ""
    log_name = f"{species}_{dataset}.log"
    new_path = os.path.join(log_dir, log_name) if log_dir else log_name
    for idx, line in enumerate(lines):
        if line.lstrip().startswith("LogFilePath"):
            lines[idx] = f"{'LogFilePath':<40}{new_path}\n"
            return True
    return False


def _find_def_files(data_dir: Path, species: str, database: str) -> List[Path]:
    species_lower = species.lower()
    matches: List[Path] = []
    db_marker = "exoatom" if database == "ExoAtom" else "exomol"
    suffix = ".adef.json" if database == "ExoAtom" else ".def.json"
    for root, _dirs, files in os.walk(data_dir, followlinks=True):
        if not files:
            continue
        root_parts = Path(root).parts

        if db_marker not in root_parts:
            continue

        try:
            mol_part = root_parts[root_parts.index(db_marker) + 1]
        except (ValueError, IndexError):
            continue
        if mol_part.lower() != species_lower:
            continue
        print("Found directory for species '{}'\n\t{}".format(species, root))
        for fname in files:
            if fname.endswith(suffix):
                matches.append(Path(root) / fname)
    return sorted(matches)


def _ensure_pyexocross_importable() -> None:
    """Make the src-layout ``pyexocross`` package importable, if not already."""
    repo_root = Path(__file__).resolve().parents[1]
    src_root = str(repo_root / "src")
    if src_root not in sys.path:
        sys.path.insert(0, src_root)


def _find_hitran_species(data_dir: Path, species: str) -> List[Tuple[str, Path]]:
    """Find downloaded HITRAN isotopologues for ``species`` under ``data_dir``.

    Mirrors the directory layout written by
    ``ExoMol_extra_tools/download_hitran_linelist_files.py``:
    ``<data_dir>/<species>/<isotopologue>/<species>__<isotopologue>.par``.
    """
    species_dir = data_dir / species
    matches: List[Tuple[str, Path]] = []
    if not species_dir.is_dir():
        return matches
    for iso_dir in sorted(species_dir.iterdir()):
        if not iso_dir.is_dir():
            continue
        par_path = iso_dir / f"{species}__{iso_dir.name}.par"
        if par_path.is_file():
            matches.append((iso_dir.name, par_path))
    return matches


def _drop_labels(labels: List[str], formats: List[str], targets: set) -> Tuple[List[str], List[str]]:
    kept_labels: List[str] = []
    kept_formats: List[str] = []
    for label, fmt in zip(labels, formats):
        if label in targets:
            continue
        kept_labels.append(label)
        kept_formats.append(fmt)
    return kept_labels, kept_formats


def _hitran_qn_labels_formats(species: str, isotopologue: str) -> Tuple[List[str], List[str]]:
    """Derive QNslabel/QNsformat for a HITRAN molecule/isotopologue.

    Reproduces -- in the same order -- the column-derivation rules
    ``pyexocross.process.hitran_qn.separate_QN_hitran`` applies at runtime to
    build the *upper*-state quantum-number columns (``QNu_col``) that
    ``QNsFilter`` entries are matched against: drop ``'none'`` placeholders
    from the global and local-upper/local-lower label sets; if the local
    group carries a literal ``'Br'`` branch label on the *lower* side,
    synthesize an upper-state ``J`` from it (mirroring
    ``LQNu_df['J'] = LQNl_df['Br'] + LQNl_df['J']``) before dropping ``'Br'``;
    drop ``'m'``/``'M'`` branch labels; and, if the molecule's raw upper
    label set has no ``'F'``, append ``F = J`` on the upper side (mirroring
    ``LQNu_df['F'] = LQNu_df['J']``). Only the resulting upper-state set is
    returned, matching the known-correct H2O example
    (``input/examples/H2O_HITRAN_xsec.inp``), where local upper and lower
    labels coincide.

    Also warns when ``species``/``isotopologue`` isn't recognized by any of
    the classification tables in ``hitran_qn.py`` (i.e. it fell back to the
    generic default labels), since the derived QNslabel/QNsformat is then a
    guess, not a lookup, and should be hand-verified against the downloaded
    .par file's Q'/Q" columns.
    """
    _ensure_pyexocross_importable()
    import pyexocross.base.input as _input
    from pyexocross.process.hitran_qn import globalQNclasses, localQNgroups

    data_info = [species, isotopologue]
    global_labels_raw, global_formats_raw = globalQNclasses(data_info)
    upper_labels_raw, lower_labels_raw, upper_formats_raw, lower_formats_raw = localQNgroups(data_info)
    global_labels: List[str] = list(global_labels_raw or [])
    global_formats: List[str] = list(global_formats_raw or [])
    orig_upper_labels: List[str] = list(upper_labels_raw or [])
    upper_labels: List[str] = list(orig_upper_labels)
    upper_formats: List[str] = list(upper_formats_raw or [])
    lower_labels: List[str] = list(lower_labels_raw or [])
    lower_formats: List[str] = list(lower_formats_raw or [])

    if (
        global_labels == list(_input.GlobalQNLabel_list)
        and orig_upper_labels == list(_input.LocalQNLabel_list)
    ):
        print(
            f"Warning: {species}/{isotopologue} is not recognized by any "
            "pyexocross.process.hitran_qn classification table; the derived "
            "QNslabel/QNsformat fell back to generic defaults and should be "
            "hand-verified against the downloaded .par file's Q'/Q\" columns."
        )

    global_labels, global_formats = _drop_labels(global_labels, global_formats, {"none"})
    upper_labels, upper_formats = _drop_labels(upper_labels, upper_formats, {"none"})
    lower_labels, lower_formats = _drop_labels(lower_labels, lower_formats, {"none"})

    if "Br" in upper_labels:  # not hit by any current table; kept for parity with separate_QN_hitran
        upper_labels, upper_formats = _drop_labels(upper_labels, upper_formats, {"Br"})
    if "Br" in lower_labels:
        if "J" not in upper_labels:
            j_format = lower_formats[lower_labels.index("J")] if "J" in lower_labels else "%3d"
            upper_labels.append("J")
            upper_formats.append(j_format)
        lower_labels, lower_formats = _drop_labels(lower_labels, lower_formats, {"Br"})

    upper_labels, upper_formats = _drop_labels(upper_labels, upper_formats, {"m", "M"})
    lower_labels, lower_formats = _drop_labels(lower_labels, lower_formats, {"m", "M"})

    if "F" not in orig_upper_labels:
        j_format = upper_formats[upper_labels.index("J")] if "J" in upper_labels else "%3d"
        upper_labels.append("F")
        upper_formats.append(j_format)

    labels = list(global_labels) + upper_labels
    formats = list(global_formats) + upper_formats
    return labels, formats


def _hitran_species_id(species: str, isotopologue: str) -> int:
    """Resolve a HITRAN SpeciesID (molecule_id*10 + local isotopologue id).

    Looks up ``molecule_id``/``local_id`` from HITRAN's iso-meta page via
    ``pyexocross.download.download_hitran``, the same resolution the
    download step already relies on -- requires internet access.
    """
    _ensure_pyexocross_importable()
    from pyexocross.download.download_hitran import _resolve_iso_meta

    meta = _resolve_iso_meta(species, isotopologue)
    return meta["molecule_id"] * 10 + meta["local_id"]


def generate_hitran_inputs(
    species: str,
    data_dir: Path,
    template_path: Path,
    output_dir: Path,
) -> List[Path]:
    entries = _find_hitran_species(data_dir, species)
    if not entries:
        raise FileNotFoundError(
            f"No downloaded HITRAN .par files found for '{species}' under {data_dir} "
            f"(expected {data_dir}/{species}/<isotopologue>/{species}__<isotopologue>.par)."
        )

    multi_iso = len(entries) > 1
    output_dir.mkdir(parents=True, exist_ok=True)
    template_lines = _read_template(template_path)
    created: List[Path] = []

    for isotopologue, _par_path in entries:
        labels, formats = _hitran_qn_labels_formats(species, isotopologue)
        species_id = _hitran_species_id(species, isotopologue)
        dataset = f"{species}-HITRAN"

        out_lines = list(template_lines)
        if not _replace_line(out_lines, "Database", "HITRAN"):
            raise ValueError("Template must include a 'Database' line.")
        if not _replace_line(out_lines, "Molecule", species):
            raise ValueError("Template must include a 'Molecule' line.")
        if not _replace_line(out_lines, "Isotopologue", isotopologue):
            raise ValueError("Template must include an 'Isotopologue' line.")
        if not _replace_line(out_lines, "Dataset", dataset):
            raise ValueError("Template must include a 'Dataset' line.")
        if not _replace_line(out_lines, "SpeciesID", str(species_id)):
            raise ValueError("Template must include a 'SpeciesID' line.")
        if not _update_log_file_path(out_lines, species, dataset, "HITRAN"):
            raise ValueError("Template must include a 'LogFilePath' line.")
        if not _replace_line(out_lines, "QNslabel", "  ".join(labels)):
            raise ValueError("Template must include a 'QNslabel' line.")
        if not _replace_line(out_lines, "QNsformat", "  ".join(formats)):
            raise ValueError("Template must include a 'QNsformat' line.")
        if not _replace_line(
            out_lines,
            "QNsFilter(Y/N)",
            f"N          {'  '.join(f'{label}[]' for label in labels)}",
        ):
            raise ValueError("Template must include a 'QNsFilter(Y/N)' line.")
        # HITRAN linelists always carry per-line uncertainty codes.
        _update_unc_filter(out_lines, has_uncertainty=True, unc_default=UNC_DEFAULT)

        filename = (
            f"{species}_{isotopologue}_HITRAN_xsec.inp" if multi_iso else f"{species}_HITRAN_xsec.inp"
        )
        out_path = output_dir / filename
        out_path.write_text("".join(out_lines))
        created.append(out_path)

    return created


def _output_name(
    database: str,
    species: str,
    isotopologue: str | None,
    dataset: str,
    multi_files: bool,
    multi_iso: bool,
    multi_dataset: bool,
) -> str:

    if database == "ExoAtom":
        return f"{species}_ExoAtom_{dataset}_xsec.inp"
    
    elif database == "ExoMol":
        return f"{species}_ExoMol_{dataset}_xsec.inp"

    else:
        raise Exception(f"Invalid database {database}")


def generate_inputs(
    database: str,
    species: str,
    data_dir: Path,
    template_path: Path,
    output_dir: Path,
    species_id: int,
) -> List[Path]:
    def_files = _find_def_files(data_dir, species, database)
    if not def_files:
        raise FileNotFoundError(f"No definition files found for '{species}' in {database}.")

    meta = []
    for def_path in def_files:
        with def_path.open() as handle:
            payload = json.load(handle)
        if database == "ExoAtom":
            dataset = payload["dataset"]["name"]
            fields = payload["dataset"]["states"]["states_file_fields"]
            has_uncertainty = bool(
                payload["dataset"]["states"].get(
                    "uncertainties_available",
                    payload["dataset"]["states"].get("uncertainty_available", False),
                )
            )
            meta.append(
                {
                    "path": def_path,
                    "isotopologue": None,
                    "dataset": dataset,
                    "fields": fields,
                    "has_uncertainty": has_uncertainty,
                }
            )
        else:
            meta.append(
                {
                    "path": def_path,
                    "isotopologue": payload["isotopologue"]["iso_slug"],
                    "dataset": payload["dataset"]["name"],
                    "fields": payload["dataset"]["states"]["states_file_fields"],
                    "has_uncertainty": bool(payload["dataset"]["states"].get("uncertainties_available", False)),
                }
            )

    isotopologues = {m["isotopologue"] for m in meta if m["isotopologue"]}
    datasets = {m["dataset"] for m in meta}
    multi_files = len(meta) > 1
    multi_iso = len(isotopologues) > 1
    multi_dataset = len(datasets) > 1

    output_dir.mkdir(parents=True, exist_ok=True)
    template_lines = _read_template(template_path)
    created: List[Path] = []

    for entry in meta:
        labels, formats = _extract_qns(entry["fields"])
        out_lines = list(template_lines)
        if database == "ExoAtom":
            if not _replace_line(out_lines, "Database", "ExoAtom"):
                raise ValueError("Template must include a 'Database' line.")
            if not _replace_line(out_lines, "Atom", species):
                raise ValueError("Template must include an 'Atom' line for ExoAtom.")
            _replace_line(out_lines, "Molecule", species)
            _replace_line(out_lines, "Isotopologue", "none")
        else:
            if not _replace_line(out_lines, "Database", "ExoMol"):
                raise ValueError("Template must include a 'Database' line.")
            if not _replace_line(out_lines, "Molecule", species):
                raise ValueError("Template must include a 'Molecule' line.")
            if entry["isotopologue"] and not _replace_line(out_lines, "Isotopologue", entry["isotopologue"]):
                raise ValueError("Template must include an 'Isotopologue' line.")

        if not _replace_line(out_lines, "Dataset", entry["dataset"]):
            raise ValueError("Template must include a 'Dataset' line.")
        if not _replace_line(out_lines, "SpeciesID", str(species_id)):
            raise ValueError("Template must include a 'SpeciesID' line.")
        if not _update_log_file_path(out_lines, species, entry["dataset"], database):
            raise ValueError("Template must include a 'LogFilePath' line.")
        if not _replace_line(out_lines, "QNslabel", "  ".join(labels)):
            raise ValueError("Template must include a 'QNslabel' line.")
        if not _replace_line(out_lines, "QNsformat", "  ".join(formats)):
            raise ValueError("Template must include a 'QNsformat' line.")
        if not _replace_line(
            out_lines,
            "QNsFilter(Y/N)",
            f"N          {'  '.join(f'{label}[]' for label in labels)}",
        ):
            raise ValueError("Template must include a 'QNsFilter(Y/N)' line.")
        _update_unc_filter(out_lines, entry["has_uncertainty"], unc_default=UNC_DEFAULT)

        filename = _output_name(
            database=database,
            species=species,
            isotopologue=entry["isotopologue"],
            dataset=entry["dataset"],
            multi_files=multi_files,
            multi_iso=multi_iso,
            multi_dataset=multi_dataset,
        )
        out_path = output_dir / filename
        out_path.write_text("".join(out_lines))
        created.append(out_path)

    return created


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate ExoMol/ExoAtom/HITRAN xsec input files."
    )
    parser.add_argument(
        "--database",
        default=None,
        choices=["ExoMol", "ExoAtom", "HITRAN"],
        help="Database to target (ExoMol, ExoAtom, or HITRAN).",
    )
    parser.add_argument(
        "--molecule",
        default=None,
        help="Molecule name for ExoMol/HITRAN (e.g., CO, H2O)."
    )
    parser.add_argument(
        "--template",
        default=None,
        help="Template input file to copy and update.",
    )
    parser.add_argument(
        "--output-dir",
        default="input",
        help="Output directory for generated .inp files.",
    )
    parser.add_argument(
        "--species-id",
        type=int,
        default=0,
        help=(
            "SpeciesID to set in generated files (defaults to 0 for ExoMol/ExoAtom). "
            "Ignored for HITRAN, where SpeciesID is resolved automatically per "
            "isotopologue."
        ),
    )
    args = parser.parse_args()

    database = args.database
    if database is None:
        raise ValueError("Database must be specified with --database.")
    
    species = args.molecule
    if species is None:
        raise ValueError("Molecule/atom name must be specified with --molecule.")

    template = args.template
    if template is None:
        if database == "ExoMol":
            template = "input/templates/CO_ExoMol_template.inp"
        elif database == "ExoAtom":
            template = "input/templates/Na_ExoAtom_template.inp"
        elif database == "HITRAN":
            template = "input/templates/H2O_HITRAN_template.inp"
        else:
            raise ValueError("Template must be specified for unknown database.")

    if database == "HITRAN" and args.species_id:
        print("Note: --species-id is ignored for HITRAN; SpeciesID is resolved per isotopologue.")

    repo_root = Path(__file__).resolve().parents[1]
    data_dir = (repo_root / 'data' / 'lbl' / database.lower()).resolve()
    template_path = (repo_root / template).resolve()
    output_dir = (repo_root / args.output_dir).resolve()

    print("Generating xsec input files with the following parameters:")
    print(f"\tDatabase: {database}")
    print(f"\tSpecies: {species}")
    print(f"\tData Directory: {data_dir}")
    print(f"\tTemplate: {template}")
    print(f"\tOutput Directory: {args.output_dir}")
    if database != "HITRAN":
        print(f"\tSpecies ID: {args.species_id}")

    if database == "HITRAN":
        created = generate_hitran_inputs(
            species=species,
            data_dir=data_dir,
            template_path=template_path,
            output_dir=output_dir,
        )
    else:
        created = generate_inputs(
            database=database,
            species=species,
            data_dir=data_dir,
            template_path=template_path,
            output_dir=output_dir,
            species_id=args.species_id,
        )
    for path in created:
        print("Created input file:\n\t{}".format(path))


if __name__ == "__main__":
    main()
