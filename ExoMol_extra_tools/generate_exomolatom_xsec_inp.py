#!/usr/bin/env python3
"""
Generate ExoMol/ExoAtom cross-section input files from definition metadata.

Usage:
  python ExoMol_extra_tools/generate_exomol_xsec_inp.py --molecule CO
  python ExoMol_extra_tools/generate_exomol_xsec_inp.py --database ExoAtom --atom Na
"""
from __future__ import annotations

import argparse
import json
import os
import re
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


def _update_unc_filter(lines: List[str], has_uncertainty: bool) -> None:
    value = "Y" if has_uncertainty else "N"
    for idx, line in enumerate(lines):
        if line.lstrip().startswith("UncFilter(Y/N)"):
            lines[idx] = f"{'UncFilter(Y/N)':<40}{value}          0.01\n"
            return


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
        _update_unc_filter(out_lines, entry["has_uncertainty"])

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
        description="Generate ExoMol/ExoAtom xsec input files from definition metadata."
    )
    parser.add_argument(
        "--database",
        default="ExoMol",
        choices=["ExoMol", "ExoAtom"],
        help="Database to target (ExoMol or ExoAtom).",
    )
    parser.add_argument("--molecule", help="Molecule name for ExoMol (e.g., CO, H2O).")
    parser.add_argument("--atom", help="Atom name for ExoAtom (e.g., Na, Al).")
    parser.add_argument(
        "--data-dir",
        default="data",
        help="Base data directory containing ExoMol/ExoAtom definition files.",
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
        help="SpeciesID to set in generated files (defaults to 0 for ExoMol/ExoAtom).",
    )
    args = parser.parse_args()

    database = args.database
    if database == "ExoAtom":
        species = args.atom or args.molecule
        if not species:
            raise ValueError("Provide --atom for ExoAtom (or use --molecule as an alias).")
    else:
        species = args.molecule
        if not species:
            raise ValueError("Provide --molecule for ExoMol.")

    template = args.template
    if template is None:
        template = "input/templates/CO_ExoMol_template.inp" if database == "ExoMol" else "input/templates/Na_ExoAtom_template.inp"

    print("Generating xsec input files with the following parameters:")
    print(f"\tDatabase: {database}")
    print(f"\tSpecies: {species}")
    print(f"\tData Directory: {args.data_dir}")
    print(f"\tTemplate: {template}")
    print(f"\tOutput Directory: {args.output_dir}")
    print(f"\tSpecies ID: {args.species_id}")

    repo_root = Path(__file__).resolve().parents[1]
    data_dir = (repo_root / args.data_dir).resolve()
    template_path = (repo_root / template).resolve()
    output_dir = (repo_root / args.output_dir).resolve()

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
