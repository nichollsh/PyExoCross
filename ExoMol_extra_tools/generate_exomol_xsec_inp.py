#!/usr/bin/env python3
"""
Generate ExoMol cross-section input files from def.json metadata.

Usage:
  python scripts/generate_exomol_xsec_inp.py --molecule CO
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


def _replace_line(lines: List[str], key: str, value: str) -> None:
    for idx, line in enumerate(lines):
        if line.lstrip().startswith(key):
            lines[idx] = f"{key:<40}{value}\n"
            return


def _update_unc_filter(lines: List[str], has_uncertainty: bool) -> None:
    value = "Y" if has_uncertainty else "N"
    for idx, line in enumerate(lines):
        if line.lstrip().startswith("UncFilter(Y/N)"):
            lines[idx] = f"{'UncFilter(Y/N)':<40}{value}          0.01\n"
            return


def _find_def_files(data_dir: Path, molecule: str) -> List[Path]:
    molecule_lower = molecule.lower()
    matches: List[Path] = []
    for root, _dirs, files in os.walk(data_dir, followlinks=True):
        if not files:
            continue
        root_parts = Path(root).parts
        
        if "exomol" not in root_parts:
            continue
        
        try:
            mol_part = root_parts[root_parts.index("exomol") + 1]
        except (ValueError, IndexError):
            continue

        if mol_part.lower() != molecule_lower:
            continue

        print("Found directory for molecule '{}'\n\t{}".format(molecule, root))

        for fname in files:
            if fname.endswith(".def.json"):
                matches.append(Path(root) / fname)
    return sorted(matches)


def _output_name(
    molecule: str,
    isotopologue: str,
    dataset: str,
    multi_files: bool,
    multi_iso: bool,
    multi_dataset: bool,
) -> str:
    if not multi_files:
        return f"{molecule}_ExoMol_xsec.inp"
    if multi_iso and multi_dataset:
        return f"{molecule}_ExoMol_{isotopologue}_{dataset}_xsec.inp"
    if multi_iso:
        return f"{molecule}_ExoMol_{isotopologue}_xsec.inp"
    return f"{molecule}_ExoMol_{dataset}_xsec.inp"


def generate_inputs(
    molecule: str,
    data_dir: Path,
    template_path: Path,
    output_dir: Path,
    species_id: int,
) -> List[Path]:
    def_files = _find_def_files(data_dir, molecule)
    if not def_files:
        raise FileNotFoundError(f"No def.json files found for molecule '{molecule}'.")

    meta = []
    for def_path in def_files:
        with def_path.open() as handle:
            payload = json.load(handle)
        meta.append(
            {
                "path": def_path,
                "isotopologue": payload["isotopologue"]["iso_slug"],
                "dataset": payload["dataset"]["name"],
                "fields": payload["dataset"]["states"]["states_file_fields"],
                "has_uncertainty": bool(payload["dataset"]["states"].get("uncertainties_available", False)),
            }
        )

    isotopologues = {m["isotopologue"] for m in meta}
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
        _replace_line(out_lines, "Database", "ExoMol")
        _replace_line(out_lines, "Molecule", molecule)
        _replace_line(out_lines, "Isotopologue", entry["isotopologue"])
        _replace_line(out_lines, "Dataset", entry["dataset"])
        _replace_line(out_lines, "SpeciesID", str(species_id))
        _replace_line(out_lines, "QNslabel", "  ".join(labels))
        _replace_line(out_lines, "QNsformat", "  ".join(formats))
        _replace_line(
            out_lines,
            "QNsFilter(Y/N)",
            f"N          {'  '.join(f'{label}[]' for label in labels)}",
        )
        _update_unc_filter(out_lines, entry["has_uncertainty"])

        filename = _output_name(
            molecule=molecule,
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
        description="Generate ExoMol xsec input files from def.json metadata."
    )
    parser.add_argument("--molecule", required=True, help="Molecule name (e.g., CO, H2O).")
    parser.add_argument(
        "--data-dir",
        default="data",
        help="Base data directory containing ExoMol def.json files.",
    )
    parser.add_argument(
        "--template",
        default="input/CO_ExoMol_xsec.inp",
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
        help="SpeciesID to set in generated files (defaults to 0 for ExoMol).",
    )
    args = parser.parse_args()

    print("Generating ExoMol xsec input files with the following parameters:")
    print(f"\tMolecule: {args.molecule}")
    print(f"\tData Directory: {args.data_dir}")
    print(f"\tTemplate: {args.template}")
    print(f"\tOutput Directory: {args.output_dir}")
    print(f"\tSpecies ID: {args.species_id}")

    repo_root = Path(__file__).resolve().parents[1]
    data_dir = (repo_root / args.data_dir).resolve()
    template_path = (repo_root / args.template).resolve()
    output_dir = (repo_root / args.output_dir).resolve()

    created = generate_inputs(
        molecule=args.molecule,
        data_dir=data_dir,
        template_path=template_path,
        output_dir=output_dir,
        species_id=args.species_id,
    )
    for path in created:
        print("Created input file:\n\t{}".format(path))


if __name__ == "__main__":
    main()
