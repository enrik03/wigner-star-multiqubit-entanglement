#!/usr/bin/env python3
"""Comprueba que el manuscrito y los resultados reproducibles siguen sincronizados."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


ROOT = Path(__file__).parents[1]


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def first_hankel_label(example: dict) -> str:
    for level in example["gme_moment_layer"]["levels"]:
        if level["certifies_gme"]:
            return f"H_{level['level']}"
    return "--"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results" / "manuscript_sync_check.json",
    )
    args = parser.parse_args()

    tex = (ROOT / "main.tex").read_text(encoding="utf-8")
    bib = (ROOT / "references.bib").read_text(encoding="utf-8")
    three = load_json(ROOT / "results" / "detector_jerarquico_checks.json")
    four = load_json(ROOT / "results" / "four_qubit_extension.json")
    examples = {item["name"]: item for item in three["examples"]}
    examples["GHZ4"] = four["example"]

    failures: list[str] = []

    expected_rows = {
        "GHZ3": ("GHZ$_3$", "GME"),
        "W3": ("W$_3$", "GME"),
        "GHZ4": ("GHZ$_4$", "GME"),
        "mezcla_biseparable_Bell": (r"$W_{\mathrm{bs}}$", "biseparable"),
        "GHZ3_ruidoso_mu_0.50": ("GHZ$_3(0.50)$", "GME por SDP"),
        "W3_ruidoso_mu_0.90": ("W$_3(0.90)$", r"GME; Hankel $\ell\leq3$ inconcluso"),
    }
    for name, (label, conclusion) in expected_rows.items():
        example = examples[name]
        f_star = example["cut_layer"][0]["f_star"]
        map_min = example["gme_moment_layer"]["map_min_eigenvalue"]
        hankel = first_hankel_label(example)
        ng = example["ppt_mixture_sdp"].get("genuine_negativity", 0.0)
        f_text = f"{f_star:.6f}"
        map_text = f"{map_min:+.6f}"
        ng_text = "$0$" if abs(ng) < 5.0e-7 else f"${ng:.6f}$"
        row_pattern = re.compile(
            re.escape(label)
            + r"\s*&\s*\$"
            + re.escape(f_text)
            + r"\$\s*&\s*\$"
            + re.escape(map_text)
            + r"\$\s*&\s*"
            + (r"--" if hankel == "--" else re.escape(f"${hankel}$"))
            + r"\s*&\s*"
            + re.escape(ng_text)
            + r"\s*&\s*"
            + re.escape(conclusion),
        )
        if not row_pattern.search(tex):
            failures.append(f"La fila de {name} no coincide con los JSON.")

    thresholds = three["noise_thresholds"]
    threshold_snippets = {
        "noisy_GHZ3_modified_map": [
            "0.733333334",
            "0.934266784",
            "0.733333335",
            "0.733333334",
        ],
        "noisy_W3_transposition_map": [
            "0.898868744",
            "no detecta",
            "0.953490622",
            "0.900474281",
        ],
    }
    for family, values in threshold_snippets.items():
        source = thresholds[family]
        numerical = [value for value in source.values() if value is not None]
        for expected, actual in zip((item for item in values if item != "no detecta"), numerical):
            if expected != f"{actual:.9f}":
                failures.append(f"Umbral no sincronizado para {family}: {expected}.")
        if any(value not in tex for value in values):
            failures.append(f"La tabla de umbrales no contiene todos los valores de {family}.")

    four_gate = four.get("validation_gate", {})
    if not four_gate.get("passed"):
        failures.append("La puerta GHZ4 no pasa.")
    if len(four["example"]["cut_layer"]) != 7:
        failures.append("La salida GHZ4 no contiene siete cortes.")

    cited = {
        key.strip()
        for group in re.findall(r"\\cite\{([^}]+)\}", tex)
        for key in group.split(",")
    }
    bib_keys = set(re.findall(r"@[A-Za-z]+\s*\{\s*([^,]+),", bib))
    missing = sorted(cited - bib_keys)
    unused = sorted(bib_keys - cited)
    if missing:
        failures.append("Citas sin entrada BibTeX: " + ", ".join(missing))
    if unused:
        failures.append("Entradas BibTeX no citadas: " + ", ".join(unused))

    forbidden = re.findall(r"TODO|FIXME|not_implemented|\\usepackage\[utf8\]\{inputenc\}", tex)
    if forbidden:
        failures.append("Marcadores o paquetes obsoletos en main.tex: " + ", ".join(forbidden))

    required = [
        r"\label{eq:W-GHZ4}",
        r"\Tfun_N[\mathcal M_T^{(N)}W]",
        r"\lambda_{\min}(H_2)=-0.6337602",
        "Zhang2022PTMoments",
        "python scripts/verify_four_qubit_extension.py",
    ]
    for snippet in required:
        if snippet not in tex:
            failures.append(f"Falta contenido requerido en main.tex: {snippet}")

    payload = {
        "passed": not failures,
        "checked_examples": sorted(expected_rows),
        "cited_keys": len(cited),
        "bib_entries": len(bib_keys),
        "ghz4_nonredundant_cuts": len(four["example"]["cut_layer"]),
        "failures": failures,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
