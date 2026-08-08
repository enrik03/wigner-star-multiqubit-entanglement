#!/usr/bin/env python3
"""Verifica la extension multiqubit con una entrada Wigner explicita GHZ4.

La formula de entrada no se obtiene reconstruyendo primero una matriz densidad.
En la base de simbolos de Pauli, los coeficientes no nulos de W_GHZ4 son:

* cadenas formadas por I y Z con un numero par de Z;
* cadenas formadas por X y Y con un numero par de Y, con signo (-1)^(n_Y/2).

Todos los coeficientes anteriores llevan el factor 1/16. La reconstruccion
matricial se usa unicamente como comprobacion independiente y para el SDP.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from verify_detector_jerarquico import (
    TOL,
    summarize_example,
)


def ghz4_wigner_coefficients() -> np.ndarray:
    """Coeficientes exactos de W_GHZ4 en la base e_0,e_x,e_y,e_z."""
    coeffs = np.zeros((4, 4, 4, 4), dtype=complex)
    for index in itertools.product(range(4), repeat=4):
        if all(axis in (0, 3) for axis in index):
            z_count = sum(axis == 3 for axis in index)
            if z_count % 2 == 0:
                coeffs[index] = 1.0 / 16.0
        elif all(axis in (1, 2) for axis in index):
            y_count = sum(axis == 2 for axis in index)
            if y_count % 2 == 0:
                coeffs[index] = (-1.0) ** (y_count // 2) / 16.0
    return coeffs


def validate_ghz4_report(report: dict) -> dict:
    failures: list[str] = []
    validation = report["validation"]
    if not validation["valid_density_operator"]:
        failures.append("La entrada Wigner GHZ4 no reconstruye un estado valido.")
    if abs(validation["purity"] - 1.0) > 1.0e-8:
        failures.append("La pureza GHZ4 no es uno.")
    if len(report["cut_layer"]) != 7:
        failures.append("No se evaluaron las siete biparticiones no redundantes de cuatro qubits.")
    if not report["all_cuts_npt"]:
        failures.append("GHZ4 no fue identificado como NPT en todos los cortes.")
    for cut in report["cut_layer"]:
        if abs(cut["f_star"] + 0.25) > 2.0e-8:
            failures.append(f"F_star inesperado en el corte {cut['cut']}.")
        if cut["max_star_matrix_moment_error"] > 2.0e-8:
            failures.append(f"Error star-matriz excesivo en el corte {cut['cut']}.")
    global_layer = report["gme_moment_layer"]
    if not global_layer["map_certifies_gme"]:
        failures.append("El mapa modificado no certifica GME para GHZ4.")
    if not any(level["certifies_gme"] for level in global_layer["levels"]):
        failures.append("Ninguna matriz de Hankel finita certifica GME para GHZ4.")
    if global_layer["max_star_matrix_moment_error"] > 2.0e-8:
        failures.append("Los momentos del mapa no coinciden entre Wigner-star y matrices.")
    sdp = report["ppt_mixture_sdp"]
    if sdp.get("available") and not sdp.get("certifies_not_ppt_mixture"):
        failures.append("El SDP disponible no certifica que GHZ4 este fuera de las mezclas PPT.")
    return {
        "passed": not failures,
        "tolerance": 2.0e-8,
        "failures": failures,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parents[1] / "results" / "four_qubit_extension.json",
    )
    parser.add_argument("--skip-sdp", action="store_true")
    args = parser.parse_args()

    report = summarize_example(
        "GHZ4",
        ghz4_wigner_coefficients(),
        "modified",
        run_sdp=not args.skip_sdp,
    )
    payload = {
        "conventions": {
            "primary_input": "explicit_GHZ4_Wigner_coefficients",
            "kernel": "Delta=(I+sqrt(3)n.sigma)/2",
            "trace_functional": "T_4[f]=16 integral f",
            "nonredundant_bipartitions": 7,
        },
        "example": report,
    }
    payload["validation_gate"] = validate_ghz4_report(report)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if not payload["validation_gate"]["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
