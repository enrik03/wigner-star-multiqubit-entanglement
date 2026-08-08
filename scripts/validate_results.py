#!/usr/bin/env python3
"""Puerta de integridad para detector_jerarquico_checks.json."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "json_path",
        nargs="?",
        type=Path,
        default=Path(__file__).parents[1] / "results" / "detector_jerarquico_checks.json",
    )
    args = parser.parse_args()
    payload = json.loads(args.json_path.read_text(encoding="utf-8"))
    by_name = {example["name"]: example for example in payload["examples"]}
    failures: list[str] = []
    max_error = 0.0

    for example in payload["examples"]:
        if not example["validation"]["valid_density_operator"]:
            failures.append(f"entrada W no fisica: {example['name']}")
        max_error = max(max_error, example["validation"]["max_coefficient_roundtrip_error"])
        for cut in example["cut_layer"]:
            max_error = max(max_error, cut["max_star_matrix_moment_error"])
        layer = example["gme_moment_layer"]
        max_error = max(
            max_error,
            layer["max_symbol_matrix_map_error"],
            layer["max_star_matrix_moment_error"],
        )

    ghz = by_name["GHZ3"]
    w_state = by_name["W3"]
    biseparable = by_name["mezcla_biseparable_Bell"]
    noisy_ghz = by_name["GHZ3_ruidoso_mu_0.50"]

    if not ghz["gme_moment_layer"]["levels"][0]["certifies_gme"]:
        failures.append("H1 debe certificar GHZ3 con el mapa modificado")
    if w_state["gme_moment_layer"]["levels"][0]["certifies_gme"]:
        failures.append("H1 no debe certificar W3")
    if not w_state["gme_moment_layer"]["levels"][1]["certifies_gme"]:
        failures.append("H2 debe certificar W3")
    if not biseparable["all_cuts_npt"]:
        failures.append("el contraejemplo biseparable debe ser NPT en todos los cortes")
    if biseparable["gme_moment_layer"]["map_certifies_gme"]:
        failures.append("falso positivo del mapa GME sobre la mezcla biseparable")
    if any(level["certifies_gme"] for level in biseparable["gme_moment_layer"]["levels"]):
        failures.append("falso positivo Hankel sobre la mezcla biseparable")
    if biseparable["ppt_mixture_sdp"].get("certifies_not_ppt_mixture", False):
        failures.append("falso positivo del SDP sobre la mezcla biseparable")
    if not noisy_ghz["ppt_mixture_sdp"].get("certifies_not_ppt_mixture", False):
        failures.append("el SDP debe certificar GHZ3 ruidoso con mu=0.50")
    if noisy_ghz["gme_moment_layer"]["map_certifies_gme"]:
        failures.append("el mapa modificado no debe detectar GHZ3 ruidoso con mu=0.50")
    for pure_name in ("GHZ3", "W3"):
        if not by_name[pure_name]["ppt_mixture_sdp"].get("certifies_not_ppt_mixture", False):
            failures.append(f"el SDP debe certificar {pure_name}")
    if max_error > 2e-10:
        failures.append(f"error star/matriz excesivo: {max_error:.3e}")

    thresholds = payload["noise_thresholds"]
    ghz_threshold = thresholds["noisy_GHZ3_modified_map"]["map_threshold"]
    w_threshold = thresholds["noisy_W3_transposition_map"]["map_threshold"]
    if not (0.732 < ghz_threshold < 0.735):
        failures.append("umbral del mapa GHZ fuera del intervalo de control")
    if not (0.897 < w_threshold < 0.901):
        failures.append("umbral del mapa W fuera del intervalo de control")

    report = {
        "passed": not failures,
        "max_independent_representation_error": max_error,
        "failures": failures,
    }
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
