#!/usr/bin/env python3
"""Verificacion reproducible del detector Wigner-star jerarquico.

La entrada primaria de cada ejemplo es el tensor de coeficientes del simbolo

    W(Omega_1,...,Omega_N) = sum_alpha c_alpha e_alpha(Omega_1,...,Omega_N),

con e_0=1 y e_i=sqrt(3)n_i.  El operador se reconstruye solo para realizar
comprobaciones independientes y para resolver el SDP de mezclas PPT.

Capas verificadas:
  L1: F_star por corte mediante reflexion, potencias star e identidades de Newton.
  L2: momentos Wigner-star de mapas GME y matrices de Hankel.
  L3: negatividad genuina via testigos completamente descomponibles (SDP).
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Iterable

import numpy as np


TOL = 1.0e-9
PAULI = (
    np.eye(2, dtype=complex),
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
)

# e_mu star e_nu = sum_gamma C[mu,nu,gamma] e_gamma.
LOCAL_STAR = np.zeros((4, 4, 4), dtype=complex)
LOCAL_STAR[0, 0, 0] = 1.0
for i in range(1, 4):
    LOCAL_STAR[0, i, i] = 1.0
    LOCAL_STAR[i, 0, i] = 1.0
    LOCAL_STAR[i, i, 0] = 1.0
for i, j, k in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
    LOCAL_STAR[i, j, k] = 1j
    LOCAL_STAR[j, i, k] = -1j


def kron_all(factors: Iterable[np.ndarray]) -> np.ndarray:
    out = np.array([[1.0]], dtype=complex)
    for factor in factors:
        out = np.kron(out, factor)
    return out


def operator_from_wigner_coefficients(coeffs: np.ndarray) -> np.ndarray:
    """Inverse SW map in the Pauli-symbol basis."""
    n_qubits = coeffs.ndim
    dim = 2**n_qubits
    out = np.zeros((dim, dim), dtype=complex)
    for index in itertools.product(range(4), repeat=n_qubits):
        if abs(coeffs[index]) > 0:
            out += coeffs[index] * kron_all(PAULI[i] for i in index)
    return out


def coefficients_from_operator(operator: np.ndarray, n_qubits: int) -> np.ndarray:
    dim = 2**n_qubits
    coeffs = np.zeros((4,) * n_qubits, dtype=complex)
    for index in itertools.product(range(4), repeat=n_qubits):
        sigma = kron_all(PAULI[i] for i in index)
        coeffs[index] = np.trace(operator @ sigma) / dim
    return coeffs


def identity_state_coefficients(n_qubits: int) -> np.ndarray:
    coeffs = np.zeros((4,) * n_qubits, dtype=complex)
    coeffs[(0,) * n_qubits] = 1.0 / 2**n_qubits
    return coeffs


def ghz3_wigner_coefficients() -> np.ndarray:
    """Coeficientes leidos directamente de la formula explicita de W_GHZ."""
    c = np.zeros((4, 4, 4), dtype=complex)
    terms = {
        (0, 0, 0): 1,
        (3, 3, 0): 1,
        (3, 0, 3): 1,
        (0, 3, 3): 1,
        (1, 1, 1): 1,
        (1, 2, 2): -1,
        (2, 1, 2): -1,
        (2, 2, 1): -1,
    }
    for index, value in terms.items():
        c[index] = value / 8.0
    return c


def w3_wigner_coefficients() -> np.ndarray:
    """Coeficientes de W_W3 en la base e_0,e_x,e_y,e_z."""
    c = np.zeros((4, 4, 4), dtype=complex)
    c[0, 0, 0] = 1.0 / 8.0
    for q in range(3):
        index = [0, 0, 0]
        index[q] = 3
        c[tuple(index)] += 1.0 / 24.0
    for i, j in itertools.combinations(range(3), 2):
        index = [0, 0, 0]
        index[i] = index[j] = 3
        c[tuple(index)] += -1.0 / 24.0
        for axis in (1, 2):
            index = [0, 0, 0]
            index[i] = index[j] = axis
            c[tuple(index)] += 1.0 / 12.0
    c[3, 3, 3] = -1.0 / 8.0
    for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        for axis in (1, 2):
            index = [0, 0, 0]
            index[i] = index[j] = axis
            index[k] = 3
            c[tuple(index)] += 1.0 / 12.0
    return c


def bell_phi_plus_coefficients() -> np.ndarray:
    c = np.zeros((4, 4), dtype=complex)
    c[0, 0] = c[1, 1] = c[3, 3] = 1.0 / 4.0
    c[2, 2] = -1.0 / 4.0
    return c


def zero_qubit_coefficients() -> np.ndarray:
    return np.array([0.5, 0.0, 0.0, 0.5], dtype=complex)


def biseparable_bell_mixture_coefficients() -> np.ndarray:
    """Mezcla simetrica de Bell en cada par y |0> en el tercer qubit."""
    bell = bell_phi_plus_coefficients()
    zero = zero_qubit_coefficients()
    out = np.zeros((4, 4, 4), dtype=complex)
    for pair in ((0, 1), (0, 2), (1, 2)):
        spectator = ({0, 1, 2} - set(pair)).pop()
        term = np.zeros_like(out)
        for i, j, k in itertools.product(range(4), repeat=3):
            index = (i, j, k)
            term[index] = bell[index[pair[0]], index[pair[1]]] * zero[index[spectator]]
        out += term / 3.0
    return out


def white_noise_mix(coeffs: np.ndarray, visibility: float) -> np.ndarray:
    return visibility * coeffs + (1.0 - visibility) * identity_state_coefficients(coeffs.ndim)


def reflect_subset(coeffs: np.ndarray, subset: tuple[int, ...]) -> np.ndarray:
    out = np.array(coeffs, copy=True)
    for index in itertools.product(range(4), repeat=coeffs.ndim):
        if sum(index[q] == 2 for q in subset) % 2:
            out[index] *= -1.0
    return out


def conjugate_x(coeffs: np.ndarray, subset: tuple[int, ...]) -> np.ndarray:
    """Conjugacion local X O X: sigma_y,z cambian de signo."""
    out = np.array(coeffs, copy=True)
    for index in itertools.product(range(4), repeat=coeffs.ndim):
        if sum(index[q] in (2, 3) for q in subset) % 2:
            out[index] *= -1.0
    return out


def star_multiply(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    n_qubits = left.ndim
    out = np.zeros_like(left)
    nonzero_left = [(idx, left[idx]) for idx in np.ndindex(left.shape) if abs(left[idx]) > 1e-15]
    nonzero_right = [(idx, right[idx]) for idx in np.ndindex(right.shape) if abs(right[idx]) > 1e-15]
    for a, ca in nonzero_left:
        for b, cb in nonzero_right:
            targets = [
                [(g, LOCAL_STAR[a[q], b[q], g]) for g in range(4) if LOCAL_STAR[a[q], b[q], g] != 0]
                for q in range(n_qubits)
            ]
            for factors in itertools.product(*targets):
                index = tuple(g for g, _ in factors)
                out[index] += ca * cb * np.prod([value for _, value in factors])
    return out


def star_power_moments(coeffs: np.ndarray, max_order: int) -> list[float]:
    current = np.array(coeffs, copy=True)
    dim = 2**coeffs.ndim
    moments: list[float] = []
    for order in range(1, max_order + 1):
        if order > 1:
            current = star_multiply(current, coeffs)
        value = np.real_if_close(dim * current[(0,) * coeffs.ndim])
        moments.append(float(np.real(value)))
    return moments


def elementary_from_power_sums(power_sums: list[float]) -> list[float]:
    e = [1.0]
    for k in range(1, len(power_sums) + 1):
        value = sum(
            (-1) ** (j - 1) * e[k - j] * power_sums[j - 1]
            for j in range(1, k + 1)
        )
        e.append(value / k)
    return e[1:]


def unique_bipartitions(n_qubits: int) -> list[tuple[int, ...]]:
    parties = set(range(n_qubits))
    cuts: list[tuple[int, ...]] = []
    seen: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    for size in range(1, n_qubits):
        for subset in itertools.combinations(range(n_qubits), size):
            complement = tuple(sorted(parties - set(subset)))
            key = tuple(sorted((tuple(subset), complement)))
            if key not in seen:
                seen.add(key)
                cuts.append(tuple(subset))
    return cuts


def partial_transpose_matrix(operator: np.ndarray, subset: tuple[int, ...], n_qubits: int) -> np.ndarray:
    tensor = operator.reshape([2] * (2 * n_qubits))
    for q in subset:
        tensor = np.swapaxes(tensor, q, n_qubits + q)
    return tensor.reshape(operator.shape)


def fstar_cut(coeffs: np.ndarray, subset: tuple[int, ...]) -> dict:
    reflected = reflect_subset(coeffs, subset)
    q_star = star_power_moments(reflected, 2**coeffs.ndim)
    e = elementary_from_power_sums(q_star)
    operator = operator_from_wigner_coefficients(coeffs)
    pt = partial_transpose_matrix(operator, subset, coeffs.ndim)
    eigenvalues = np.linalg.eigvalsh(pt)
    q_matrix = [float(np.sum(eigenvalues**n)) for n in range(1, len(eigenvalues) + 1)]
    return {
        "cut": "".join(chr(ord("A") + q) for q in subset),
        "q_star": q_star,
        "e_newton": e,
        "f_star": float(min(e)),
        "pt_eigenvalues": [float(x) for x in eigenvalues],
        "negativity": float(-np.sum(eigenvalues[eigenvalues < -TOL])),
        "max_star_matrix_moment_error": float(np.max(np.abs(np.array(q_star) - q_matrix))),
    }


def gme_transposition_map_coefficients(coeffs: np.ndarray, modified_x: bool = False) -> np.ndarray:
    n_qubits = coeffs.ndim
    out = np.zeros_like(coeffs)
    for subset in unique_bipartitions(n_qubits):
        term = reflect_subset(coeffs, subset)
        if modified_x:
            term = conjugate_x(term, subset)
        out += term
    c_n = (2 ** (n_qubits - 1) - 2) / 2.0
    trace_input = 2**n_qubits * coeffs[(0,) * n_qubits].real
    out[(0,) * n_qubits] += c_n * trace_input
    return out


def gme_transposition_map_matrix(operator: np.ndarray, modified_x: bool = False) -> np.ndarray:
    n_qubits = int(round(np.log2(operator.shape[0])))
    out = np.zeros_like(operator)
    for subset in unique_bipartitions(n_qubits):
        term = partial_transpose_matrix(operator, subset, n_qubits)
        if modified_x:
            local = [PAULI[1] if q in subset else PAULI[0] for q in range(n_qubits)]
            unitary = kron_all(local)
            term = unitary @ term @ unitary.conj().T
        out += term
    c_n = (2 ** (n_qubits - 1) - 2) / 2.0
    out += c_n * np.trace(operator) * np.eye(operator.shape[0])
    return out


def hankel_matrix(moments: list[float], level: int) -> np.ndarray:
    if len(moments) < 2 * level + 1:
        raise ValueError("Se requieren momentos hasta orden 2*level+1")
    return np.array(
        [[moments[i + j] for j in range(level + 1)] for i in range(level + 1)],
        dtype=float,
    )


def gme_moment_layer(coeffs: np.ndarray, modified_x: bool, max_level: int = 3) -> dict:
    mapped = gme_transposition_map_coefficients(coeffs, modified_x=modified_x)
    moments = star_power_moments(mapped, 2 * max_level + 1)
    mapped_operator_from_w = operator_from_wigner_coefficients(mapped)
    source_operator = operator_from_wigner_coefficients(coeffs)
    mapped_operator_direct = gme_transposition_map_matrix(source_operator, modified_x=modified_x)
    levels = []
    for level in range(1, max_level + 1):
        hankel = hankel_matrix(moments, level)
        eigvals = np.linalg.eigvalsh(hankel)
        levels.append(
            {
                "level": level,
                "moments_used": 2 * level + 1,
                "min_hankel_eigenvalue": float(eigvals[0]),
                "hankel_determinant": float(np.linalg.det(hankel)),
                "g_star_hankel": float(max(0.0, -eigvals[0])),
                "certifies_gme": bool(eigvals[0] < -TOL),
            }
        )
    mapped_eigvals = np.linalg.eigvalsh(mapped_operator_direct)
    matrix_moments = [
        float(np.trace(np.linalg.matrix_power(mapped_operator_direct, n)).real)
        for n in range(1, 2 * max_level + 2)
    ]
    return {
        "map": "modified_transposition" if modified_x else "transposition",
        "map_min_eigenvalue": float(mapped_eigvals[0]),
        "map_certifies_gme": bool(mapped_eigvals[0] < -TOL),
        "map_moments_star": moments,
        "levels": levels,
        "max_symbol_matrix_map_error": float(np.max(np.abs(mapped_operator_from_w - mapped_operator_direct))),
        "max_star_matrix_moment_error": float(np.max(np.abs(np.array(moments) - matrix_moments))),
    }


def _index_bits(index: int, n_qubits: int) -> list[int]:
    return [int(bit) for bit in f"{index:0{n_qubits}b}"]


def _bits_index(bits: list[int]) -> int:
    return int("".join(str(bit) for bit in bits), 2)


def partial_transpose_expression(variable, subset: tuple[int, ...], n_qubits: int, cp):
    dim = 2**n_qubits
    rows = []
    for row in range(dim):
        row_bits = _index_bits(row, n_qubits)
        entries = []
        for col in range(dim):
            col_bits = _index_bits(col, n_qubits)
            source_row = row_bits.copy()
            source_col = col_bits.copy()
            for q in subset:
                source_row[q], source_col[q] = source_col[q], source_row[q]
            entries.append(variable[_bits_index(source_row), _bits_index(source_col)])
        rows.append(entries)
    return cp.bmat(rows)


def genuine_negativity_sdp(operator: np.ndarray) -> dict:
    try:
        import cvxpy as cp
    except ImportError as exc:
        return {"available": False, "reason": str(exc)}

    n_qubits = int(round(np.log2(operator.shape[0])))
    dim = operator.shape[0]
    witness = cp.Variable((dim, dim), hermitian=True)
    constraints = []
    eye = np.eye(dim)
    for subset in unique_bipartitions(n_qubits):
        p_var = cp.Variable((dim, dim), hermitian=True)
        q_var = cp.Variable((dim, dim), hermitian=True)
        constraints.extend(
            [
                p_var >> 0,
                eye - p_var >> 0,
                q_var >> 0,
                eye - q_var >> 0,
                witness == p_var + partial_transpose_expression(q_var, subset, n_qubits, cp),
            ]
        )
    objective = cp.Minimize(cp.real(cp.trace(witness @ operator)))
    problem = cp.Problem(objective, constraints)
    solver_used = None
    last_error = None
    for solver, options in (("CLARABEL", {}), ("SCS", {"eps": 1e-7, "max_iters": 100000})):
        try:
            problem.solve(solver=solver, verbose=False, **options)
            solver_used = solver
            if problem.status in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
                break
        except Exception as exc:  # pragma: no cover - solver fallback
            last_error = str(exc)
    if problem.value is None:
        return {"available": True, "status": problem.status, "error": last_error}
    value = float(np.real(problem.value))
    return {
        "available": True,
        "status": problem.status,
        "solver": solver_used,
        "optimal_witness_expectation": value,
        "genuine_negativity": float(max(0.0, -value)),
        "certifies_not_ppt_mixture": bool(value < -2e-6),
    }


def validate_wigner_input(coeffs: np.ndarray) -> dict:
    operator = operator_from_wigner_coefficients(coeffs)
    roundtrip = coefficients_from_operator(operator, coeffs.ndim)
    eigvals = np.linalg.eigvalsh(operator)
    return {
        "trace": float(np.trace(operator).real),
        "min_eigenvalue": float(eigvals[0]),
        "purity": float(np.trace(operator @ operator).real),
        "max_coefficient_roundtrip_error": float(np.max(np.abs(coeffs - roundtrip))),
        "valid_density_operator": bool(abs(np.trace(operator).real - 1.0) < TOL and eigvals[0] > -TOL),
    }


def summarize_example(name: str, coeffs: np.ndarray, map_kind: str, run_sdp: bool) -> dict:
    operator = operator_from_wigner_coefficients(coeffs)
    cuts = [fstar_cut(coeffs, subset) for subset in unique_bipartitions(coeffs.ndim)]
    global_layer = gme_moment_layer(coeffs, modified_x=(map_kind == "modified"))
    return {
        "name": name,
        "primary_input": "explicit_Wigner_coefficients",
        "validation": validate_wigner_input(coeffs),
        "cut_layer": cuts,
        "all_cuts_npt": bool(all(cut["f_star"] < -TOL for cut in cuts)),
        "gme_moment_layer": global_layer,
        "ppt_mixture_sdp": genuine_negativity_sdp(operator) if run_sdp else {"available": False, "reason": "--skip-sdp"},
    }


def first_detection_threshold(
    pure_coeffs: np.ndarray,
    modified_x: bool,
    detector,
    iterations: int = 55,
) -> float | None:
    low, high = 0.0, 1.0
    if not detector(gme_moment_layer(white_noise_mix(pure_coeffs, high), modified_x)):
        return None
    for _ in range(iterations):
        middle = (low + high) / 2.0
        result = gme_moment_layer(white_noise_mix(pure_coeffs, middle), modified_x)
        if detector(result):
            high = middle
        else:
            low = middle
    return high


def threshold_report() -> dict:
    def map_detected(result: dict) -> bool:
        return result["map_min_eigenvalue"] < -TOL

    def hankel_detected(level: int):
        return lambda result: result["levels"][level - 1]["min_hankel_eigenvalue"] < -TOL

    families = {
        "noisy_GHZ3_modified_map": (ghz3_wigner_coefficients(), True),
        "noisy_W3_transposition_map": (w3_wigner_coefficients(), False),
    }
    out = {}
    for name, (coeffs, modified) in families.items():
        out[name] = {
            "map_threshold": first_detection_threshold(coeffs, modified, map_detected),
            "H1_threshold": first_detection_threshold(coeffs, modified, hankel_detected(1)),
            "H2_threshold": first_detection_threshold(coeffs, modified, hankel_detected(2)),
            "H3_threshold": first_detection_threshold(coeffs, modified, hankel_detected(3)),
        }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path(__file__).parents[1] / "results" / "detector_jerarquico_checks.json")
    parser.add_argument("--skip-sdp", action="store_true")
    parser.add_argument("--skip-thresholds", action="store_true")
    args = parser.parse_args()

    examples = [
        ("GHZ3", ghz3_wigner_coefficients(), "modified"),
        ("W3", w3_wigner_coefficients(), "standard"),
        ("mezcla_biseparable_Bell", biseparable_bell_mixture_coefficients(), "standard"),
        ("GHZ3_ruidoso_mu_0.50", white_noise_mix(ghz3_wigner_coefficients(), 0.50), "modified"),
        ("W3_ruidoso_mu_0.90", white_noise_mix(w3_wigner_coefficients(), 0.90), "standard"),
    ]
    payload = {
        "conventions": {
            "kernel": "Delta=(I+sqrt(3)n.sigma)/2",
            "measure": "dmu=sin(theta)dtheta dphi/(4pi)",
            "trace_functional": "T_N[f]=2^N integral f",
            "hankel_score": "G_star,l=max(0,-lambda_min(H_l))",
            "ppt_mixture_score": "genuine negativity from fully decomposable witnesses",
        },
        "examples": [
            summarize_example(name, coeffs, map_kind, run_sdp=not args.skip_sdp)
            for name, coeffs, map_kind in examples
        ],
        "noise_thresholds": None if args.skip_thresholds else threshold_report(),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
