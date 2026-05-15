# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import sys

import numpy as np

from pyGateProperties.core import gate_singular_values
from pyGateProperties.qiskit_utils import circuit_singular_values
from pyGateProperties.qubit_circuits import composite_circuit, cz_circuit, iswap_circuit, product_gate_circuit, swap_circuit
from pyGateProperties.qubit_matrices import rxx_gate_matrix
from pyGateProperties.qudit_matrices import sum_gate_matrix, cz_gate_matrix, swap_gate_matrix


# %%
def parse_dimension(args: list[str]) -> int:
    if not args:
        return 2

    if len(args) != 1:
        raise ValueError("usage: python pyGateProperties/Gate_Schmidt_Decomposition.py [d]")

    try:
        d = int(args[0])
    except ValueError as exc:
        raise ValueError(f"expected integer local dimension d, got {args[0]}") from exc

    if d < 2:
        raise ValueError(f"expected local dimension d >= 2, got {d}")

    return d


# %%
def print_singular_values(label: str, singular_values: np.ndarray) -> None:
    print(label)
    print(np.array2string(np.real_if_close(singular_values), precision=10))
    print()

# %%
def qubit_fixed_gate_cases() -> list[tuple[str, np.ndarray]]:
    return [
        ("T I CXab I T CXba", circuit_singular_values(composite_circuit())),
        ("CZ", circuit_singular_values(cz_circuit())),
        ("iSWAP", circuit_singular_values(iswap_circuit())),
        ("Product gate", circuit_singular_values(product_gate_circuit())),
        ("SWAP", circuit_singular_values(swap_circuit())),
    ]


# %%
def qudit_fixed_gate_cases(d: int) -> list[tuple[str, np.ndarray]]:
    return [
        ("SUM (CX_d)", gate_singular_values(sum_gate_matrix(d), d)),
        ("CZ_d", gate_singular_values(cz_gate_matrix(d), d)),
        ("SWAP_d", gate_singular_values(swap_gate_matrix(d), d)),
    ]


# %%
def run_cases(cases: list[tuple[str, np.ndarray]]) -> None:
    for label, singular_values in cases:
        print_singular_values(label, singular_values)


# %%
def main(d: int) -> None:
    if d == 2:
        run_cases(qubit_fixed_gate_cases())

        for label, theta in [
            ("pi", np.pi),
            ("pi/2", np.pi / 2),
            ("pi/3", np.pi / 3),
            ("pi/4", np.pi / 4),
            ("pi/8", np.pi / 8),
        ]:
            print_singular_values(f"RXX({label})", gate_singular_values(rxx_gate_matrix(theta), 2))
        return

    print(f"d = {d}")
    print("Using standard qudit gate generalizations only.")
    print()
    run_cases(qudit_fixed_gate_cases(d))


# %%
if __name__ == "__main__":
    main(parse_dimension(sys.argv[1:]))
