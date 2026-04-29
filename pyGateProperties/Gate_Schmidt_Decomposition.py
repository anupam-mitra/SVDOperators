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
from qiskit import QuantumCircuit
from qiskit.circuit.library import RXXGate
from qiskit.quantum_info import Operator

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
def gate_singular_values(gate_matrix: np.ndarray, d: int) -> np.ndarray:
    """Singular values of a two-site gate under the Schmidt bipartition."""

    expected_shape = (d * d, d * d)
    if gate_matrix.shape != expected_shape:
        raise ValueError(f"expected a {expected_shape} two-site gate matrix, got {gate_matrix.shape}")

    gate_tensor = gate_matrix.reshape((d, d, d, d))
    gate_tensor_transposed = gate_tensor.transpose((0, 2, 1, 3))
    gate_tensor_transposed_matrix = gate_tensor_transposed.reshape(d * d, d * d)
    return np.linalg.svdvals(gate_tensor_transposed_matrix)

# %%
def circuit_matrix(circuit: QuantumCircuit) -> np.ndarray:
    return Operator.from_circuit(circuit).to_matrix()

# %%
def circuit_singular_values(circuit: QuantumCircuit) -> np.ndarray:
    return gate_singular_values(circuit_matrix(circuit), 2)

# %%
def composite_circuit() -> QuantumCircuit:
    circuit = QuantumCircuit(2)
    circuit.t(0)
    circuit.cx(0, 1)
    circuit.t(1)
    circuit.cx(1, 0)
    return circuit

# %%
def cz_circuit() -> QuantumCircuit:
    circuit = QuantumCircuit(2)
    circuit.cz(0, 1)
    return circuit

# %%
def iswap_circuit() -> QuantumCircuit:
    circuit = QuantumCircuit(2)
    circuit.iswap(0, 1)
    return circuit

# %%
def product_gate_circuit() -> QuantumCircuit:
    circuit = QuantumCircuit(2)
    circuit.t(0)
    circuit.s(1)
    return circuit

# %%
def swap_circuit() -> QuantumCircuit:
    circuit = QuantumCircuit(2)
    circuit.swap(0, 1)
    return circuit

# %%
def rxx_gate_matrix(theta: float) -> np.ndarray:
    return Operator(RXXGate(theta)).data

# %%
def computational_basis_index(a: int, b: int, d: int) -> int:
    return a * d + b


# %%
def sum_gate_matrix(d: int) -> np.ndarray:
    matrix = np.zeros((d * d, d * d), dtype=complex)
    for a in range(d):
        for b in range(d):
            input_index = computational_basis_index(a, b, d)
            output_index = computational_basis_index(a, (a + b) % d, d)
            matrix[output_index, input_index] = 1.0
    return matrix


# %%
def cz_gate_matrix(d: int) -> np.ndarray:
    matrix = np.zeros((d * d, d * d), dtype=complex)
    omega = np.exp(2j * np.pi / d)
    for a in range(d):
        for b in range(d):
            index = computational_basis_index(a, b, d)
            matrix[index, index] = omega ** (a * b)
    return matrix


# %%
def swap_gate_matrix(d: int) -> np.ndarray:
    matrix = np.zeros((d * d, d * d), dtype=complex)
    for a in range(d):
        for b in range(d):
            input_index = computational_basis_index(a, b, d)
            output_index = computational_basis_index(b, a, d)
            matrix[output_index, input_index] = 1.0
    return matrix


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
