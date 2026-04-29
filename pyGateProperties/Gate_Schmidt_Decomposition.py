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
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import RXXGate
from qiskit.quantum_info import Operator

# %%
def gate_singular_values(gate_matrix: np.ndarray) -> np.ndarray:
    """Singular values of a two-qubit gate under the Schmidt bipartition."""

    if gate_matrix.shape != (4, 4):
        raise ValueError(f"expected a 4x4 two-qubit gate matrix, got {gate_matrix.shape}")

    gate_tensor = gate_matrix.reshape((2, 2, 2, 2))
    gate_tensor_transposed = gate_tensor.transpose((0, 2, 1, 3))
    gate_tensor_transposed_matrix = gate_tensor_transposed.reshape(4, 4)
    return np.linalg.svdvals(gate_tensor_transposed_matrix)

# %%
def circuit_matrix(circuit: QuantumCircuit) -> np.ndarray:
    return Operator.from_circuit(circuit).to_matrix()

# %%
def circuit_singular_values(circuit: QuantumCircuit) -> np.ndarray:
    return gate_singular_values(circuit_matrix(circuit))

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
def print_singular_values(label: str, singular_values: np.ndarray) -> None:
    print(label)
    print(np.array2string(np.real_if_close(singular_values), precision=10))
    print()

# %%
def main() -> None:
    fixed_gate_cases = [
        ("T I CXab I T CXba", circuit_singular_values(composite_circuit())),
        ("CZ", circuit_singular_values(cz_circuit())),
        ("iSWAP", circuit_singular_values(iswap_circuit())),
        ("Product gate", circuit_singular_values(product_gate_circuit())),
        ("SWAP", circuit_singular_values(swap_circuit())),
    ]

    for label, singular_values in fixed_gate_cases:
        print_singular_values(label, singular_values)

    for label, theta in [
        ("pi", np.pi),
        ("pi/2", np.pi / 2),
        ("pi/3", np.pi / 3),
        ("pi/4", np.pi / 4),
        ("pi/8", np.pi / 8),
    ]:
        print_singular_values(f"RXX({label})", gate_singular_values(rxx_gate_matrix(theta)))


# %%
if __name__ == "__main__":
    main()
