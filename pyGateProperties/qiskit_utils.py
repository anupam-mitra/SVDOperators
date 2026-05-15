import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator
from pyGateProperties.core import gate_singular_values

def circuit_matrix(circuit: QuantumCircuit) -> np.ndarray:
    return Operator.from_circuit(circuit).to_matrix()

def circuit_singular_values(circuit: QuantumCircuit) -> np.ndarray:
    return gate_singular_values(circuit_matrix(circuit), 2)
