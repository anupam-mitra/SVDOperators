import numpy as np
from qiskit.circuit.library import RXXGate
from qiskit.quantum_info import Operator

def rxx_gate_matrix(theta: float) -> np.ndarray:
    return Operator(RXXGate(theta)).data
