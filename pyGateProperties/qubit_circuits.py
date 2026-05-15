from qiskit import QuantumCircuit

def composite_circuit() -> QuantumCircuit:
    circuit = QuantumCircuit(2)
    circuit.t(0)
    circuit.cx(0, 1)
    circuit.t(1)
    circuit.cx(1, 0)
    return circuit

def cz_circuit() -> QuantumCircuit:
    circuit = QuantumCircuit(2)
    circuit.cz(0, 1)
    return circuit

def iswap_circuit() -> QuantumCircuit:
    circuit = QuantumCircuit(2)
    circuit.iswap(0, 1)
    return circuit

def product_gate_circuit() -> QuantumCircuit:
    circuit = QuantumCircuit(2)
    circuit.t(0)
    circuit.s(1)
    return circuit

def swap_circuit() -> QuantumCircuit:
    circuit = QuantumCircuit(2)
    circuit.swap(0, 1)
    return circuit
