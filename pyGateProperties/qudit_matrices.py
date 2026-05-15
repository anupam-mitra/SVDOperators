import numpy as np

def computational_basis_index(a: int, b: int, d: int) -> int:
    return a * d + b


def sum_gate_matrix(d: int) -> np.ndarray:
    matrix = np.zeros((d * d, d * d), dtype=complex)
    for a in range(d):
        for b in range(d):
            input_index = computational_basis_index(a, b, d)
            output_index = computational_basis_index(a, (a + b) % d, d)
            matrix[output_index, input_index] = 1.0
    return matrix


def cz_gate_matrix(d: int) -> np.ndarray:
    matrix = np.zeros((d * d, d * d), dtype=complex)
    omega = np.exp(2j * np.pi / d)
    for a in range(d):
        for b in range(d):
            index = computational_basis_index(a, b, d)
            matrix[index, index] = omega ** (a * b)
    return matrix


def swap_gate_matrix(d: int) -> np.ndarray:
    matrix = np.zeros((d * d, d * d), dtype=complex)
    for a in range(d):
        for b in range(d):
            input_index = computational_basis_index(a, b, d)
            output_index = computational_basis_index(b, a, d)
            matrix[output_index, input_index] = 1.0
    return matrix
