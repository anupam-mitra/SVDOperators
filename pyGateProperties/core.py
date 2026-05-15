import numpy as np

def gate_singular_values(gate_matrix: np.ndarray, d: int) -> np.ndarray:
    """Singular values of a two-site gate under the Schmidt bipartition."""

    expected_shape = (d * d, d * d)
    if gate_matrix.shape != expected_shape:
        raise ValueError(f"expected a {expected_shape} two-site gate matrix, got {gate_matrix.shape}")

    gate_tensor = gate_matrix.reshape((d, d, d, d))
    gate_tensor_transposed = gate_tensor.transpose((0, 2, 1, 3))
    gate_tensor_transposed_matrix = gate_tensor_transposed.reshape(d * d, d * d)
    return np.linalg.svdvals(gate_tensor_transposed_matrix)
