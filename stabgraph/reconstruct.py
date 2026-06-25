#!/usr/bin/env python3

import numpy as np

from . import gf2
from .stab_to_graph import _binary_symplectic_matrix


_PAULI_PRODUCT = {
    ("I", "I"): (1, "I"),
    ("I", "X"): (1, "X"),
    ("I", "Y"): (1, "Y"),
    ("I", "Z"): (1, "Z"),
    ("X", "I"): (1, "X"),
    ("X", "X"): (1, "I"),
    ("X", "Y"): (1, "Z"),
    ("X", "Z"): (-1, "Y"),
    ("Y", "I"): (1, "Y"),
    ("Y", "X"): (-1, "Z"),
    ("Y", "Y"): (1, "I"),
    ("Y", "Z"): (1, "X"),
    ("Z", "I"): (1, "Z"),
    ("Z", "X"): (1, "Y"),
    ("Z", "Y"): (-1, "X"),
    ("Z", "Z"): (1, "I"),
}


def graph_state_generators(adjacency):
    n = adjacency.shape[0]
    stabs = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append("X")
            elif adjacency[i, j]:
                row.append("Z")
            else:
                row.append("I")
        stabs.append("".join(row))
    return stabs


def recombine_generators(generators, recombinations):
    """
    Recombines a generator list according to a binary matrix.

    ``generators`` may be either plain Pauli strings or signed generators of the
    form ``(sign, pauli_string)``. The return value always uses signed generator
    pairs ``(sign, pauli_string)``.
    """

    if generators and isinstance(generators[0], str):
        signed_generators = [(1, stab) for stab in generators]
    else:
        signed_generators = list(generators)

    recombinations = np.asarray(recombinations, dtype=np.uint8) % 2
    identity = "I" * len(signed_generators[0][1])
    out = []
    for column in range(recombinations.shape[1]):
        sign = 1
        stab = identity
        for row in np.flatnonzero(recombinations[:, column]):
            local_sign, stab = _multiply_pauli_strings(stab, signed_generators[row][1])
            sign *= local_sign * signed_generators[row][0]
        out.append((sign, stab))
    return out


def _multiply_pauli_strings(left, right):
    sign = 1
    out = []
    for pauli_left, pauli_right in zip(left, right):
        local_sign, pauli = _PAULI_PRODUCT[(pauli_left, pauli_right)]
        sign *= local_sign
        out.append(pauli)
    return sign, "".join(out)


def _apply_hadamard(stab, qubit):
    chars = list(stab)
    pauli = chars[qubit]
    chars[qubit] = {"I": "I", "X": "Z", "Y": "Y", "Z": "X"}[pauli]
    sign = -1 if pauli == "Y" else 1
    return sign, "".join(chars)


def _apply_s_gate(stab, qubit):
    chars = list(stab)
    pauli = chars[qubit]
    if pauli == "X":
        chars[qubit] = "Y"
        sign = 1
    elif pauli == "Y":
        chars[qubit] = "X"
        sign = -1
    else:
        sign = 1
    return sign, "".join(chars)


def reconstruct_generators(graph, target=(), z=()):
    """
    Returns signed stabilizer generators obtained from the graph-state generators
    after applying H on the target qubits and S on the qubits in z.

    The return value is a list of ``(sign, pauli_string)`` pairs where
    ``sign`` is either ``+1`` or ``-1``.
    """

    signed_stabs = [(1, stab) for stab in graph_state_generators(graph)]

    for qubit in target:
        updated = []
        for sign, stab in signed_stabs:
            local_sign, new_stab = _apply_hadamard(stab, qubit)
            updated.append((sign * local_sign, new_stab))
        signed_stabs = updated

    for qubit in z:
        updated = []
        for sign, stab in signed_stabs:
            local_sign, new_stab = _apply_s_gate(stab, qubit)
            updated.append((sign * local_sign, new_stab))
        signed_stabs = updated

    return signed_stabs


def same_binary_stabilizer_group(left, right):
    left_matrix = _binary_symplectic_matrix(left, len(left))
    right_matrix = _binary_symplectic_matrix(right, len(right))
    combined = np.concatenate((left_matrix, right_matrix), axis=1)
    return gf2.rank(combined) == gf2.rank(left_matrix) == gf2.rank(right_matrix)


def _solve_full_rank_linear_system(matrix, rhs):
    matrix = np.asarray(matrix, dtype=np.uint8) % 2
    rhs = np.asarray(rhs, dtype=np.uint8) % 2
    if rhs.ndim == 1:
        rhs = rhs[:, np.newaxis]

    m, n = matrix.shape
    augmented = np.concatenate((matrix.copy(), rhs.copy()), axis=1)
    pivot_rows = []
    pivot_row = 0

    for col in range(n):
        pivot = next((row for row in range(pivot_row, m) if augmented[row, col]), None)
        if pivot is None:
            raise ValueError("matrix must have full column rank")
        if pivot != pivot_row:
            augmented[[pivot_row, pivot]] = augmented[[pivot, pivot_row]]
        for row in range(m):
            if row != pivot_row and augmented[row, col]:
                augmented[row] ^= augmented[pivot_row]
        pivot_rows.append(pivot_row)
        pivot_row += 1

    return augmented[pivot_rows, n:]


def infer_generator_phase_signs(reference_stabs, signed_generators):
    """
    Infers the +/-1 phase of each generator in ``reference_stabs`` relative to a
    signed generator set that spans the same binary stabilizer group.

    This function does not rely on phase data produced inside ``convert()``.
    Instead, it:
    1. compares binary symplectic representations,
    2. solves for each reference generator as a GF(2) combination of the
       provided signed generators, and
    3. recovers the resulting sign from signed Pauli multiplication.
    """

    unsigned_basis = [stab for _, stab in signed_generators]
    basis_matrix = _binary_symplectic_matrix(unsigned_basis, len(unsigned_basis))
    target_matrix = _binary_symplectic_matrix(reference_stabs, len(reference_stabs))
    coefficients = _solve_full_rank_linear_system(basis_matrix, target_matrix)

    resolved = recombine_generators(signed_generators, coefficients)

    phase_signs = []
    for column, (sign, stab) in enumerate(resolved):
        if stab != reference_stabs[column]:
            raise ValueError("reference generators are not in the binary span of the signed generators")
        phase_signs.append(sign)
    return phase_signs
