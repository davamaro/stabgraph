#!/usr/bin/env python3

import numpy as np

try:
    import galois
except Exception as exc:  # pragma: no cover - exercised in environments without a working galois install
    galois = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None


BACKEND = "galois" if galois is not None else "numpy"

if galois is not None:
    GF2 = galois.GF(2)


def _as_uint8_matrix(A):
    matrix = np.asarray(A, dtype=np.uint8)
    if matrix.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    return matrix


def _row_reduce_numpy(A):
    A = _as_uint8_matrix(A).copy()
    m, n = A.shape
    if m == 0 or n == 0:
        raise ValueError("empty matrix")

    pivot_row = 0
    for col in range(n):
        pivot = None
        for row in range(pivot_row, m):
            if A[row, col]:
                pivot = row
                break
        if pivot is None:
            continue
        if pivot != pivot_row:
            A[[pivot_row, pivot]] = A[[pivot, pivot_row]]
        mask = A[:, col].astype(bool)
        mask[pivot_row] = False
        A[mask] ^= A[pivot_row]
        pivot_row += 1
        if pivot_row == m:
            break
    return A


def _rank_numpy(A):
    A = _as_uint8_matrix(A)
    if A.size == 0 or A.shape[0] == 0 or A.shape[1] == 0:
        return 0
    reduced = _row_reduce_numpy(A)
    return int(np.count_nonzero(np.any(reduced, axis=1)))


def _inverse_numpy(A):
    A = _as_uint8_matrix(A)
    m, n = A.shape
    if m != n:
        raise ValueError("a non-square matrix does not have an inverse")

    augmented = np.concatenate((A.copy(), np.eye(n, dtype=np.uint8)), axis=1)
    reduced = _row_reduce_numpy(augmented)
    left = reduced[:, :n]
    if not np.array_equal(left, np.eye(n, dtype=np.uint8)):
        raise ValueError("singular matrix")
    return reduced[:, n:]


def gauss(A):
    if galois is None:
        return _row_reduce_numpy(A)

    field = GF2(_as_uint8_matrix(A))
    m, n = field.shape
    if m == 0 or n == 0:
        raise ValueError("empty matrix")
    return np.asarray(field.row_reduce(), dtype=np.uint8)


def rank(A):
    if galois is None:
        return _rank_numpy(A)

    field = GF2(_as_uint8_matrix(A))
    return int(field.row_space().shape[0])


def inverse(A):
    if galois is None:
        return _inverse_numpy(A)

    field = GF2(_as_uint8_matrix(A))
    m, n = field.shape
    if m != n:
        raise ValueError("a non-square matrix does not have an inverse")
    try:
        return np.asarray(np.linalg.inv(field), dtype=np.uint8)
    except np.linalg.LinAlgError as exc:
        raise ValueError("singular matrix") from exc
