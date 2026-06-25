# Acceleration Notes

`stabgraph` spends most of its time in binary linear algebra over GF(2), especially
Gaussian elimination, rank, and matrix inversion.

## Current direction

The current backend is aimed at correctness and maintainability first:

- On supported Python versions, `galois` handles dense GF(2) row reduction,
  rank, and inversion.
- When `galois` is unavailable or not importable in the active interpreter,
  `stabgraph` falls back to a bundled NumPy-based GF(2) implementation.
- `stabgraph` keeps the higher-level graph-state logic and validation.

This is already a cleaner baseline than maintaining the original
`gauss_binary.py` module while still keeping the package usable in more runtime
environments.

## Other speed opportunities

- Reduce repeated Gaussian elimination on closely related matrices.
- Cache intermediate row-echelon forms when the same matrices are reused.
- Use boolean arrays consistently through the internal algebra pipeline.
- Benchmark the `galois` path against alternative dense GF(2) packages if
  Python 3.13 support becomes important for the accelerated path.
- Explore a sparse backend such as `ldpc` if very large low-density checks become
  the dominant workload.
- Benchmark large random graph-state and CSS-state inputs before and after any
  backend change.
