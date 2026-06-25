# Acceleration Notes

`stabgraph` spends most of its time in binary linear algebra over GF(2), especially
Gaussian elimination, rank, and matrix inversion. The current implementation keeps
that logic in pure NumPy plus Python loops, which is compact but not ideal for
large stabilizer states.

## Recommended direction

For a future speed-focused release, prefer an existing library backend for GF(2)
linear algebra instead of hand-optimizing the elimination code:

- `galois` is a good candidate for dense GF(2) matrices when it is available.
- `ldpc` is another candidate when sparse parity-check style workloads matter more.

The package should keep the current backend as the portable fallback and offer an
optional accelerated backend when a dedicated GF(2) library is installed.

## Why not in this patch

The current workspace does not have a dedicated GF(2) linear algebra package such
as `galois` installed, so this patch does not add an unverified optional backend.
Instead, it lays down tests and packaging so a backend swap can be made safely in
the next pass.

## Other speed opportunities

- Reduce repeated Gaussian elimination on closely related matrices.
- Cache intermediate row-echelon forms when the same matrices are reused.
- Use boolean arrays consistently through the internal algebra pipeline.
- Benchmark large random graph-state and CSS-state inputs before and after any
  backend change.
