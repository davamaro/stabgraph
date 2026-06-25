# Changelog

## 0.1.4

- Added a GF(2) backend layer with an optional `galois` acceleration path and a
  safe NumPy fallback.
- Removed the old `gauss_binary.py` module.
- Reduced Python-level overhead in Pauli-to-binary conversion and qubit
  reordering inside `convert()`.
- Expanded tests with repeated randomized calls and backend coverage.
- Added initial PyPI publishing documentation and workflow scaffolding.

## 0.1.3

- Added explicit input validation for `convert()`.
- Replaced an `IndexError` on over-specified generator inputs with clear `ValueError` messages.
- Added automated pytest coverage for examples, validation errors, and large generated regression cases.
- Added modern build metadata and GitHub Actions test automation.
- Improved README documentation and installation guidance.
