# StabGraph

Every stabilizer state can be transformed into a graph state by means of a local
Clifford unitary. `stabgraph` provides a single public function, `convert`, that
takes a stabilizer description of an `n`-qubit stabilizer state and returns a
graph-state representation together with the local Clifford bookkeeping needed to
recover the original state.

The package is useful when you want a graph-state form as an intermediate step
for circuit synthesis or state-preparation workflows.
Internally, the binary linear algebra is carried out over GF(2). On supported
Python versions, `stabgraph` uses the maintained `galois` package for this
backend. If that backend is unavailable, it falls back to a bundled NumPy
implementation so the package still works.

[1] Design and Experimental Performance of Local Entanglement Witness Operators
https://arxiv.org/abs/1911.01144

## Installation

Install from PyPI:

```bash
python -m pip install stabgraph
```

For the accelerated backend on Python versions where `galois` is supported:

```bash
python -m pip install "stabgraph[accel]"
```

For local development:

```bash
git clone https://github.com/davamaro/stabgraph.git
cd stabgraph
python -m pip install -e .[test]
```

## Usage

```python
import stabgraph

G, c, t, z, R = stabgraph.convert(stabs, control=None, target=None, shuffle=False)
```

Reconstruct the signed post-processing generators and compare them with another
stabilizer description:

```python
signed_output = stabgraph.reconstruct_generators(G, t, z)
same_group = stabgraph.same_binary_stabilizer_group(
    [stab for _, stab in signed_output],
    stabs,
)
phase_signs = stabgraph.infer_generator_phase_signs(stabs, signed_output)
```

## Input contract

- `stabs` must be a non-empty list of exactly `n` independent commuting Pauli
  strings on `n` qubits.
- Every stabilizer string must have the same length and use only the symbols
  `I`, `X`, `Y`, and `Z`.
- `control` and `target` are optional lists of qubit labels in `0..n-1`.
- `control` and `target` must be disjoint and contain no duplicates.
- `shuffle=True` randomizes the internal qubit ordering before selecting a valid
  control/target partition. This can be useful when searching for lower-degree
  graph representations.

At the moment, `convert()` expects exactly `n` independent commuting generators on
`n` qubits. If you provide extra generators, the package now raises a clear error
instead of failing later in the algorithm.

## Return values

- `G` is the adjacency matrix of the resulting graph state.
- `c` is the completed list of control qubits.
- `t` is the completed list of target qubits.
- `z` is the subset of `c` where a pi/2 z-rotation is applied.
- `R` is the binary matrix describing the stabilizer recombinations used during
  the transformation.

## Reconstruction Helpers

The package also provides a few lightweight helper functions for checking and
post-processing the output of `convert()`:

- `reconstruct_generators(G, t, z)` builds the signed generators obtained after
  applying the returned local Cliffords to the graph-state generators.
- `same_binary_stabilizer_group(left, right)` checks whether two generator lists
  span the same binary stabilizer group over GF(2).
- `infer_generator_phase_signs(reference_stabs, signed_generators)` infers the
  `+/-1` phase of each reference generator relative to a signed generator basis,
  without requiring phase tracking inside `convert()`.
- `recombine_generators(generators, R)` applies a binary recombination matrix to
  either unsigned or signed generators.

## Examples

Bell pair
```python
>>> stabs = ['XX','ZZ']
>>> G , c , t , z , R = stabgraph.convert(stabs)
>>> G
np.array([[0,1],[1,0]])
>>> c
[0]
>>> t
[1]
>>> z
[]
>>> R
np.array([[1,0],[0,1]])
```

GHZ state fixing 0 as a control qubit and 1 as a target qubit
```python
>>> stabs = ['XXX','ZZI','IZZ']
>>> G , c , t , z , R = stabgraph.convert(stabs,[0],[1])
>>> c
[0]
>>> t
[1,2]
```

Steane code in the |0> logical state. 
Multiple graphs can be obtained, so put `shuffle=True` to obtain one of them randomly chosen.
The result respects the selection of control and target qubits
```python
>>> stabs = ['XXXXIII','IXXIXXI','IIXXIXX','ZZZZIII','IZZIZZI','IIZZIZZ','ZZZZZZZ']
>>> G , c , t , z , R = stabgraph.convert(stabs, control = [0], shuffle=True)
>>> c
[0,1,2]
>>> G , c , t , z , R = stabgraph.convert(stabs, control = [2, 5], shuffle=True)
>>> c
[2,5,6]
```

## Testing

Run the unit tests with:

```bash
python -m pytest
```

Continuous integration runs the test suite on supported Python versions through
GitHub Actions.

## Performance

Large stabilizer states can be slow because the heavy work is dense binary linear
algebra over GF(2). This release also removes a few obvious Python-level
bottlenecks in qubit reordering and Pauli-to-binary conversion. Notes for future
backend work live in
[docs/ACCELERATION_NOTES.md](docs/ACCELERATION_NOTES.md).

## Release notes

Recent user-facing changes are tracked in [CHANGELOG.md](CHANGELOG.md).

## Citation
```
@misc{amaro2019,
  author = 	"David Amaro",
  title = 	"StabGraph",
  year = 	"2019",
  month =   "July",
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/davamaro/stabgraph}}
}
```
