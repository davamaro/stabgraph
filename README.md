# StabGraph

Every stabilizer state can be transformed into a graph state by means of a local
Clifford unitary. `stabgraph` provides a single public function, `convert`, that
takes a stabilizer description of an `n`-qubit stabilizer state and returns a
graph-state representation together with the local Clifford bookkeeping needed to
recover the original state.

The package is useful when you want a graph-state form as an intermediate step
for circuit synthesis or state-preparation workflows.

[1] Manuscript under preparation.

## Installation

Install from PyPI:

```bash
python -m pip install stabgraph
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
algebra over GF(2). Notes for a future accelerated backend live in
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
