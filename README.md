# StabGraph

Every stabilizer state can be transformed into a graph state by means of a local 
Clifford unitary. stabgraph is contains the function `convert` that, given a 
stabilizer state, it finds such a graph state and the local Clifford unitary 
that transforms the stabilizer state into the graph state. This function follows
the steps described in the article [1].

[1] Manuscript under preparation. 

## Installation

To install using [pip](https://pip.pypa.io/en/stable/):

```bash
python -m pip install --upgrade pip
python -m pip install stabgraph
```

Replace `python` with `python3` as appropriate.


## Usage

```python
import stabgraph

G , c , t , z , R = stabgraph.convert(stabs,control=None,target=None,shuffle=False) 
```

INPUT

`stabs`     contains the N stabilizer operators defining a stabilizer state of N 
            qubits. It is a list of N strings 'PPPP...', one for each stabilizer 
            operator. Every string has N elements 'P' from the set 'I', 'X', 'Y', 
            'Z', that represents the set of Pauli matrices.

OPTIONAL INPUTS

`control`   is a list of control qubits. Gives the option to set some qubits as 
            control qubits. Qubits are labelled from 0 to N-1. It is an empy list by
            default.
        
`target`    is a list of target qubits. Gives the option to set some qubits as
            target qubits. Qubits are labelled from 0 to N-1. It is an empty list by
            default.
        
`shuffle`   can be set to be True. For a given stabilizer state there are multiple 
            local Clifford equivalent graph states that can be obtained by this 
            program. If shuffle=True the output graph is one of these graphs chosen 
            randomly. shuffle=False by default.
        
OUTPUTS

`G`       adjacency matrix defining the underlying graph of the graph state. It is
        a NxN numpy array composed by 0 and 1.
        
`c`       list of control qubits labelled from 0 to N-1. The program completes the
        list of control qubits given as an input.
        
`t`       list of target qubits labelled from 0 to N-1. The program completes the
        list of target qubits given as an input. A Hadamard gate is applied on
        every target qubit.
        
`z`       list of control qubits labelled from 0 to N-1 where a pi/2 z-rotation is
        applied.
        
`R`       invertible binary matrix representing the recombination of stabilizers
        performed to obtain the stabilizers of the graph state. It is a NxN 
        numpy array composed of 0 and 1.
        
EXAMPLES

Bell pair
```
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
```
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
```
>>> stabs = ['XXXXIII','IXXIXXI','IIXXIXX','ZZZZIII','IZZIZZI','IIZZIZZ','ZZZZZZZ']
>>> G , c , t , z , R = stabgraph.convert(stabs, control = [0], shuffle=True)
>>> c
[0,1,2]
>>> G , c , t , z , R = stabgraph.convert(stabs, control = [2, 5], shuffle=True)
>>> c
[2,5,6]

```

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
