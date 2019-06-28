# stabgraph

Every stabilizer state can be transformed into a graph state by means of a local 
Clifford unitary. stabgraph is a python package that, given a stabilizer state,
it finds such a graph state and the local Clifford unitary that transforms one
into the other.  

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install foobar.

```bash
pip install foobar
```

## Usage

```python
import stabgraph

G , c , t , z , R = stabgraph.convert(S,c=[],t=[],shuffle=False) 
```
INPUTS
S       contains the N stabilizer operators defining a stabilizer state of N 
        qubits. It is a list of N strings 'PPPP...', one for each stabilizer 
        operator. Every string has N elements 'P' from the set 'I', 'X', 'Y', 
        'Z', that represents the set of Pauli matrices.

OPTIONAL INPUTS
c       is a list of control qubits. Gives the option to set some qubits as 
        control qubits. Qubits are labelled from 0 to N-1. It is an empy list by
        default.
t       is a list of target qubits. Gives the option to set some qubits as
        target qubits. Qubits are labelled from 0 to N-1. It is an empty list by
        default.
shuffle can be set to be True. For a given stabilizer state there are multiple 
        local Clifford equivalent graph states that can be obtained by this 
        program. If shuffle=True the output graph is one of these graphs chosen 
        randomly. shuffle=False by default.
        
OUTPUTS
G       adjacency matrix defining the underlying graph of the graph state. It is
        a NxN numpy array composed by 0 and 1.
c       list of control qubits labelled from 0 to N-1. The program completes the
        list of control qubits given as an input.
t       list of target qubits labelled from 0 to N-1. The program completes the
        list of target qubits given as an input. A Hadamard gate is applied on
        every target qubit.
z       list of control qubits labelled from 0 to N-1 where a pi/2 z-rotation is
        applied.
R       invertible binary matrix representing the recombination of stabilizers
        performed to obtain the stabilizers of the graph state. It is a NxN 
        numpy array composed of 0 and 1.
        
EXAMPLES
Bell pair
```
S = ['XX','ZZ']
G , c , t , z , R = stabgraph.convert(S)
```
GHZ state fixing 0 as a control qubit and 1 as a target qubit
```
S = ['XXX','ZZI','IZZ']
G , c , t , z , R = stabgraph.convert(S,[0],[1])
```
Steane code in the |0> logical state. Shuffle=True to obtain different outputs
```
S = ['XXXXIII','IXXIXXI','IIXXIXX','ZZZZIII','IZZIZZI','IIZZIZZ','ZZZZZZZ']
G , c , t , z , R = stabgraph.convert(S,shuffle=True)
```




        

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)