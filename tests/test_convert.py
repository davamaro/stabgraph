import random

import numpy as np
import pytest

import stabgraph
from stabgraph import gf2


def _graph_state_stabilizers(n, seed, edge_probability=0.25):
    rng = random.Random(seed)
    adjacency = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < edge_probability:
                adjacency[i][j] = 1
                adjacency[j][i] = 1

    stabs = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append("X")
            elif adjacency[i][j]:
                row.append("Z")
            else:
                row.append("I")
        stabs.append("".join(row))
    return stabs


def _assert_valid_output(result, n):
    G, control, target, z_qubits, R = result
    assert G.shape == (n, n)
    assert np.array_equal(G, G.T)
    assert np.array_equal(np.diag(G), np.zeros(n, dtype=int))
    assert set(control).isdisjoint(set(target))
    assert sorted(control + target) == list(range(n))
    assert set(z_qubits).issubset(set(control))
    assert R.shape == (n, n)


def test_bell_example():
    stabs = ["XX", "ZZ"]
    result = stabgraph.convert(stabs)
    _assert_valid_output(result, 2)
    assert result[1] == [0]
    assert result[2] == [1]
    assert result[3] == []


def test_ghz_example():
    stabs = ["XXX", "ZZI", "IZZ"]
    result = stabgraph.convert(stabs, [0], [1])
    _assert_valid_output(result, 3)
    assert result[1] == [0]
    assert result[2] == [1, 2]


def test_steane_examples():
    stabs = [
        "XXXXIII",
        "IXXIXXI",
        "IIXXIXX",
        "ZZZZIII",
        "IZZIZZI",
        "IIZZIZZ",
        "ZZZZZZZ",
    ]
    result = stabgraph.convert(stabs, control=[0], shuffle=True)
    _assert_valid_output(result, 7)
    assert set(result[1]).issuperset({0})


def test_gf2_backend_is_known():
    assert gf2.BACKEND in {"galois", "numpy"}


@pytest.mark.parametrize(
    ("stabs", "message"),
    [
        ([], "non-empty"),
        (["", ""], "at least one qubit"),
        (["XX", "ZZ", "XX"], "exactly N independent generators"),
        (["XX", "Z"], "same length"),
        (["X_", "ZZ"], "Pauli symbols"),
        (["XX", "XZ"], "do not commute"),
    ],
)
def test_validation_errors(stabs, message):
    with pytest.raises(ValueError, match=message):
        stabgraph.convert(stabs)


def test_bad_partitions_raise_clean_errors():
    stabs = ["XXX", "ZZI", "IZZ"]
    with pytest.raises(ValueError, match="empty intersection"):
        stabgraph.convert(stabs, control=[0], target=[0])
    with pytest.raises(ValueError, match="duplicates"):
        stabgraph.convert(stabs, control=[0, 0])
    with pytest.raises(ValueError, match="0 to N-1"):
        stabgraph.convert(stabs, control=[-1])


def test_randomized_bell_runs_repeat_cleanly():
    stabs = ["XX", "ZZ"]
    for seed in range(50):
        random.seed(seed)
        result = stabgraph.convert(stabs, shuffle=True)
        _assert_valid_output(result, 2)


def test_graph_state_round_trip_shape():
    stabs = _graph_state_stabilizers(12, seed=9, edge_probability=0.35)
    result = stabgraph.convert(stabs)
    _assert_valid_output(result, 12)


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_large_generated_graph_state_90(seed):
    stabs = _graph_state_stabilizers(90, seed=14)
    random.seed(seed)
    result = stabgraph.convert(stabs, shuffle=True)
    _assert_valid_output(result, len(stabs))


@pytest.mark.parametrize("seed", [0, 17])
def test_large_generated_graph_state_225(seed):
    stabs = _graph_state_stabilizers(225, seed=49, edge_probability=0.18)
    random.seed(seed)
    result = stabgraph.convert(stabs, shuffle=True)
    _assert_valid_output(result, len(stabs))


def test_many_randomized_calls_on_medium_instance():
    stabs = _graph_state_stabilizers(25, seed=21, edge_probability=0.2)
    for seed in range(100):
        random.seed(seed)
        result = stabgraph.convert(stabs, shuffle=True)
        _assert_valid_output(result, len(stabs))
