#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 15:32:08 2019

@author: davamaro
"""

import random
import numpy as np
from . import gf2


VALID_PAULIS = {'I', 'X', 'Y', 'Z'}


def _validate_stabilizers(stabs):
    if not stabs:
        raise ValueError('stabs must be a non-empty list of Pauli strings')

    N = len(stabs[0])
    if N == 0:
        raise ValueError('stabilizers must act on at least one qubit')
    if len(stabs) != N:
        raise ValueError('stabs must contain exactly N independent generators for N qubits')

    for stab in stabs:
        if len(stab) != N:
            raise ValueError('all stabilizers must have the same length')
        invalid = set(stab) - VALID_PAULIS
        if invalid:
            raise ValueError(
                "stabilizers may only contain Pauli symbols from {'I', 'X', 'Y', 'Z'}"
            )
    return N


def _binary_symplectic_matrix(stabs, N):
    encoded = np.frombuffer("".join(stabs).encode("ascii"), dtype=np.uint8).reshape(N, N)
    z_block = np.isin(encoded, (ord("Y"), ord("Z"))).T.astype(np.uint8, copy=False)
    x_block = np.isin(encoded, (ord("Y"), ord("X"))).T.astype(np.uint8, copy=False)
    return np.vstack((z_block, x_block))


def _validate_partition(control, target, N):
    if len(set(control)) != len(control):
        raise ValueError('control qubits must not contain duplicates')
    if len(set(target)) != len(target):
        raise ValueError('target qubits must not contain duplicates')
    if set(control).intersection(set(target)):
        raise ValueError('c and t must have empty intersection')

    for idx in control:
        if idx < 0 or idx >= N:
            raise ValueError('control qubits must be labelled from 0 to N-1')
    for idx in target:
        if idx < 0 or idx >= N:
            raise ValueError('target qubits must be labelled from 0 to N-1')


def convert(stabs, control=None, target=None, shuffle=False):
    """
    Returns a graph state that is local Clifford unitary equivalent to the given stabilizer state.

    :param stabs: list of N strings ['XYYXZ...','IIXXZ...',...,'stabilizerN'] defining the stabilizer operators
                  'stabilizer': string ['IXIYZ...'] of N elements from the set of Pauli operators 'I', 'X', 'Y', 'Z'
    :param control: list of control qubits [0,43,2,1,...] labelled from 0 to N-1
    :param target: list of target qubits [3,23,42,4...] labelled from 0 to N-1
    :param shuffle: True for obtaining a random output from the set of valid outputs.

    G, c, t, z, R = convert(s, control=None, target=None, shuffle=False)

    :return G: NxN numpy array composed only by 0s and 1s, that represents the adjacency matrix of the resulting graph
    :return c: final list of control qubits labelled from 0 to N-1. Contains the control qubits given in the input
               `control` and other control qubits assigned by the function.
    :return t: final list of target qubits labelled from 0 to N-1. Contains the target qubits given in the input
               `control` and other target qubits assigned by the function.
    :return z: list of control qubits where a pi/2 z-rotation is applied. It is a subset of the output c
    :return R: NxN numpy array composed only by 0s and 1s, that represents the recombinations of the stabilizers
               performed to obtain the stabilizers of the final graph state.
    """

    if control is None:
        control = []
    if target is None:
        target = []
    else:
        target = list(target)
    control = list(control)

    # number of qubits N and number of stabilizers Ns
    N = _validate_stabilizers(stabs)
    _validate_partition(control, target, N)

    # binary representation
    A = _binary_symplectic_matrix(stabs, N)
    # raise Exception if rank(p) is not N
    if gf2.rank(A) != N:
        raise ValueError('S must contain the same number of independent stabilizers as qubits')
    # raise Exception if stabilizers do not commute
    commutator = ((A[:N].T @ A[N:]) + (A[N:].T @ A[:N])) % 2
    if np.any(commutator):
        rows, cols = np.argwhere(commutator)
        i, j = next((i, j) for i, j in zip(rows.tolist(), cols.tolist()) if i < j)
        raise ValueError('generators ' + stabs[i] + ' and ' + stabs[j] + ' do not commute')
    # shuffle rows
    selected = set(control).union(target)
    remaining = [q for q in range(N) if q not in selected]
    if shuffle:
        qubits = control + random.sample(remaining, len(remaining)) + target
    else:
        qubits = control + remaining + target
    # reorder rows in A
    permutation = np.asarray(qubits, dtype=int)
    A = np.concatenate((A[permutation], A[permutation + N]), axis=0)
    # put A in the form that allows to make Gauss elimination in the right way add the identity to monitor the
    # recombinations performed by the Gaussian elimination
    A = np.concatenate((A[N:, :], A[:N, :], np.eye(N, dtype=np.uint8)), axis=0).T
    # perform the Gaussian elimination and identify the matrix R that performs that recombination
    A = gf2.gauss(A)
    R = A[:, -N:].T
    A = A[:, :-N]
    n = gf2.rank(A[:, :N])
    if len(control) > n:
        raise ValueError('too many control qubits selected')
    if len(target) > N - n:
        raise ValueError('too many target qubits selected')
    if gf2.rank(A[:n, :len(control)]) < len(control):
        raise ValueError('wrong selection of control qubits')
    # select control and target qubits
    for i in range(len(control), n):
        for j in range(i, N):
            if A[i, j] and qubits[j] not in target and qubits[j] not in control:
                control.append(qubits[j])
                break
            # elif qubits[j] not in target and qubits[j] not in control:
            #     target.append(qubits[j])
    target = target + [q for q in range(N) if q not in set(control).union(target)]  # add missing qubits
    if len(control) != n:
        raise ValueError('wrong selection of control and/or target qubits')
    # put A back to the original form
    A = A.T
    A = np.concatenate((A[N:, :], A[:N, :]), axis=0)
    qubit_positions = np.empty(N, dtype=int)
    qubit_positions[permutation] = np.arange(N)
    ordered_qubits = control + target
    ordered_positions = qubit_positions[ordered_qubits]
    A = np.concatenate(
        (
            A[ordered_positions],
            A[ordered_positions + N],
        ),
        axis=0,
    )
    # update the order of the qubits in the rows
    qubits = ordered_qubits
    # build xlc and the inverse
    xlc = A[N:n+N, :n]
    if n != 0:
        xlc_inv = gf2.inverse(xlc)
    else:
        xlc_inv = np.zeros((0, 0), dtype=int)
    # identify the rest of blocks
    xlt = A[n+N:, :n]
    zlc = A[:n, :n]
    zlt = A[n:N, :n]
    zrt = A[n:N, n:]
    # compute r
    if n != N:
        zrt_inv = gf2.inverse(zrt)
    else:
        zrt_inv = np.zeros((0, 0), dtype=int)
    R = R.dot(
        np.block(
            [
                [xlc_inv, np.zeros((n, N - n), dtype=np.uint8)],
                [zrt_inv.dot(zlt).dot(xlc_inv) % 2, zrt_inv],
            ]
        )
    ) % 2
    # compute C and B and obtain the list z of the qubits where the z-rotation is performed
    B = xlt.dot(xlc_inv) % 2
    C = (zlc+B.T.dot(zlt)).dot(xlc_inv) % 2
    z = [qubits[i] for i in range(n) if C[i, i]]
    C = (C + np.diag(np.diagonal(C))) % 2
    # Adjacency matrix
    G = np.block([[C, B.T], [B, np.zeros((N - n, N - n), dtype=np.uint8)]])
    final_positions = np.empty(N, dtype=int)
    final_positions[qubits] = np.arange(N)
    G = G[final_positions]
    G = G[:, final_positions]
    return G, sorted(control), sorted(target), sorted(z), R
