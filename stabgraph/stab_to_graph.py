#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 15:32:08 2019

@author: davamaro
"""

import numpy as np
import random
from . import gauss_binary as gb


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
    # number of qubits N and number of stabilizers Ns    
    N = len(stabs[0])
    Ns = len(stabs)
    # binary representation
    A = np.array([[stabs[j][i] in {'Y', 'Z'} for j in range(N)] for i in range(Ns)] +
                 [[stabs[j][i] in {'Y', 'X'} for j in range(N)] for i in range(Ns)], dtype=int)
    # raise Exception if there are not enough stabilizers
    if Ns < N:
        raise Exception('The number of stabilizers in S can not be smaller than the number of qubits')
    # raise Exception if rank(p) is not N
    if gb.rank(gb.gauss(A)) != N:
        raise Exception('S must contain the same number of independent stabilizers than qubits')
    # raise Exception if stabilizers do not commute
    for i in range(Ns-1):
        for j in range(i+1, Ns):
            if A[:N, i].dot(A[N:, j]) % 2 != A[:N, j].dot(A[N:, i]) % 2:
                raise Exception('generators ' + stabs[i] + ' and ' + stabs[j] + ' do not commute')
    # raise Exception if control and target qubits have non empty intersection
    if set(control).intersection(set(target)):
        raise Exception('c and t must have empty intersection')
    # if a control or target qubit has a label bigger than N
    for i in control:
        if i >= N:
            raise Exception('control qubits must be labelled from 0 to N-1')
    for i in target:
        if i >= N:
            raise Exception('target qubits must be labelled from 0 to N-1')
    # shuffle rows
    if shuffle:
        qubits = control + random.sample(set(range(N))-(set(control).union(set(target))), N-len(control)-len(target)) +\
                                                                                          target
    else:
        qubits = control + list(set(range(N))-(set(control).union(set(target)))) + target
    # reorder rows in A
    A = np.array([A[qubits[i]] for i in range(N)] + [A[qubits[i]+N] for i in range(N)])
    # put A in the form that allows to make Gauss elimination in the right way add the identity to monitor the
    # recombinations performed by the Gaussian elimination
    A = np.array(list(A[N:, :])+list(A[:N, :])+list(np.eye(N, dtype=int))).T
    # perform the Gaussian elimination and identify the matrix R that performs that recombination
    A = gb.gauss(A)
    R = A[:, -N:].T
    A = A[:, :-N]
    n = gb.rank(A[:, :N])
    if len(control) > n:
        raise Exception('too many control qubits selected')
    if len(target) > N - n:
        raise Exception('too many target qubits selected')
    if gb.rank(A[:n, :len(control)]) < len(control):
        raise Exception('wrong selection of control qubits')    
    # select control and target qubits 
    for i in range(len(control), n):
        for j in range(i, N):
            if A[i, j] and qubits[j] not in target and qubits[j] not in control:
                control.append(qubits[j])
                break
            elif qubits[j] not in target:
                target.append(qubits[j])
    target = target + list(set(range(N))-set(control).union(set(target)))
    if len(control) != n:
        raise Exception('wrong selection of control and/or target qubits')    
    # put A back to the original form
    A = A.T
    A = np.array(list(A[N:, :])+list(A[:N, :]))
    A = np.array([A[qubits.index(i)] for i in control] + [A[qubits.index(i)] for i in target] +
                 [A[qubits.index(i)+N] for i in control] + [A[qubits.index(i)+N] for i in target])
    # update the order of the qubits in the rows
    qubits = control+target
    # build xlc and the inverse
    xlc = A[N:n+N, :n]
    if n != 0:
        xlc_inv = gb.inverse(xlc)
    else:
        xlc_inv = np.zeros((0, 0), dtype=int)
    # identify the rest of blocks
    xlt = A[n+N:, :n]
    zlc = A[:n, :n]
    zlt = A[n:N, :n]
    zrt = A[n:N, n:]
    # compute r
    if n != N:
        zrt_inv = gb.inverse(zrt)
    else:
        zrt_inv = np.zeros((0, 0), dtype=int)
    R = R.dot(np.block([[xlc_inv, np.zeros((n, N-n), dtype=int)], [zrt_inv.dot(zlt).dot(xlc_inv) % 2, zrt_inv]])) % 2
    # compute C and B and obtain the list z of the qubits where the z-rotation is performed
    B = xlt.dot(xlc_inv) % 2
    C = (zlc+B.T.dot(zlt)).dot(xlc_inv) % 2
    z = [qubits[i] for i in range(n) if C[i, i]]
    C = (C + np.diag(np.diagonal(C))) % 2
    # Adjacency matrix
    G = np.block([[C, B.T], [B, np.zeros((N-n, N-n), dtype=int)]])
    G = np.array([G[qubits.index(i)] for i in range(N)])
    G = np.array([G[:, qubits.index(i)] for i in range(N)]).T
    return G, sorted(control), sorted(target), sorted(z), R
