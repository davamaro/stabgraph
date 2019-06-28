#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 15:32:08 2019

@author: davamaro
"""

import numpy as np
import random
from . import gauss_binary as gb


##### FINDS A GRAPH STATE THAT IS LOCAL CLIFFORD EQUIVALENT TO A GIVEN 
##### STABILIZER STATE ########################################################

# EXAMPLES
# Bell pair
#S = ['XX','ZZ']
# state |+0>
#S = ['XI','IZ']
# tensor product of Bell pair and |+>
#S = ['XXI','ZZI','IIX']
# GHZ state of three qubits
#S = ['XXX','ZZI','IZZ']
# three-qubit state local Clifford equivalent to the GHZ
#S = ['ZXX','XZX','XXZ']
# other three-qubit state local Clifford equivalent to the GHZ
#S = ['YZX','ZXX','ZZY']
# five-qubit code in the |0> logical state
#S = ['XZZXI','IXZZX','XIXZZ','ZXIXZ','ZZZZZ']
# steane-code in the |0> logical state
#S = ['XXXXIII','IXXIXXI','IIXXIXX','ZZZZIII','IZZIZZI','IIZZIZZ','ZZZZZZZ']

def convert(S , c=[] , t=[] , shuffle=False):
    # number of qubits N and number of stabilizers Ns
    N = len(S[0])
    Ns = len(S)
    # binary representation
    A = np.array([[S[j][i] in {'Y','Z'} for j in range(N)] for i in range(Ns)]+
                 [[S[j][i] in {'Y','X'} for j in range(N)] for i in range(Ns)],
                 dtype=int)
    # raise Exception if there are not enough stabilizers
    if Ns < N:
        raise Exception('The number of stabilizers in S can not be smaller\
                         than the number of qubits')
    # raise Exception if rank(p) is not N
    if gb.rank(gb.gauss(A))!=N:
        raise Exception('S must contain the same number of independent\
                         stabilizers than qubits')
    # raise Exception if stabilizers do not commute
    for i in range(Ns-1):
        for j in range(i+1,Ns):
            if A[:N,i].dot(A[N:,j])%2 !=  A[:N,j].dot(A[N:,i])%2:
                raise Exception('generators '+S[i]+' and '+S[j]+' do not\
                                 commute')
    # raise Exception if control and target qubits have non empty intersection
    if set(c).intersection(set(t)):
        raise Exception('c and t must have empty intersection')
    # if a control or target qubit has a label bigger than N
    for i in c:
        if i>=N:
            raise Exception('control qubits must be labelled from 0 to N-1')
    for i in t:
        if i>=N:
            raise Exception('target qubits must be labelled from 0 to N-1')
    # shuffle rows
    if shuffle:
        qubits = c + random.sample(set(range(N))-(set(c).union(set(t))),\
                                   N-len(c)-len(t)) + t
    else:
        qubits = c + list(set(range(N))-(set(c).union(set(t)))) + t
    # reorder rows in A
    A = np.array([A[qubits[i]] for i in range(N)] +\
                 [A[qubits[i]+N] for i in range(N)])                
    # put A in the form that allows to make Gauss elimination in the right way
    # add the identity to monitor the recombinations performed by the Gaussian 
    #elimination
    A = np.array(list(A[N:,:])+list(A[:N,:])+list(np.eye(N,dtype=int)))
    A = A.T
    # perform the Gaussian elimination and identify the matrix R that performs
    # that recombination
    A = gb.gauss(A)
    R = A[:,-N:].T
    A = A[:,:-N]
    n = gb.rank(A[:,:N])    
    if len(c) > n:
        raise Exception('too many control qubits selected')
    if len(t) > N - n:
        raise Exception('too many target qubits selected')
    if gb.rank(A[:n,:len(c)]) < len(c):
        raise Exception('wrong selection of control qubits')    
    # select control and target qubits 
    for i in range(len(c),n):
        for j in range(i,N):
            if A[i,j] and qubits[j] not in t and qubits[j] not in c:
                c.append(qubits[j])
                break
            elif qubits[j] not in t:
                t.append(qubits[j])                
    t = t + list(set(range(N))-set(c).union(set(t)))
    if len(c)!=n:
        raise Exception('wrong selection of control and/or target qubits')    
    # put A back to the original form
    A = A.T
    A = np.array(list(A[N:,:])+list(A[:N,:]))
    A = np.array([A[qubits.index(i)] for i in c] +\
                 [A[qubits.index(i)] for i in t]+
                 [A[qubits.index(i)+N] for i in c] +\
                 [A[qubits.index(i)+N] for i in t])
    # update the order of the qubits in the rows
    qubits = c+t
    # build Xlc and the inverse
    Xlc = A[N:n+N,:n]
    if n!=0:
        Xlc_inv = gb.inverse(Xlc)
    else:
        Xlc_inv = np.zeros((0,0),dtype=int)
    # identify the rest of blocks
    Xlt = A[n+N:,:n]
    Zlc = A[:n,:n]
    Zlt = A[n:N,:n]
    Zrt = A[n:N,n:]    
    # compute R
    if n!=N:
        Zrt_inv = gb.inverse(Zrt)
    else:
        Zrt_inv = np.zeros((0,0),dtype=int)
    R = R.dot(np.block([[Xlc_inv,np.zeros((n,N-n),dtype=int)],
                        [Zrt_inv.dot(Zlt).dot(Xlc_inv)%2,Zrt_inv]]))%2
    # compute C and B and obtain the list z of the qubits where the z-rotation
    # is performed
    B = Xlt.dot(Xlc_inv)%2
    C = (Zlc+B.T.dot(Zlt)).dot(Xlc_inv)%2
    z = [qubits[i] for i in range(n) if C[i,i]]
    C = (C + np.diag(np.diagonal(C)))%2    
    # Adjacency matrix
    G = np.block([[C,B.T],[B,np.zeros((N-n,N-n),dtype=int)]])
    G = np.array([G[qubits.index(i)] for i in range(N)])
    G = np.array([G[:,qubits.index(i)] for i in range(N)]).T
    return G , sorted(c) , sorted(t) , sorted(z) , R    