#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 15:33:58 2019

@author: davamaro
"""

import numpy as np
from copy import deepcopy

##### GAUSSIAN ELIMINATION ON BINARY MATRICES WITH SUM MODULO 2 ###############

# gauss elimination of rows
def gauss(A):
    output=deepcopy(A)
    (m,n) = np.shape(A)
    if m == 0 or n==0:
        raise NameError('empty matrix')
        output = np.zeros((m,n))
    row=-1
    for j in range(n):
        pivot=-1
        for i in range(row+1,m):
            if output[i][j]:
                if pivot==-1:
                    pivot=i
                    row+=1
                    output[[row,pivot]]=output[[pivot,row]]
                else:
                    output[i]=np.logical_xor(output[i],output[row])
    return output

# rank of a matrix A. A must be in gauss form
def rank(A):
    m , n =np.shape(A)
    rank=0
    for i in range(m):
        for j in range(i,n):
            if A[i][j]:
                rank+=1
                break
    return rank
  
# inverse of a matrix
def inverse(A):
    (m,n)=np.shape(A)
    if m!=n:
        raise NameError('a not square matrix does not have inverse')
        output=np.empty((m,0),bool)
    A = np.block([[A , np.eye(m,dtype=bool)]])
    A = gauss(A)
    if rank(A[:,:m]) < m:
        raise Exception('singular matrix')
    A = np.array([A[-i-1] for i in range(m)]).T
    A = np.array([A[m-i-1] for i in range(m)] + [A[-i-1] for i in range(m)]).T
    output = gauss(A)[:,m:]
    output =np.array([[output[-i-1,-j-1] for j in range(m)] for i in range(m)])
    return output