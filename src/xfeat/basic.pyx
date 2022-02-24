# distutils: language = c++
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test module
Created on Sat Feb  5 17:26:13 2022

@author: alexander
"""
import numpy as np

def primes(int nb_primes):
    cdef int n, i, len_p
    cdef int p[1000]

    if nb_primes > 1000:
        nb_primes = 1000




    len_p = 0  # The current number of elements in p.
    n = 2
    while len_p < nb_primes:
        # Is n prime?
        for i in p[:len_p]:
            if n % i == 0:
                break

        # If no break occurred in the loop, we have a prime.
        else:
            p[len_p] = n
            len_p += 1
        n += 1

    # Let's copy the result into a Python list:
    result_as_list = [prime for prime in p[:len_p]]
    return result_as_list

def rot_elast_tens(double C11, double C12, double C44, x):
    """
    Calculate elastic tensor rotated into orientations given by Miller indices.

    Parameters
    ----------
    C11, C12, C44 : float
        Elastic constants
    x : (3, 3) array
        Miller indeces of crystal orientation

    Returns
    -------
    EE : (6, 6) array
        Rotated elasticity tensor

    """
    cdef int i, j, k, l, p, q, s, t
    cdef double[3] norm
    cdef double[3][3] R
    x = np.array(x, dtype=np.double)
    
    Dt = np.array([
        [C11, C12, C12, 0.,    0.,  0.], \
        [C12, C11, C12, 0.,    0.,  0.], \
        [C12, C12, C11, 0.,    0.,  0.], \
        [0.,  0.,  0.,  2*C44, 0.,  0.], \
        [0.,  0.,  0.,  0., 2*C44,  0.], \
        [0.,  0.,  0.,  0.,    0., 2*C44]], dtype=np.double)

    norm = np.linalg.norm(x, axis=1)
    for i in range(3):
        for j in range(3):
            R[i][j] = x[i][j]/norm[i]

    TransD6toD9 = np.zeros((9, 6), dtype=np.double)
    TransD6toD9[0, 0] = TransD6toD9[4, 1] = TransD6toD9[8, 2] = 1.0
    TransD6toD9[1, 3] = TransD6toD9[3, 3] = TransD6toD9[2, 4] =\
    TransD6toD9[6, 4] = TransD6toD9[5, 5] = TransD6toD9[7, 5] = np.sqrt(0.5)

    Dmat = np.zeros((9,9), dtype=np.double)
    for i in range(9):
        for j in range(9):
            for k in range(6):
                for l in range(6):
                    Dmat[i, j] += TransD6toD9[i, k] * Dt[k, l] * TransD6toD9[j, l]

    QDmat = np.zeros((9, 9), dtype=np.double)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for p in range(3):
                        for q in range(3):
                            for s in range(3):
                                for t in range(3):
                                    QDmat[j + 3*i][l + 3*k] +=\
    									 R[i][p] * R[j][q] * R[k][s] * R[l][t] *\
    									 Dmat[q + 3*p][t + 3*s]
                       
    # Rotated Elasticity Tensor
    EE = np.zeros((6, 6), dtype=np.double)
    EE[0, 0] = QDmat[0, 0]
    EE[0, 1] = QDmat[0, 4]
    EE[0, 2] = QDmat[0, 8]
    EE[0, 3] = QDmat[0, 1]
    EE[0, 4] = QDmat[0, 5]
    EE[0, 5] = QDmat[0, 6]

    EE[1, 1] = QDmat[4, 4]
    EE[1, 2] = QDmat[4, 8]
    EE[1, 3] = QDmat[4, 1]
    EE[1, 4] = QDmat[4, 5]
    EE[1, 5] = QDmat[4, 6]

    EE[2, 2] = QDmat[8, 8]
    EE[2, 3] = QDmat[8, 1]
    EE[2, 4] = QDmat[8, 5]
    EE[2, 5] = QDmat[8, 6]

    EE[3, 3] = QDmat[1, 1]
    EE[3, 4] = QDmat[1, 5]
    EE[3, 5] = QDmat[1, 6]

    EE[4, 4] = QDmat[5, 5]
    EE[4, 5] = QDmat[5, 6]

    EE[5, 5] = QDmat[6, 6]   
    for i in range(6):
        for j in range(0, i):
            EE[i, j] = EE[j, i]
    return EE
