# distutils: language = c++
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test module
Created on Sat Feb  5 17:26:13 2022

@author: alexander
"""
import numpy as np
import pyvista as pv

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

def pv_mesh(nodes, elmts):
    faces = set()
    lines = set()
    for vox in elmts:
        vert = nodes[vox]
        x0 = np.amin(vert[:, 0])
        x1 = np.amax(vert[:, 0])
        y0 = np.amin(vert[:, 1])
        y1 = np.amax(vert[:, 1])
        z0 = np.amin(vert[:, 2])
        z1 = np.amax(vert[:, 2])
        ind = np.nonzero(np.abs(vert[:, 0] - x0) < 1.e-4)[0]
        assert(len(ind) == 4)
        key = '{0}_{1}_{2}_{3}'.format(*[x for x in sorted(vox[ind])])
        faces.add(key)
        #find 4 lines
        il = np.nonzero(np.abs(vert[ind, 1] - y0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 1] - y1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 2] - z0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 2] - z1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        
        # find face
        ind = np.nonzero(np.abs(vert[:, 0] - x1) < 1.e-4)[0]
        assert(len(ind) == 4)
        key = '{0}_{1}_{2}_{3}'.format(*[x for x in sorted(vox[ind])])
        faces.add(key)
        #find lines
        il = np.nonzero(np.abs(vert[ind, 1] - y0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 1] - y1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 2] - z0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 2] - z1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        
        ind = np.nonzero(np.abs(vert[:, 1] - y0) < 1.e-4)[0]
        assert(len(ind) == 4)
        key = '{0}_{1}_{2}_{3}'.format(*[x for x in sorted(vox[ind])])
        faces.add(key)
        #find lines
        il = np.nonzero(np.abs(vert[ind, 0] - x0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 0] - x1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 2] - z0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 2] - z1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        
        ind = np.nonzero(np.abs(vert[:, 1] - y1) < 1.e-4)[0]
        assert(len(ind) == 4)
        key = '{0}_{1}_{2}_{3}'.format(*[x for x in sorted(vox[ind])])
        faces.add(key)
        #find lines
        il = np.nonzero(np.abs(vert[ind, 0] - x0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 0] - x1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 2] - z0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 2] - z1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        
        ind = np.nonzero(np.abs(vert[:, 2] - z0) < 1.e-4)[0]
        assert(len(ind) == 4)
        key = '{0}_{1}_{2}_{3}'.format(*[x for x in sorted(vox[ind])])
        faces.add(key)
        #find lines
        il = np.nonzero(np.abs(vert[ind, 1] - y0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 1] - y1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 0] - x0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 0] - x1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        
        ind = np.nonzero(np.abs(vert[:, 2] - z1) < 1.e-4)[0]
        assert(len(ind) == 4)
        key = '{0}_{1}_{2}_{3}'.format(*[x for x in sorted(vox[ind])])
        faces.add(key)
        #find lines
        il = np.nonzero(np.abs(vert[ind, 1] - y0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 1] - y1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 0] - x0) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)
        il = np.nonzero(np.abs(vert[ind, 0] - x1) < 1.e-4)[0]
        assert(len(il) == 2)
        key = '{0}_{1}'.format(*[x for x in sorted(vox[ind[il]])])
        lines.add(key)

    nf = len(faces)*2
    flist = []
    for key in faces:
        hh = key.split('_')
        flist.append(3)
        [flist.append(int(x)) for x in hh[0:3]]
        flist.append(3)
        [flist.append(int(x)) for x in hh[1:4]]


    nl = len(lines)
    llist = []
    for key in lines:
        llist.append(2)
        [llist.append(int(x)) for x in key.split('_')]

    mesh = pv.PolyData(var_inp=nodes, faces=flist, n_faces=nf,
                       lines=llist, n_lines=nl)
    return mesh