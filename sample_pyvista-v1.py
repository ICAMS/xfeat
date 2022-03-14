#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 17:26:13 2022

@author: alexander
"""
import numpy as np
import pyvista as pv
import xfeat


# define material of single crystal as dictionary
mat = {
       'name' : 'iron',
       'cs'   : 'bcc',  # crystal structure
       'lp'   : 2.8553,  # lattice parameter in Angstrom
       'mass' : 55.845,  # atomic mass in u
       # define anisotropic elastic constants in GPa
       'C11'  : 243.4,
       'C12'  : 145.0,
       'C44'  : 116.0,
       # W-Bop 6
       # C11. C12, C44 = 2.837825, 1.5317, 0.6256225
       # define crystallograhic orientation of crystal
       'ori_x' : [-1, 2, -1],
       'ori_y' : [-1, 0,  1],
       'ori_z' : [ 1, 1,  1],
       }

# create XFEM model 
mod = xfeat.Model(mat, size=200)
# create atomic core
mod.atoms()
# create mesh and set up system stiffness matrix
mod.mesh()

faces = []
lines = []
for vox in mod.elmts:
    vert = mod.nodes[vox]
    x0 = np.amin(vert[:, 0])
    x1 = np.amax(vert[:, 0])
    y0 = np.amin(vert[:, 1])
    y1 = np.amax(vert[:, 1])
    z0 = np.amin(vert[:, 2])
    z1 = np.amax(vert[:, 2])
    ind = np.nonzero(np.abs(vert[:, 0] - x0) < 1.e-4)[0]
    faces.append(4)
    assert(len(ind) == 4)
    [faces.append(vox[x]) for x in ind]
    #find lines
    il = np.nonzero(np.abs(vert[ind, 1] - y0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 1] - y1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 2] - z0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 2] - z1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    
    # find face
    ind = np.nonzero(np.abs(vert[:, 0] - x1) < 1.e-4)[0]
    faces.append(4)
    assert(len(ind) == 4)
    [faces.append(vox[x]) for x in ind]
    #find lines
    il = np.nonzero(np.abs(vert[ind, 1] - y0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 1] - y1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 2] - z0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 2] - z1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    
    ind = np.nonzero(np.abs(vert[:, 1] - y0) < 1.e-4)[0]
    faces.append(4)
    assert(len(ind) == 4)
    [faces.append(vox[x]) for x in ind]
    #find lines
    il = np.nonzero(np.abs(vert[ind, 0] - x0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 0] - x1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 2] - z0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 2] - z1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    
    ind = np.nonzero(np.abs(vert[:, 1] - y1) < 1.e-4)[0]
    faces.append(4)
    assert(len(ind) == 4)
    [faces.append(vox[x]) for x in ind]
    #find lines
    il = np.nonzero(np.abs(vert[ind, 0] - x0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 0] - x1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 2] - z0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 2] - z1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    
    ind = np.nonzero(np.abs(vert[:, 2] - z0) < 1.e-4)[0]
    faces.append(4)
    assert(len(ind) == 4)
    [faces.append(vox[x]) for x in ind]
    #find lines
    il = np.nonzero(np.abs(vert[ind, 1] - y0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 1] - y1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 0] - x0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 0] - x1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    
    ind = np.nonzero(np.abs(vert[:, 2] - z1) < 1.e-4)[0]
    faces.append(4)
    assert(len(ind) == 4)
    [faces.append(vox[x]) for x in ind]
    #find lines
    il = np.nonzero(np.abs(vert[ind, 1] - y0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 1] - y1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 0] - x0) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    il = np.nonzero(np.abs(vert[ind, 0] - x1) < 1.e-4)[0]
    assert(len(il) == 2)
    lines.append(2)
    lines.append(vox[ind[il[0]]])
    lines.append(vox[ind[il[1]]])
    
nf = int(len(faces)/5)
nl = int(len(lines)/3)
mesh = pv.PolyData(var_inp=mod.nodes, faces=faces, n_faces=nf,
                   lines=lines, n_lines=nl)

pl = pv.Plotter()
pl.add_mesh(mesh, show_edges=True, color='white')
pl.add_points(mesh.points, color='red', point_size=10)

single_cell = mesh.extract_cells(1)
pl.add_mesh(single_cell, color='pink', edge_color='blue',
            line_width=5, show_edges=True)

pl.camera_position = [(6.20, 3.00, 7.50),
                      (0.16, 0.13, 2.65),
                      (-0.28, 0.94, -0.21)]
pl.show()