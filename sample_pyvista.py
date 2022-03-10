#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 17:26:13 2022

@author: alexander
"""
import numpy as np
import pyvista as pv
from vtk import VTK_HEXAHEDRON
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
mod.atoms()
mod.mesh()
mod.init_dislo()
#mod.plot('ubcz')
# iterate into equilibrium configuration
for i in range(1):
    mod.atom_bc()  # apply relaxed atom positions as BC to XFEM 
    mod.solve()  # colculate nodal displacements for mechanical equilibrium
    mod.shift_atoms()  # move boundary atoms according to strain field
    mod.relax_atoms(i)  # relax atomic structure with fixed boundary atoms
    #mod.plot('ux')

#mod.plot('sigyz')
#mod.plot('uz')
#mod.plot('dispz')
#mod.plot('epot')

