#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
XFEAT 
Create a screw dislocation in a bcc FE crystal and apply a shear stress until the
dislocation becomes mobile.

Author: Alexander Hartmaier
Institution: ICAMS / Ruhr University Bochum, Germany
Date: March 2022

Published under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
"""

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
mod = xfeat.Model(mat, size=500, verbose=True)
# create atomic core
mod.atoms([10, 17, 3])
# create mesh and set up system stiffness matrix
mod.mesh()
# create screw dislocation
mod.init_dislo([0, 0, 1])
# plot nodes with boundary conditions
#mod.plot('ubcz')

# iterate into relaxed configuration
for i in range(10):
    # relax atomic region under shear strain
    mod.atom_bc()  # apply relaxed atom positions as BC to XFEM 
    mod.solve()  # colculate nodal displacements for mechanical equilibrium
    mod.shift_atoms()  # move boundary atoms according to strain field
    mod.relax_atoms(i)  # relax atomic structure with fixed boundary atoms
    mod.plot('uy')

# plot nodes with boundary conditions
mod.plot('ubcz')
# plot nodal displacement
mod.plot('uz')
# plot potential energy of atoms
mod.plot('epot')
# plot stresses
mod.plot('sigyz')

# Apply sub-critical shear stress on boundary
mod.apply_bc(1.55)
# iterate into relaxed configuration
for i in range(10):
    mod.atom_bc()
    mod.solve()
    mod.shift_atoms()
    mod.relax_atoms(i, name='applied_stress_155')

# plot stresses
mod.plot('sigyz')
mod.plot('epot')

# Apply critical shear stress on boundary
mod.apply_bc(1.65)
# iterate into relaxed configuration
for i in range(10):
    mod.atom_bc()
    mod.solve()
    mod.shift_atoms()
    mod.relax_atoms(i, name='applied_stress_165')

# plot stresses
mod.plot('sigyz')
mod.plot('epot')

