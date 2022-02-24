#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 17:26:13 2022

@author: alexander
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
mod = xfeat.Model(mat, size=500)
# create atomic core
mod.atoms()
# create mesh and set up system stiffness matrix
mod.mesh()
# create screw dislocation
mod.init_dislo()
# plot nodes with boundary conditions
mod.plot_nodal('ubc')

# iterate into relaxed configuration
for i in range(5):
    # relax atomic region under shear strain
    mod.relax_atoms(i)  # relax atomic structure with fixed boundary atoms
    mod.atom_bc()  # apply relaxed atom positions as BC to XFEM 
    mod.solve()  # colculate nodal displacements for mechanical equilibrium
    mod.shift_atoms()  # move boundary atoms according to strain field

# plot nodes with boundary conditions
mod.plot_nodal('ubc')
#plot nodal displacement
mod.plot_nodal('u')
#plot nodal force
mod.plot_nodal('f')
# plot stresses
sig0 = mod.calc_stress()
mod.plot_el('sig', sig=sig0)

# Apply sub-critical shear stress on boundary
mod.apply_bc(1.55)
# iterate into relaxed configuration
for i in range(10):
    mod.relax_atoms(i, name='applied_stress_155')
    mod.atom_bc()
    mod.solve()
    mod.shift_atoms()

# plot stresses
sig1 = mod.calc_stress()
mod.plot_el('sig', sig=sig1)

# Apply critical shear stress on boundary
mod.apply_bc(1.65)
# iterate into relaxed configuration
for i in range(10):
    mod.relax_atoms(i, name='applied_stress_165')
    mod.atom_bc()
    mod.solve()
    mod.shift_atoms()

# plot stresses
sig2 = mod.calc_stress()
mod.plot_el('sig', sig=sig2)
