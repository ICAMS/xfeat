#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 17:26:13 2022

@author: alexander
"""
import numpy as np
import matplotlib.pyplot as plt
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
       'ori_x' : [ 1, 1,  1],
       'ori_y' : [-1, 0,  1],
       'ori_z' : [-1, 2, -1],
       }

# create XFEM model 
mod = xfeat.Model(mat, size=200)
mod.atoms([15, 17, 2])
mod.mesh()
#mod.grid.plot(show_edges=True)
mod.init_dislo([1, 0, 0])
#mod.plot('ubcx')
#mod.plot('ubcy')
# iterate into equilibrium configuration
for i in range(10):
    mod.atom_bc()  # apply relaxed atom positions as BC to XFEM 
    mod.solve()  # colculate nodal displacements for mechanical equilibrium
    mod.shift_atoms()  # move boundary atoms according to strain field
    mod.relax_atoms(i)  # relax atomic structure with fixed boundary atoms
    #mod.plot('uy')

plt.scatter(mod.apos[mod.a2n[:,1], 0], mod.apos[mod.a2n[:,1], 1])
plt.show()
plt.scatter(mod.nodes[mod.a2n[:,0], 0], mod.nodes[mod.a2n[:,0], 1])
plt.show()
dvec = mod.nodes[mod.a2n[:,0], :] - mod.apos[mod.a2n[:,1], :]
dist = np.linalg.norm(dvec, axis=1)
print('Minimum distance b/w atom and node: {}, maximum distance: {}'
      .format(np.amin(dist), np.amax(dist)))
mod.plot('sigxy')
mod.plot('ux')
mod.plot('dispx')
mod.plot('epot')

# Apply shear stress on boundary
mod.apply_bc(0.55, comp='xy')
# iterate into relaxed configuration
for i in range(10):
    mod.atom_bc()
    mod.solve()
    mod.shift_atoms()
    mod.relax_atoms(i, name='applied_stress_055')


