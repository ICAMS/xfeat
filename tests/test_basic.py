#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 17:26:13 2022

@author: alexander
"""
import os
import numpy as np
import xfeat

def test_rot():
    ref = np.array([
        [1.93618, 0.766065, 0.627083, 0., 0., -0.19655],
        [0.766065, 1.93618, 0.627083, 0., 0., 0.19655], 
        [0.627083, 0.627083, 2.07516 , 0, 0, 0],
        [0, 0, 0, 0.585055, 0.19655, 0],
        [0, 0, 0, 0.19655, 0.446073, 0],
        [-0.19655, 0.19655, 0, 0, 0, 0.446073]])
    assert(np.all(np.abs(mod.Cel - ref) < 1.e-5))
    
def test_atoms():
    hh = mod.shift-np.array([-0.0003508092037449728, -0.004243414970321879])
    assert(np.abs(np.amax(hh)) < 1.e-5)
    assert(np.linalg.norm(mod.dist - 
                          np.array([6.994057, 4.038021, 2.472772])) < 1.e-5)
    assert(os.path.isfile('{}/relaxed_perfect_crystal_with_atom_type.imd'
                          .format(temp_dir)))
    
def test_mesh():
    assert(mod.xdof == 20490)
    assert(np.amax(mod.a_csr.diagonal()) - 18.516398876358256 < 1.e-5)
    
def test_init_dislo():
    assert(np.abs(mod.u[1440] + 0.0078246745716456) < 1.e-5)
    assert(np.abs(mod.u[1441] + 0.006752027037675507) < 1.e-5)
    assert(np.abs(mod.u[7440] - 0.0003205111755714418) < 1.e-5)
    assert(os.path.isfile('{}/atomistic_dislocation_with_fem_solution.imd'
                          .format(temp_dir)))

def test_iteration():
    # relax atomic region under shear strain
    mod.atom_bc()  # apply relaxed atom positions as BC to XFEM 
    mod.solve()  # colculate nodal displacements for mechanical equilibrium
    
    # calculate stresses
    sig0 = mod.calc_stress()
    assert(np.abs(np.amax(sig0) - 2.2458018354577853) < 1.e-5)
    assert(np.abs(np.amin(sig0) + 2.1907796579448022) < 1.e-5)

def test_apply_bc():
    # Apply shear stress on boundary
    mod.apply_bc(1.55)
    assert(np.abs(np.amax(mod.ubc) - 2.3061843499574897) < 1.e-5)
    assert(np.abs(np.amin(mod.ubc) + 2.306184349957489) < 1.e-5)
    assert(np.abs(mod.u[1440] + 0.012036538734587841) < 1.e-5)
    assert(np.abs(mod.u[1439] + 1.623329283079101) < 1.e-5)
    assert(os.path.isfile('{}/atomistic_dislocation_with_fem_solution.imd'
                          .format(temp_dir)))

    #remove temp directory
    os.system('rm -rf {}'.format(temp_dir))
    

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

cwd = os.getcwd()
temp_dir = cwd + '/temp'
if os.path.exists(temp_dir):
        os.system('rm -rf {}'.format(temp_dir))
os.makedirs(temp_dir)
#os.chdir(cwd[:-6])
#print('CWD: ', os.getcwd())
#print(os.path.isfile('JMakeConfig.jar'))
# create XFEM model 
mod = xfeat.Model(mat, size=200)
# create atomic core
mod.atoms()
# create mesh and set up system stiffness matrix
mod.mesh()
# create screw dislocation
mod.init_dislo()
