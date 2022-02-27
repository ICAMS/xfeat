# distutils: language = c++
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 17:46:17 2022

@author: alexander
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import vtk
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from xfeat cimport cpp_wrapper as xfc
from xfeat.basic import rot_elast_tens

class Model(object):
    def __init__(self, mat, size=500):
        self.fem_size = size  # size of FEM part
        self.bv = None  # Burgers vector
        self.nodes = None  # nodal positions
        self.Nnode = None  # Number of nodes
        self.elmts = None  # elements
        self.r_in = None  # inner coordinate of fem domain
        self.r_out = None  # outer coordinate of fem domain
        self.b_vec = None  # Burgers vector
        self.nd_in = None  # inner nodes
        self.nd_out = None  # outer nodes
        self.nd_dis = None  # nodal pairs along dislocation line
        self.r_dis = None  # nodal pairs coordinates along dislocation line
        self.ubc = None  # displacement boundary conditions
        self.a_csr = None  # system stiffness matrix in CSR format
        self.F_vec = None  # nodal forces in CSR format
        self.JJglob = None  # Pointer to equation number from full nodal DOF
        self.u = None  # nodal displacements in dense notation
        self.xdof = None  # number of xtended degrees of freedom
        
        self.mat = mat
        self. lp = self.mat['lp']  # lattice parameter
        ori_x = mat['ori_x']
        ori_y = mat['ori_y']
        ori_z = mat['ori_z']
        # define orientation of crystal axes by Miller indices
        self.ori = np.array([ori_x, ori_y, ori_z], dtype=int)
        self.C11 = mat['C11']
        self.C12 = mat['C12']
        self.C44 = mat['C44']
        self.Cel = rot_elast_tens(self.C11, self.C12, self.C44, self.ori)\
                   * 0.0062417  # elastic tensor in eV/A^3

        # create temp directory for atomistic files
        cwd = os.getcwd()
        path = ''
        pel = cwd.split('/')
        for hs in pel[0:pel.index('XFEAt')+1]:
            path += hs+'/'
        self.main = path
        self.temp = path + '/temp'
        self.md_dir = path + 'Fe_MD'
        self.libs = path + 'libs'
        if not os.path.exists(self.temp):
            os.makedirs(self.temp)
        if not os.path.exists(self.md_dir):
            raise RuntimeError('Directory "{}" is not existing.'
                               .format(self.md_dir))
        if not os.path.isfile(self.md_dir+'/imd_eam_fire_homdef_stress_nbl'):
            raise RuntimeError('IMD executable "imd_eam_fire_homdef_stress_nbl" is not existing in path {}.'
                               .format(self.md_dir))
        tdir = self.temp.encode()
        xfc.temp_dir = tdir
        xfc.fem_size = self.fem_size
        xfc.PI = np.pi
        
    def atoms(self):
        '''
        Set up atomic core. 
        Produces IMD input file with atom types for boundary conditions.

        Returns
        -------
        res : TYPE
            DESCRIPTION.
            
        Attributes
        ----------
        bv : float
            Norm of Brugers vector after relaxation
        lp : float
            Lattice parameter after relaxtion
        shift : (2,) array
            Shift between FEM and atomistic origins

        '''
        # create relaxed perfect crystal as IMD input (Steps 2 and 3)
        name_urel = self.temp + '/unrelaxed_perfect_crystal.imd'
        sx = 10 * self.lp * np.linalg.norm(self.ori[0,:])
        sy = 17 * self.lp * np.linalg.norm(self.ori[1,:])
        sz =  3 * self.lp * np.linalg.norm(self.ori[2,:])
        with open(self.temp+'/imd_create_cryst.param', 'w') as fd:
            fd.write('structure    {}\n'.format(self.mat['cs']))
            fd.write('new_z    {} {} {}\n'\
                     .format(self.ori[2, 0], self.ori[2, 1], self.ori[2, 2]))
            fd.write('new_y    {} {} {}\n'\
                     .format(self.ori[1, 0], self.ori[1, 1], self.ori[1, 2]))
            fd.write('outfile    {}\n'.format(name_urel))
            fd.write('startnr    1\n')
            fd.write('lattice_const    {}\n'.format(self.lp))
            fd.write('mass    {}\n'.format(self.mat['mass']))
            fd.write('pbc    1 1 1\n')
            fd.write('box_x    {}\n'.format(sx))
            fd.write('box_y    {}\n'.format(sy))
            fd.write('box_z    {}\n'.format(sz))
        
        # create crystal structure with java tool (Step 2)
        os.system(' cd {}'.format(self.main))
        os.system('java -jar {0}/JMakeConfig.jar {1}'\
                  .format(self.libs, self.temp+'/imd_create_cryst.param'))
        os.system('mv crystal.conf {}'.format(self.temp))
        # relax crystal with fire algorithm in IMD (Step 3)
        os.system('cd {}; ./imd_eam_fire_homdef_stress_nbl -p Fe101-glok-sample.param'
                  .format(self.md_dir))
        os.system('mv {0}/relaxed_perfect_crystal.imd.00000.ss {0}/relaxed_perfect_crystal.imd'\
                  .format(self.temp))
        os.system('rm {0}/*eng {0}/*itr {0}/*ssdef {0}/*chkpt; cd ..'\
                  .format(self.temp))
        #xfc.get_n_atoms()  # get number of atoms and allocate memory
        xfc.atom_set_up()  # define types for atoms on which BC are applied
        self.bv = xfc.bv  # Burgers vector
        if np.abs(2.*self.bv/np.sqrt(3.) - self.lp) > 1.e-9:
            print('Correcting lattice parameter in material definition after atomic relaxation.')
            print('Given value: {}'.format(self.lp))
            self.lp = 2.*self.bv/np.sqrt(3.)
            print('New value: {}'.format(self.lp))
        self.shift = xfc.shift  # shift along slip plane
        self.Lx = xfc.Lx  # dimensions of atomic box
        self.Ly = xfc.Ly
        self.Lz = xfc.Lz
        self.natom = xfc.natom  # number of atoms
        # get reference atomic positions in coordinates of XFEM model
        hh = np.array(xfc.coords, dtype=np.double)
        self.apos = hh[1:self.natom+1, 0:3]
        self.apos[:, 0] -= 0.5*self.Lx - self.shift[0]
        self.apos[:, 1] -= 0.5*self.Ly - self.shift[1]
        #self.apos[:, 2] -= 0.5*self.Lz
        return
    
    
    def mesh(self):
        '''
        Create mesh of XFEAt model and setup system stiffness matrix
        in CSR form

        Returns
        -------
        None.
        
        Attributes
        ----------
        Nnode
        nodes
        NEL
        elmts
        el_ctr
        NEL_x, NEL_y, NEL_z
        LX, Ly, LZ
        JJglob
        xdof
        NLAG
        NNODE
        NDF
        Fglob
        a_csr
        f
        ubc

        '''
        cdef int i, ii, iicon, cf
        
        if self.bv is None:
            raise AttributeError('No atomic box found. Run "get_atoms" first.')
        # setup elastic tensor for XFEM
        EE3D = np.zeros((7, 7), dtype=np.double)
        EE3D[1:7, 1:7] = self.Cel
        xfc.EE3D = EE3D
        
        # cerate XFEM mesh
        xfc.create_mesh()
        self.Nnode = xfc.NNODE  # number of nodes
        self.nodes = np.zeros((self.Nnode, 3), dtype=np.double)
        assert(len(xfc.Xglob) == self.Nnode)
        self.nodes[:, 0] = xfc.Xglob
        self.nodes[:, 1] = xfc.Yglob
        self.nodes[:, 2] = xfc.Zglob
        #self.n_cnstr = xfc.NCONT
        #self.nodes_c = np.array(xfc.ContraintNodes, dtype=int)
        LOTOGO = np.array(xfc.LOTOGO, dtype=int) - 1
        self.NEL = int(len(LOTOGO) / 8)
        self.elmts = LOTOGO.reshape((self.NEL, 8))
        self.el_ctr = np.zeros((self.NEL, 3), dtype=np.double)
        for IEL in range(self.NEL):
            self.el_ctr[IEL,0] = np.mean(self.nodes[self.elmts[IEL,:], 0])
            self.el_ctr[IEL,1] = np.mean(self.nodes[self.elmts[IEL,:], 1])
            self.el_ctr[IEL,2] = np.mean(self.nodes[self.elmts[IEL,:], 2])
        self.NEL_x = xfc.nelem_full_big_box_x
        self.NEL_y = xfc.nelem_full_big_box_y
        self.NEL_z = xfc.nelem_full_big_box_z
        
        # setup system stiffness matrix for XFEM
        xfc.init_matrix()  # calculates Aglob matrix in CSR format
        self.JJglob = np.zeros(3*self.Nnode, dtype=int)
        self.JJglob[:]  = xfc.JJglob[3:]
        self.xdof = xfc.xdof
        assert(self.NEL == xfc.NEL)
        self.NLAG = xfc.NLAG  # number of augmented elmts (Lagrangians)
        self.NDF = xfc.NDF
        a_csr = csr_matrix((xfc.NEWAglob, xfc.Jglob, xfc.Iglob))
        assert(a_csr.shape == (self.xdof, self.xdof))
        # scipy.sparse does not support symmetric matrices
        # full matrix is constructed by adding lower triangle
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
        a_csr.setdiag(0.5*a_csr.diagonal())
        self.a_csr = a_csr + a_csr.transpose()
        
        # create pyVista mesh
        cells = np.ones((self.NEL, 9), dtype=int)*8
        cells[:,1:9] = self.elmts
        celltypes = np.empty(self.NEL, dtype=np.uint8)
        celltypes[:] = vtk.VTK_HEXAHEDRON
        self.grid = pv.UnstructuredGrid(cells.ravel(), celltypes, self.nodes)


    def solve(self):
        '''
        Solve system of equations A.disp = F_vec for unknown displacements 

        Returns
        -------
        None.
        
        Attributes
        ----------
        u

        '''
        xfc.Fglob = np.zeros(self.NDF, dtype=np.double)
        xfc.update_fglob()  # update Fglob
        self.Fglob = np.array(xfc.Fglob, dtype=np.double)
        F_vec = np.zeros(self.xdof, dtype=np.double)
        F_vec[0:self.NDF] = self.Fglob

        print('\n********************\nstart solver\n')
        t0 = time.time()
        disp = spsolve(self.a_csr, F_vec)
        t1 = time.time()
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve.html#scipy.sparse.linalg.spsolve
        assert(np.allclose(self.a_csr.dot(disp), F_vec))
        print('Solution successful in {:6.5}s.'.format(t1-t0))
        xfc.Dglob = disp
        
        # store solution with nodal indices
        self.u = np.zeros(3*self.Nnode, dtype=np.double)
        for i in range(3*self.Nnode):
            iicon = self.JJglob[i]
            if iicon >= 0:
                self.u[i] = disp[iicon]
        # add boundary conditions
        self.ubc = np.array(xfc.ENFRDISPglob[3:], dtype=np.double)
        self.u  += self.ubc
        xfc.node_dis = self.u
        
        # compose force vector
        self.f = np.zeros(3*self.Nnode)
        for i in range(3*self.Nnode):
            iicon = self.JJglob[i]
            if iicon >= 0: 
                self.f[i] = self.Fglob[iicon]
        return
    
    def init_dislo(self):
        '''
        Introduce screw dislocation and write updated atomic configuration

        Returns
        -------
        None.

        '''
        # create Volterra screw dislocation by introducing a shift in
        # XFEM elements and atomic core
        xfc.create_xfem_dis()  # create displacements in XFEM region
        self.ubc = np.array(xfc.ENFRDISPglob[3:], dtype=np.double)
        self.solve()  # calculate nodal displacements under BC of shift
        xfc.atom_element()  # assign type 4 atoms to elements
        xfc.atom_node()  # assign type 2 atoms to nodes
        xfc.displacement_interpolation()  # calculate strain on type 4 atoms
        # create IMD input file with updates positions of type 4 atoms
        xfc.create_atom_dis()  # create dislocation in atomic core
        self.relax_atoms(name='init_dislocation')
        print('\n Created screw dislocation in model.\n')
        
    def atom_bc(self):
        '''
        Get displacement of relaxed atoms and apply as BC to XFEM model
        Reads atomistic configuration from file: 
            relaxed_atomistic_dislocation_structure.00000.ss
        Atomic displacements are calculated relative to:
            relaxed_perfect_crystal_with_atom_type.imd
            
        Returns
        -------
        None.

        '''
        xfc.nodal_displacement()  # calculate BC from relaxed atomic positions
        hh =  np.array(xfc.at_disp, dtype=np.double)
        self.at_disp =hh[1:self.natom+1]
        
    def shift_atoms(self):
        '''
        Apply nodal positons to type 4 atoms and write updated atomic
        configuration
        
        Reads atomistic configuration from file: 
            relaxed_atomistic_dislocation_structure.00000.ss
        Produces file:
            atomistic_dislocation_with_fem_solution.imd

        Returns
        -------
        None.

        '''
        # calculate shift of type 4 atoms due to strain in elements
        xfc.displacement_interpolation()
        xfc.atom_configuration()  # produce updated IMD input file
        
    def relax_atoms(self, i=0, name='relaxed_atoms'):
        '''
        Run IMD-FIRE to relax atomic structure with fixed boundary atoms
        Reads atomistic configuration from file: 
            atomistic_dislocation_with_fem_solution.imd
        Produces file: relaxed_atomistic_dislocation_structure.00000.ss
        Parameters
        ----------
        i : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        t0 = time.time()
        os.system('cd {0}; ./imd_eam_fire_homdef_stress_nbl -p Fe101-fire-disloc-sample.param'
                  .format(self.md_dir))
        t1 = time.time()
        if os.path.isfile('{}/relaxed_atomistic_dislocation_structure.00000.ssitr'
                     .format(self.temp)):
            print('Relaxation of atoms successful in {:6.5}s.'.format(t1-t0))
            i1 = int(i/10)
            i2 = i % 10
            # store relaxed config
            os.system('cp {0}/relaxed_atomistic_dislocation_structure.00000.ss {0}/{1}.{2}{3}.ss'\
                      .format(self.temp, name, i1, i2))
        else:
            raise RuntimeError('Relaxation of atoms not successful in iteration {}.'
                               .format(i))
        os.system('rm {0}/*eng {0}/*itr {0}/*ssdef {0}/*chkpt; cd ..'\
                  .format(self.temp))
            
    def apply_bc(self, val, comp='yz', bc_type='stress'):
        '''
        

        Parameters
        ----------
        val : TYPE
            DESCRIPTION.
        comp : TYPE, optional
            DESCRIPTION. The default is 'yz'.
        bc_type : TYPE, optional
            DESCRIPTION. The default is 'stress'.

        Returns
        -------
        None.

        '''
        cdef double e23
        
        if bc_type == 'stress':
            e23 = val/self.C44
        elif bc_type == 'strain':
            e23 = val
        else:
            raise ValueError('Unsupported parameter bc_type: {}'
                             .format(bc_type))
        xfc.apply_e23_outer(e23)
        self.solve()
        self.shift_atoms()
        self.ubc = np.array(xfc.ENFRDISPglob[3:], dtype=np.double)
        
    
    def calc_stress(self):
        '''
        Evaluate element stresses.

        Returns
        -------
        s_field : (NEL, 3) array
            Components of stress field in each element

        '''
        cdef int IEL, j
        cdef double XLOC[9]
        cdef double YLOC[9]
        cdef double ZLOC[9]
        XLOC[0] = 0.0
        YLOC[0] = 0.0
        ZLOC[0] = 0.0
        s_field = np.zeros((self.NEL, 3), dtype=np.double)
        for IEL in range(self.NEL):
            for j in range(8):
                XLOC[j+1] = self.nodes[self.elmts[IEL, j], 0]
                YLOC[j+1] = self.nodes[self.elmts[IEL, j], 1]
                ZLOC[j+1] = self.nodes[self.elmts[IEL, j], 2]
            xfc.STRESS(XLOC, YLOC, ZLOC, IEL)
            s_field[IEL,:] = xfc.sigma
        return s_field
        
            
    def plot_nodal(self, tag, comp='z', layer=None, atoms=True):
        '''
        Plot given component of nodal siolution.

        Parameters
        ----------
        tag : TYPE
            DESCRIPTION.
        comp : TYPE, optional
            DESCRIPTION. The default is 'z'.

        Returns
        -------
        None.

        '''
        # select component to be plotted
        if type(comp)==int:
            if comp<0 or comp>2:
                raise ValueError('Parameter comp mus be in range 0...2, not {}.'
                                 .format(comp))
            istart = comp
        elif comp.lower() == 'z':
            istart = 2
        elif comp.lower() == 'y':
            istart = 1
        elif comp.lower() == 'x':
            istart = 0
        else:
            raise ValueError('Unknown value for parameter comp: {}'
                             .format(comp))
        # select layer of nodes to be plotted
        zc = np.unique(self.nodes[:, 2])
        if layer is None:
            layer = int(len(zc)/2)
        ind = np.nonzero(np.abs(self.nodes[:,2] - zc[layer]) < 1.e-4)[0]
            
        # select quantity to be plotted
        if tag == 'ubc':
            self.ubc = np.array(xfc.ENFRDISPglob[3:], dtype=np.double)
            quant = self.ubc[istart + ind*3]
            title = 'Boundary conditions for u_{}'.format(comp)
        elif tag == 'u':
            quant = self.u[istart + ind*3]
            title = 'Nodal solution for displacement u_{}'.format(comp)
        elif tag == 'f':
            quant = self.f[istart + ind*3]
            title = 'Nodal solution for force f_{}'.format(comp)
        else:
            raise ValueError('Unknown value for parameter tag: {}')

        plt.scatter(self.nodes[ind, 0], self.nodes[ind, 1],
                    c=quant, marker=',')
        plt.colorbar()
        plt.title(title)
        
        if atoms and tag == 'u':
            ind = np.nonzero(np.abs(self.apos[:,2] - zc[layer]) < 0.3)[0]
            if len(ind) == 0:
                raise ValueError('XFEM mesh and nodal positions do not conform.')
            plt.scatter(self.apos[ind, 0], self.apos[ind, 1],
                        c=self.at_disp[ind, istart], marker=',')
            
        plt.show()
        
    def plot_el(self, tag, comp='yz', sig=None):
        if comp == 'yz':
            sc = 2
        elif comp == 'xx':
            sc = 0
        elif comp == 'yy':
            sc = 1
        else:
            raise ValueError('Value for parameter comp not valid: {}'
                             .format(comp))
        if sig is None:
            sig = self.calc_stress()  # evaluate element stresses
        iz = 1

        plt.figure()
        plt.scatter(self.el_ctr[iz::3, 0],  self.el_ctr[iz::3, 1],
                    c=sig[iz::3, sc], marker=',')
        plt.colorbar()
        plt.title('Stress component {} at z={:5.4}'
                  .format(sc, self.el_ctr[iz, 2]))
        plt.show()
