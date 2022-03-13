# distutils: language = c++
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 17:46:17 2022

@author: alexander
"""

import os, warnings
import time
import numpy as np
import pyvista as pv
from vtk import VTK_HEXAHEDRON
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from xfeat cimport cpp_wrapper as xfc
from xfeat.basic import rot_elast_tens

class Model(object):
    def __init__(self, mat, size=500):
        self.fem_size = size  # size of FEM part
        self.bv = None  # Burgers vector
        self.dist = None  # reference distances of atoms along box edges
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
        self.update_stress = True
        
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
        self.nu = self.C12 / (self.C11 + self.C12)
        xfc.nu = self.nu
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
        
        # set global constants
        xfc.fem_size = self.fem_size
        xfc.PI = np.pi
        
        DIR1 = np.zeros(9, dtype=np.double)
        DIR2 = np.zeros(9, dtype=np.double)
        DIR3 = np.zeros(9, dtype=np.double)
        DIR1[1] = DIR2[4] = DIR3[3] = -1.0
        DIR1[2] = DIR2[8] = DIR3[4] = -1.0
        DIR1[3] = DIR2[3] = DIR3[7] = 1.0
        DIR1[4] = DIR2[7] = DIR3[8] = 1.0
        DIR1[5] = DIR2[1] = DIR3[2] = -1.0
        DIR1[6] = DIR2[5] = DIR3[1] = -1.0
        DIR1[7] = DIR2[2] = DIR3[5] = 1.0
        DIR1[8] = DIR2[6] = DIR3[6] = 1.0
        xfc.DIR1 = DIR1
        xfc.DIR2 = DIR2
        xfc.DIR3 = DIR3
	    	    
        A = np.zeros((7, 10), dtype=np.double)
        A[1][1] = A[2][5] = A[3][9] = 1.0
        A[4][2] = A[4][4] = 1.0
        A[5][6] = A[5][8] = 1.0
        A[6][3] = A[6][7] = 1.0
        xfc.A = A
        
        # set global parameters for pyVista
        pv.global_theme.font.title_size = 24
        pv.global_theme.font.label_size = 20
        self.e_names = ['eps_xx', 'eps_yy', 'eps_zz', 'eps_xy', 'eps_yz', 'eps_xz']
        self.s_names = ['sig_xx', 'sig_yy', 'sig_zz', 'sig_xy', 'sig_yz', 'sig_xz']
        self.comps = ['xx', 'yy', 'zz', 'xy', 'yz', 'xz']
        
    def atoms(self, size=None):
        '''
        Set up atomic core. 
        Produces IMD input file with atom types for boundary conditions.

        Returns
        -------
        res : TYPE
            DESCRIPTION.
            
        Attributes
        ----------
        lp : float
            Lattice parameter after relaxtion
        shift : (2,) array
            Shift between FEM and atomistic origins

        '''
        # create relaxed perfect crystal as IMD input (Steps 2 and 3)
        if size is None:
            size = [10, 17, 3]
        name_urel = self.temp + '/unrelaxed_perfect_crystal.imd'
        sx = size[0] * self.lp * np.linalg.norm(self.ori[0,:])  # 10 
        sy = size[1] * self.lp * np.linalg.norm(self.ori[1,:])  # 17
        sz = size[2] * self.lp * np.linalg.norm(self.ori[2,:])
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
        self.dist = xfc.dist  # atomic distances along Cartesian coordinates
        """self.bv = self.dist[2]  # Burgers vector

        if np.abs(2.*self.bv/np.sqrt(3.) - self.lp) > 1.e-9:
            print('Correcting lattice parameter in material definition after atomic relaxation.')
            print('Given value: {}'.format(self.lp))
            self.lp = 2.*self.bv/np.sqrt(3.)
            print('New value: {}'.format(self.lp))"""
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
        
        # create pyVista atomistic grid
        self.atom_grid = pv.PolyData(self.apos)
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
        
        if self.dist is None:
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
        
        # create pyVista XFEM grid
        cells = np.ones((self.NEL, 9), dtype=int)*8
        cells[:,1:9] = self.elmts
        celltypes = np.empty(self.NEL, dtype=np.uint8)
        celltypes[:] = VTK_HEXAHEDRON
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
        self.grid.point_data['u_x'] = self.u[0::3]
        self.grid.point_data['u_y'] = self.u[1::3]
        self.grid.point_data['u_z'] = self.u[2::3]
        self.grid.point_data['f_x'] = self.f[0::3]
        self.grid.point_data['f_y'] = self.f[1::3]
        self.grid.point_data['f_z'] = self.f[2::3]
        self.grid.point_data['ubc_x'] = self.ubc[0::3]
        self.grid.point_data['ubc_y'] = self.ubc[1::3]
        self.grid.point_data['ubc_z'] = self.ubc[2::3]
        self.update_stress = True
        return
    
    def init_dislo(self, bdir=None):
        '''
        Introduce screw dislocation and write updated atomic configuration

        Returns
        -------
        None.

        '''
        # create Volterra screw dislocation by introducing a shift in
        # XFEM elements and atomic core
        if bdir is None:
            print('No Burgers vector given, creating screw dislocation.')
            self.b_vec = np.zeros(3, dtype=np.double)
            self.b_vec[2] = 1.
            self.bv = self.dist[2]
        else:
            if len(bdir) != 3:
                raise ValueError('3 Components must be given for Burgers vector direction.')
            self.b_vec = np.array(bdir, dtype=np.double)
            ind = np.nonzero(self.b_vec)[0]
            if len(ind) != 1:
                raise ValueError('Only one component for Burgers vector supported.')
            if not np.isclose(np.linalg.norm(self.b_vec), 1.):
                self.b_vec /= np.linalg.norm(self.b_vec)
                warnings.warn('Burgers vector director (bdir) should be given as unit vector.')
            self.bv = self.dist[ind[0]] 
        xfc.et = self.b_vec
        xfc.bv = self.bv
        xfc.create_xfem_dis()  # create displacements in XFEM region
        self.ubc = np.array(xfc.ENFRDISPglob[3:], dtype=np.double)
        self.grid.point_data['ubc_x'] = self.ubc[0::3]
        self.grid.point_data['ubc_y'] = self.ubc[1::3]
        self.grid.point_data['ubc_z'] = self.ubc[2::3]
        self.solve()  # calculate nodal displacements under BC of shift
        xfc.atom_element()  # assign type 4 atoms to elements
        xfc.atom_node()  # assign type 2 atoms to nodes
        xfc.displacement_interpolation()  # calculate strain on type 4 atoms
        # create IMD input file with updates positions of type 4 atoms
        xfc.create_atom_dis()  # create dislocation in atomic core
        self.relax_atoms(name='init_dislocation')
        print(f'\n Created dislocation with Burgers vector ({self.b_vec}) in model.')
        print(f'Norm of Burgers vector is {self.bv} A')
        
    def atom_bc(self):
        '''
        Get displacement of relaxed atoms and apply as BC to XFEM model
        Get potential energies of atoms
        Reads atomistic configuration from file: 
            relaxed_atomistic_dislocation_structure.00000.ss
        Atomic displacements are calculated relative to:
            relaxed_perfect_crystal_with_atom_type.imd
            
        Returns
        -------
        None.

        '''
        xfc.nodal_displacement()  # calculate BC from relaxed atomic positions
        self.ubc = np.array(xfc.ENFRDISPglob[3:], dtype=np.double)
        self.grid.point_data['ubc_x'] = self.ubc[0::3]
        self.grid.point_data['ubc_y'] = self.ubc[1::3]
        self.grid.point_data['ubc_z'] = self.ubc[2::3]
        hh =  np.array(xfc.at_disp, dtype=np.double)
        self.atom_grid.point_data['u_x'] = hh[1:self.natom+1, 0]
        self.atom_grid.point_data['u_y'] = hh[1:self.natom+1, 1]
        self.atom_grid.point_data['u_z'] = hh[1:self.natom+1, 2]
        hh = np.array(xfc.at_energy, dtype=np.double)
        self.atom_grid.point_data['epot'] = hh[1:self.natom+1]
        self.int_at_node = xfc.interaction_atom_node
        self.Ntype2 = xfc.NCONSNODE
        
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
        cdef double eps
        
        if bc_type == 'stress':
            eps = val/self.C44
        elif bc_type == 'strain':
            eps = val
        else:
            raise ValueError('Unsupported parameter bc_type: {}'
                             .format(bc_type))
        if comp == 'yz':
            xfc.apply_e23_outer(eps)
        elif comp == 'xy':
            xfc.apply_e12_outer(eps)
        self.ubc = np.array(xfc.ENFRDISPglob[3:], dtype=np.double)
        self.grid.point_data['ubc_x'] = self.ubc[0::3]
        self.grid.point_data['ubc_y'] = self.ubc[1::3]
        self.grid.point_data['ubc_z'] = self.ubc[2::3]
        #self.solve()
        #self.shift_atoms()

    def calc_stress(self):
        '''
        Evaluate element stresses.

        Returns
        -------
        s_field : (NEL, 3) array
            Components of stress field in each element

        '''
        cdef double xloc[9], yloc[9], zloc[9]
        cdef int iel, j
        xloc = np.zeros(9, dtype=np.double)
        yloc = np.zeros(9, dtype=np.double)
        zloc = np.zeros(9, dtype=np.double)
        s_field = np.zeros((self.NEL, 6), dtype=np.double)
        e_field = np.zeros((self.NEL, 6), dtype=np.double)
        for iel in range(self.NEL):
            for j in range(8):
                xloc[j+1] = self.nodes[self.elmts[iel, j], 0]
                yloc[j+1] = self.nodes[self.elmts[iel, j], 1]
                zloc[j+1] = self.nodes[self.elmts[iel, j], 2]
            xfc.calc_stress(xloc, yloc, zloc, iel)
            s_field[iel,:] = xfc.stress
            e_field[iel,:] = xfc.strain
        s_field /= 0.0062417  # convert from eV/A^3 into GPa
        for i in range(6):
            self.grid.cell_data[self.s_names[i]] = s_field[:, i]
            self.grid.cell_data[self.e_names[i]] = e_field[:, i]
        self.update_stress = False
        return s_field
        
            
    def plot_nodal(self, tag, comp='z', atoms=True, deformed=None):
        '''
        Plot given component of nodal solution.

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
        comp = comp.lower()
        tag = tag.lower()
        if tag[-1] == 'z' or tag[-1] == 'y' or tag[-1] == 'x':
            comp = tag[-1]
            tag = tag[0:-1]
        elif comp != 'z' and comp != 'y' and comp != 'x':
            raise ValueError('Unknown value for parameter comp: {}'
                             .format(comp))
            
        # select quantity to be plotted
        if tag == 'ubc':
            title = r'BC u_{} (A)'.format(comp)
            sc = 'ubc_{}'.format(comp)
        elif tag == 'u':
            title = r'u_{} (A)'.format(comp)
            sc = 'u_{}'.format(comp)
        elif tag == 'f':
            title = r'f_{} (1e-11 N)'.format(comp)
            sc = 'f_{}'.format(comp)
        else:
            raise ValueError('Unknown value for parameter tag: {}')
        
        pl = pv.Plotter()
        pl.camera_position = (1.0, 0.0, 2.1*self.fem_size)
        pl.camera.azimuth = 270
        self.grid.set_active_scalars(sc)
        umin = np.amin(self.grid.active_scalars)
        umax = np.amax(self.grid.active_scalars)
        if atoms and tag == 'u':
            self.atom_grid.set_active_scalars(sc)
            umin = np.minimum(umin, np.amin(self.atom_grid.active_scalars))
            umax = np.maximum(umax, np.amax(self.atom_grid.active_scalars))
            pl.add_mesh(self.atom_grid, style='points', point_size=10,
                        render_points_as_spheres=True, show_scalar_bar=False,
                        clim=(umin, umax))
        pl.add_mesh(self.grid, show_edges=True, clim=(umin, umax),
                    scalar_bar_args=dict(vertical=True, position_y=0.25,
                                         title=title))
        pl.show()
        return
        
    def plot_el(self, tag, comp='yz'):
        tag = tag.lower()
        comp = comp.lower()
        if tag[-2:] in self.comps:
            comp = tag[-2:]
            tag = tag[0:-2]
        if tag != 'sig' and tag != 'eps':
            raise ValueError('Value for parameter "tag" not valid: {}. Must be eiter "sig" or "eps".'
                             .format(tag))
            
        if not comp in self.comps:
            raise ValueError('Value for parameter "comp" not valid: {}'
                             .format(comp))
        sc = '{}_{}'.format(tag, comp)
        title = '{}_{}'.format(tag, comp)
        if tag == 'sig':
            title += ' (GPa)'
        if self.update_stress:
            sig = self.calc_stress()  # evaluate element stresses

        self.grid.set_active_scalars(sc)
        self.grid.plot(cpos='xy', show_edges=True,
                    scalar_bar_args=dict(vertical=True, position_y=0.25,
                                         title=title))
        return
    
    def plot_at(self, tag, deformed=None):
        tag = tag.lower()
        if tag == 'epot':
            sc = 'epot'
            title = r'Epot (eV)'
        elif tag == 'uz' or tag == 'dispz':
            sc = 'u_z'
            title = r'u_z (A)'
        elif tag == 'uy' or tag == 'dispy':
            sc = 'u_y'
            title = r'u_y (A)'
        elif tag == 'ux' or tag == 'dispx':
            sc = 'u_x'
            title = r'u_x (A)'
        else:
            raise ValueError('Unknown value in parameter tag: {}'.format(tag))
            
        self.atom_grid.set_active_scalars(sc)
        vmin = np.amin(self.atom_grid.active_scalars)
        vmax = np.amax(self.atom_grid.active_scalars)
        if sc == 'epot':
            vmax = vmin*0.95
        # pl = pv.Plotter()
        # pl.camera_position = (1.0, 0.0, 2*self.fem_size)
        # pl.camera.azimuth = 270
        # pl.add_mesh(self.atom_grid, style='points', point_size=20,
        #             render_points_as_spheres=True)
        # pl.show()
        self.atom_grid.plot(cpos='xy', style='points', point_size=20,
                    render_points_as_spheres=True, clim=(vmin, vmax),
                    scalar_bar_args=dict(vertical=True, position_y=0.25,
                                 title=title))

    def plot(self, tag, atoms=None, deformed=None):
        tag = tag.lower()
        if atoms is None:
            atoms = True
        if tag == 'epot':
            self.plot_at(tag, deformed=deformed)
        elif tag == 'uz' or tag == 'uy' or tag == 'ux':
            self.plot_nodal(tag, atoms=atoms, deformed=deformed)
        elif tag == 'fz' or tag == 'fy' or tag == 'fx':
            self.plot_nodal(tag, atoms=atoms, deformed=deformed)
        elif tag == 'ubcz' or tag == 'ubcy' or tag == 'ubcx':
            self.plot_nodal(tag, atoms=atoms, deformed=deformed)
        elif tag == 'dispx' or tag == 'dispy' or tag == 'dispz':
            self.plot_at(tag, deformed=deformed)
        elif tag == 'sigxx' or tag == 'sigyy' or tag == 'sigzz' or\
             tag == 'sigxy' or tag == 'sigyz' or tag == 'sigxz':
            self.plot_el(tag='sig', comp=tag[-2:])
        elif tag == 'epsxx' or tag == 'epsyy' or tag == 'epszz' or\
             tag == 'epsxy' or tag == 'epsyz' or tag == 'epsxz':
            self.plot_el(tag='eps', comp=tag[-2:])
        else:
            raise ValueError('Unknown value in parameter tag.')
