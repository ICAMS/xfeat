# distutils: language = c++
# distutils: libraries = c++0x
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 17:46:17 2022

@author: alexander

"""
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "cpp_wrapper.cpp":
    pass

cdef extern from "cpp_wrapper.h":
    # Declaration of all variables and functions used jointly in CPP and Cython parts.
    
    # constants and model parameters
    double PI
    string temp_dir  # absolute path to temp directory, defined in Model.__init__
    double fem_size  # size of FEM geometry, defined in Model.__init__
    double EE3D[7][7]  # rotatated material stiffness tensor, defined in Model.mesh
    
    # atomistic quantities
    int natom  # number of atoms, defined in set_up_atoms
    double bv  # magnitude of Burgers vector, defined in atom_set_up
    double shift[2]  # shift between atmomistic coords and nodal positions
    double coords[20000][3]  # atomic positions in reference configuration, defined in set_up_atoms
    double at_disp[20000][3]  # atomic displacements, defined in nodal_displacement
    double at_energy[20000]  # potential energy of atoms
    
    # XFEM quantities
    int NNODE  # number of nodes, defined in create_mesh
    int NLAG  # number of augmented elmts (Lagrangians)
    int NEL  # number of elements
    int NDF  # no. of degrees of freedom
    int xdof  # no. of extended DOF
    int nelem_full_big_box_x  # number of elements along x-axis
    int nelem_full_big_box_y
    int nelem_full_big_box_z
    double DIR1[9], DIR2[9], DIR3[9]
    double A[7][10]
    double et[4]
    double Lx, Ly, Lz  # lengths of axes of atomic box
    double stress[6]  # Voigt stress tensor in element, calulated in calc_stress
    double strain[6]  # Voigt strain tensor in element, calulated in calc_stress
    vector[int] Iglob  # indices of compressed stiffness matrix
    vector[int] Jglob  # columnes of compressed stiffness matrix
    vector[int] JJglob  # number of equation in global system per node (negative of no DOF for node)
    vector[int] LOTOGO  # nodes associated to elements (8 per element)
    vector[double] NEWAglob  # system stiffness matrix
    vector[double] Dglob  # global displacement vector
    vector[double] Fglob  # global force vector
    vector[double] ENFRDISPglob  # field of enforced displacements (boundary conditions)
    vector[double] Xglob, Yglob, Zglob  # nodal positions
    vector[double] node_dis  # nodal displacements, defined in Model.solve

    # functions
    void atom_set_up()  # read initial atomic structure
    void create_mesh()  # define FE mesh
    void init_matrix()  # initialize system stiffness matrix in compressed format
    void update_fglob()  # update global force vector
    void atom_node()  # assign type 2 atoms to nodes on inner boundary
    void atom_element()  # assign type 4 atoms to elements in overlap region
    void displacement_interpolation()  # interpolate displacements in elements in overlap region
    void atom_configuration()  # transfer strains from XFEM solution to type 4 atoms
    void create_atom_dis()  # create a screw dislocation in atomic core
    void nodal_displacement()  # calculate displacements on inner boundary nodes from type 2 atoms
    void create_xfem_dis()  # create displacements for Volterra dislocation in XFEM part
    void apply_e23_outer(double e23)  # apply shear strain on outer XFEM boundary
    void calc_stress(double XLOC[9], double YLOC[9], double ZLOC[9], int IEL)  # evaluate element stress
