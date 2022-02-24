# distutils: language = c++
# distutils: libraries = c++0x
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 17:46:17 2022

@author: alexander
"""
cdef extern from "<vector>" namespace "std":
    cdef cppclass vector[T]:
        cppclass iterator:
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        vector()
        void push_back(T&)
        T& operator[](int)
        T& at(int)
        iterator begin()
        iterator end()

cdef extern from "cpp_wrapper.cpp":
    pass

cdef extern from "cpp_wrapper.h":
    double PI
    # Atom_set_up definitios
    double shift[2]
    double dist[3]
    double Lx
    double Ly
    double Lz
    double bv
    int natom
    int *atom_id
    double *coords
    double *masses
    double minwidth_x
    double minwidth_y
    double fem_size
    double eps
    #int nodecount
    int nelem_full_big_box_x
    int nelem_full_big_box_y
    int nelem_full_big_box_z
    void get_n_atoms()
    void atom_set_up()

    # FEM_Mesh definitions
    void create_mesh()
    void init_matrix()
    void iter_matrix()
    void apply_e23_outer(double e23)
    void STRESS(double XLOC[9], double YLOC[9], double ZLOC[9], int IEL)
    int xdof
    #int size_Aglob
    int NLAG
    int NEL
    int NDF  # no. of degrees of freedom
    long int NBIG  # no. of co-effients in the skyline matrix
    int NNODE, NCF, NENFD, NCONT
    
    double EE3D[7][7]
    double sigma[3]
    vector[int] JJglob
    vector[double] NEWAglob
    vector[double] Dglob
    vector[int] Iglob
    vector[int] Jglob
    vector[int] LOTOGO
    vector[double] Fglob
    vector[double] ENFRDISPglob
    vector[double] Xglob, Yglob, Zglob
    vector[double] node_dis
    vector[int] ContraintNodes

    # Init_Dislocation definitions
    void atom_element()
    void displacement_interpolation()
    void atom_configuration()
    void init_atom_config()
    void nodal_displacement()
    void create_volterra_dis()
    void atom_node()
    
    # Iter_Step definitions
    void displacement_interpolation()
    void atom_configuration()

