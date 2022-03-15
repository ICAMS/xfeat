# XFEAt

### eXtended Finite Element Analysis with Atomistics

  - Authors: Karthikeyan Chockalingam, Alexander Hartmaier
  - Organization: ICAMS, Ruhr University Bochum, Germany
  - Contact: <alexander.hartmaier@rub.de>


XFEAt is a tool to combine atomistic simulation with boundary conditions from eXteneded Finite Element Method (XFEM). Useful for studying properties of dislocations under well-defined mechanical boundary conditions.

## Installation
XFEAt is written in [Cython](https://cython.org) and C++ and requires a C++ compiler. Cython will be installed together with all other dependencies. The Python interface requires an [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) environment with a recent Python version. If a Conda environment and C++ compiler are available, the XFEAt package can be installed from its GitHub repository via the following steps:

```
$ git clone https://github.com/ICAMS/xfeat ./XFEAt
$ cd XFEAt
$ conda env create -f environment.yml
$ conda activate xfeat 
$ make install
```

Furthermore, an executable of the [ITAP Molecular Dynamics Program](http://imd.itap.physik.uni-stuttgart.de) (IMD) is required. This software is licensed under the GNU General Public License GPLv3 and a copy of IMD is included in the XFEAt distribution in "libs/imd". See the [IMD user guide] (http://imd.itap.physik.uni-stuttgart.de/userguide/compiling.html) for ways to build the executable "imd\_eam\_fire\_homdef\_stress\_nbl", which should be placed in the directory "XFEAt/Fe\_MD". On many systems something like

```
$ cd libs/imd/src
$ make IMDSYS=P4-gcc3 BIN_DIR="../../../Fe_MD" imd_eam_fire_homdef_stress_nbl
$ cd ../../..
```

might work, although it produces a rather slow executable. For larger systems the optimization of compiler flags is strongly recommended.

To test if the installation has been successful, run

```
$ pytest tests
```

## Examples
Examples that generate a screw or an edge dislocation in a bcc iron crystal and exposes them to a shear stress can be found under "examples/FE\_dislocations". 

## Dependencies
XFEAt requires the following packages as imports:

 - [Cython](https://cython.org) as programming language
 - [NumPy](http://numpy.scipy.org) for array handling
 - [Scipy](https://www.scipy.org/) for numerical solutions, in particular for sparse matrix solving
 - [pyVista](https://docs.pyvista.org) for visualization of numerical reuslults on XFEM grid and atomistic lattice

## License
The XFEAt package comes with ABSOLUTELY NO WARRANTY. This is free
software, and you are welcome to redistribute it under the conditions of
the GNU General Public License
([GPLv3](http://www.fsf.org/licensing/licenses/gpl.html)).

The contents of the examples and notebooks are published under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
([CC BY-NC-SA 4.0](http://creativecommons.org/licenses/by-nc-sa/4.0/)).
