# XFEAt

### eXtended Finite Element Analysis with Atomistics

  - Authors: Karthikeyan Chockalingam, Alexander Hartmaier
  - Organization: ICAMS, Ruhr University Bochum, Germany
  - Contact: <alexander.hartmaier@rub.de>


XFEAt is a tool to combine atomistic simulation with boundary conditions from eXteneded Finite Element Method (XFEM). Useful for studying properties of dislocations under well-defined mechanical boundary conditions.

## Installation
The XFEAt package requires an [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) environment with a recent Python version. It can be installed from its GitHub repository via the following steps:

```
$ git clone https://github.com/ICAMS/xfeat ./XFEAt
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

might work.

To test if the installation has been successful, run

```
$ pytest tests
```

## Examples
An example that generates a screw dislocation in a bcc iron crystal and exposes it to a shear stress can be found under "examples/FE\_screw\_dislocation". 

## License

The XFEAt package comes with ABSOLUTELY NO WARRANTY. This is free
software, and you are welcome to redistribute it under the conditions of
the GNU General Public License
([GPLv3](http://www.fsf.org/licensing/licenses/gpl.html))

The contents of the examples and notebooks are published under the 
Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
([CC BY-NC-SA 4.0](http://creativecommons.org/licenses/by-nc-sa/4.0/))
