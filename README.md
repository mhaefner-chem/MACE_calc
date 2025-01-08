# MACE_calc
![GitHub Release](https://img.shields.io/github/v/release/mhaefner-chem/MACE_calc?include_prereleases) ![GitHub License](https://img.shields.io/github/license/mhaefner-chem/MACE_calc)

MACE_calc is a tool based on ASE that uses MACE models as calculators to carry out various basic kinds of atomistic calculations.

## Table of contents
- [How to use?](#how-to-use)
- [Requirements & Installation](#requirements-and-installation)

## How to use?
MACE_calc can be invoked with 
```console
MACE_calc.py -c <COMPOUND> -s <SETTINGS_FILE> <OPTIONS>
```
or
```console
MACE_calc.py -l <LIST_OF_COMPOUNDS> -s <SETTINGS_FILE> <OPTIONS>
```
The most basic options are:
-  -h [ --help ]                Prints a help message.
-  -c [ --compound ] <ARG>      Input structure <ARG> formatted as CIF, VASP, or (EXT)XYZ.
-  -l [ --list ] <ARG>          <ARG> is a list of valid input structures.
-  -m [ --model ] <ARG>         <ARG> specifies the MACE model used for the calculation. Default = macemp0.
-  -p [ --proc ] <ARG>          <ARG> specifies the mode that will be used.  
  Modes:
    - sp: single point  
    - opt: optimization  
    - phon: phonon calculation  
    - md: MD simulation  
    - bulk: calculation of bulk modulus  
    - neb: NEB calculation  

-  -s [ --settings ] <ARG>      <ARG> contains more detailed settings for the calculation.
-  -d [ --d3bj ]                Activates Grimme D3(BJ) dispersion correction.
-  -v [ --verbosity ] <ARG>     Changes the verbosity of the output. Default = 1.
-  -g [ --gpu ] <ARG>           Use a GPU for calculations instead of CPU. Default = False.
-  --nosym                      Symmetry is not enforced in calculations.
-  --mol                        Deactivates periodic boundary conditions for calculation of molecules.

E.g.,
```console
MACE_calc.py -c NaCl.cif -p opt -d
```
carries out a structure optimization of the file `NaCl.cif` with  DFTD3 dispersion correction.

Running the program without any specifications, i.e.,
```console
MACE_calc.py
```
yields a full list of program settings that can be specified in a settings file, subdivided by general settings and settings for the various procedures.

## Requirements and installation

So far, the program has been successfully tested with `python 3.11`, `ase 3.23.0`, and `mace-torch 0.3.4`. 
The source code for the program is obtained with the command

```console
$ git clone https://github.com/mhaefner-chem/MACE_calc
```

Running the program with python requires at least the python packages [`ase`](https://wiki.fysik.dtu.dk/ase/install.html), [`mace-torch`](https://github.com/ACEsuit/mace), and all their dependencies.

### Additional Notes

The python package `spglib` is required to take into account the symmetry of the investigated structures. Instructions for installing spglib can be found at its [github page](https://github.com/spglib/spglib) or the [documentation](https://spglib.readthedocs.io) of spglib.

Except for the MACE-MP foundational model, dispersion corrections are not included in the base installation. However, it supports two kinds of DFTD3 dispersion correction:
(`torch-dftd`)[https://github.com/pfnet-research/torch-dftd] (Caveat: no longer actively maintained)
(`dftd3-python`)[https://dftd3.readthedocs.io/en/latest/installation.html]

MACE_calc also allows for computations to be run on GPUs, but in order to do this, the CUDA version of pytorch has to be installed during the installation of `mace-torch`.


