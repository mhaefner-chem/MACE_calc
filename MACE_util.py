#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:41:06 2024

@author: bt308570

utilities for MACE_calc
"""

def print_separator(symbol,n=32,skip=True):
    if skip == True:
        print("")
    for i in range(n):
        print(symbol,end="")
    print("")
    if skip == True:
        print("")
        
        
def write_results(file,parameters,mode):
    with open(file, mode=mode) as f:
        f.write(parameters)
        
# gets the length of a vector
import numpy as np
def get_len_vec(vec):
    vec_len = 0
    for i in range(len(vec)):
        vec_len += (vec[i])**2
    vec_len = np.sqrt(vec_len)
    return vec_len
        
# analyzes symmetry
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry
def get_symmetry(structure,verbosity = 0, prec=0.01):
    if verbosity > 0:
        print("Given symmetry at precision {}:".format(prec))
        sym = check_symmetry(structure, prec, verbose=False)
        print("SG {} Hermann-Maughin {}, Hall {}".format(sym["number"],sym["international"],sym["hall"]))
    structure_sym = structure.copy()
    structure_sym.set_constraint(FixSymmetry(structure))
    return structure_sym

# prints results from calculation

def print_results(settings,compound,structure,verbosity=1):
    
    results = structure.calc.results
    E_0 = structure.calc.get_potential_energy() #results["energy"]
    n_atoms = len(structure.symbols)
    positions = structure.get_positions()
    frac_pos = structure.get_scaled_positions()
    elements = structure.get_chemical_symbols()
    
    print("Energy: {:9.3f} eV, {:6.3f} eV/atom".format(E_0,E_0/n_atoms))
    
    cart = ["x","y","z"]
    
    if settings.print_geometries == True and verbosity > -1 or verbosity == 2:
        print_separator("-")
        format = "{:>8} {:>8} {:>8} {:>8} {:>8} {:>8}"
        
        print(format.format("a/Å","b/Å","c/Å","α/°","β/°","γ/°"))
        for i in range(6):
            print("{:8.4f}".format(structure.cell.cellpar()[i]),end=" ")
        print()    
        print()
        
        cubic_div = get_cubic_divergence(structure.cell.volume,structure.cell.cellpar()[0:2])
        
        print("V/Å³: "+"{:10.3f}".format(structure.cell.volume))
        print("Δ_cube:  ","{:6.3f}".format(cubic_div))
        print()
    
        print("Positions:")
        print("")
        print("{:8} {:^32} {:^32}".format("","Cartesian Coordinates","Fractional Coordinates"))
        print("{:5} {:2} {:^10} {:^10} {:^10} {:^10} {:^10} {:^10}".format("","","x","y","z","x","y","z"))
        print("{:>5} {:2} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}".format("ID","El","Å","Å","Å","Å","Å","Å"))
        format = "{:5} {:2} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f}"
        for i in range(len(positions)):
            p = positions[i]
            f = frac_pos[i]
            e = elements[i]
            
            print(format.format(i,e,p[0],p[1],p[2],f[0],f[1],f[2]))
       
    
    if settings.print_stresses == True and verbosity > -1 or verbosity == 2:
        stress = results["stress"]
        
        
        print_separator("-")
        print("Stress Tensor [eV/Å²]:")
        print("")
        print("{:>8} {:>8} {:>8} {:>8}".format("",cart[0],cart[1],cart[2]))
        print("{:>8} {:8.5f} {:8.5f} {:8.5f}".format(cart[0],stress[0],stress[5],stress[4]))
        print("{:>8} {:8.5f} {:8.5f} {:8.5f}".format(cart[1],stress[5],stress[1],stress[3]))
        print("{:>8} {:8.5f} {:8.5f} {:8.5f}".format(cart[2],stress[4],stress[3],stress[2]))
        
    if settings.print_forces == True and verbosity > -1 or verbosity == 2:
        forces = results["forces"]
        
        print_separator("-")
        print("Forces:")
        print("")
        print("{:5} {:2} {:^10} {:^10} {:^10} {:^10} {:^10} {:^10} {:^10}".format("","","x","y","z","F(x)","F(y)","F(z)","F_tot"))
        print("{:>5} {:2} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}".format("ID","El","Å","Å","Å","eV/Å","eV/Å","eV/Å","eV/Å"))
        format = "{:5} {:2} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f} {:10.4f}"
        for i in range(len(forces)):
            p = positions[i]
            f = forces[i]
            e = elements[i]
            
            
            print(format.format(i,e,p[0],p[1],p[2],f[0],f[1],f[2],get_len_vec(f)))
            

    print_separator("-")
    
    if settings.write_energies == True and verbosity > -1:
        write_results(settings.file_results, "{:16} {:8.3f} {:5}\n".format(compound.name,E_0,compound.multiplicity), "a")



def get_cubic_divergence(V,latt_param):
    # calculate cubic optimum
    a_opt = V**(1/3)
    latt_dev = 0.0

    for x in latt_param:
        latt_dev += (x-a_opt)**2
        
    latt_dev = (latt_dev)**(1/2)

    return latt_dev/a_opt

# calculates the RMSD for a given list of parameters and estimators
import sys
def get_rmsd(parameter,estimator,normalize=""):
    
    if not len(parameter) == len(estimator):
        print("Error in evaluation of RMSD. Parameter and estimator aren't of equal length.")
        sys.exit()
    
    # RMSD calculation plus average and delta of the parameters for normalization to NRMSD (if desired, check wikipedia article for RMSD for details)
    rmsd = 0
    av_parameter = 0
    for i in range(len(parameter)):
        rmsd += (estimator[i]-parameter[i])**2
        av_parameter += parameter[i]/len(parameter)
    
    rmsd *= 1/len(parameter)
    rmsd = rmsd**0.5
    
    delta_parameter = max(parameter)-min(parameter)
    
    if "av" in normalize.lower():
        rmsd *= 1/av_parameter
    elif "de" in normalize.lower():
        rmsd *= 1/delta_parameter
    
    return rmsd

# calculates the distance between to n-dimensional vectors
def get_distance(pos1, pos2):
    
    if not len(pos1) == len(pos2):
        print("Error in evaluation of distances. Vectors aren't of equal length.")
        sys.exit()
    
    distance = 0
    for i in range(len(pos1)):
        distance += (pos1[i]-pos2[i])**2
    distance = np.sqrt(distance)
        
    return distance
    





















