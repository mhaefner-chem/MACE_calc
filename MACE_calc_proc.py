#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 09:43:58 2024

@author: bt308570

MACE procedures
"""
import MACE_util as util

# function for single point calculation
def perform_singlepoint(compound,settings):
    compound.structure.get_potential_energy()
    print("")
    util.print_results(settings,compound,compound.structure,settings.verbosity)
    print("")
    
from ase.constraints import UnitCellFilter
from ase.optimize import BFGSLineSearch as BFGS_LS

def perform_optimization(compound,settings):
    
    # relaxation
    compound.structure.get_potential_energy()
    util.print_results(settings,compound,compound.structure,verbosity=-1)
    
    struc_sym = util.get_symmetry(compound.structure)
    struc_sym.calc = settings.calculator
    if settings.sym == True:
        ucf = UnitCellFilter(struc_sym)
    else:
        ucf = struc_sym
    print("Optimization with F_max = {} eV/Å".format(settings.f_max))
    print("")
    relax = BFGS_LS(atoms=ucf)
    relax.run(fmax=settings.f_max,steps=settings.max_steps)
    
    print("Remaining pressure: {:9.3f} kbar".format(-np.mean(struc_sym.get_stress()[0:3])* util.unit_conversion("p", "eV", "kbar")))
    
    # relax.run(fmax=0.1)
    
    print("")
    util.print_results(settings,compound,struc_sym,settings.verbosity)
    
    # compound.structure = struc_sym
    if settings.write_geometries == True:
        struc_sym.write(compound.name+"_OPT.xyz",format="extxyz")
    

# function for calculation of bulk modulus

def perform_bulkmod(compound,settings):
    opt_struc = ase_read(compound.name+"_OPT.xyz",format="extxyz")
    
    strucs = {"big":{"struc":opt_struc.copy(),"E":0,"V":0,"p":0},"small":{"struc":opt_struc.copy(),"E":0,"V":0,"p":0}}

    
    
    for size in strucs:
        
        util.print_separator("-")
        if size == "big":
            factor = 1
            label = "Bigger"
        elif size == "small":
            factor = -1
            label = "Smaller"
        else:
            factor = 0
                
        strucs[size]["struc"].cell = strucs[size]["struc"].cell * (1 + factor * settings.bulk_mod_delta)
        strucs[size]["struc"].calc = settings.calculator
        
        
        relax = BFGS_LS(atoms=strucs[size]["struc"])
        relax.run(fmax=settings.f_max,steps=settings.max_steps)
        strucs[size]["E"] = strucs[size]["struc"].get_potential_energy()
        strucs[size]["V"] = strucs[size]["struc"].get_volume()
        strucs[size]["p"] = np.mean(strucs[size]["struc"].get_stress()[0:3])
        
        print("{} cell".format(label))
        print("Volume: {:9.3f} Å³, pressure: {:6.3f} GPa".format(strucs[size]["V"],-strucs[size]["p"]* util.unit_conversion("p", "eV", "GPa")))
    
    util.print_separator("-")
    
    delta_V = strucs["big"]["V"]-strucs["small"]["V"]
    delta_p = strucs["big"]["p"]-strucs["small"]["p"]
    
    bulk_mod = (strucs["big"]["V"]+strucs["small"]["V"])/2 * delta_p/delta_V
    
    print("Bulk modulus: {:6.3f} GPa".format(bulk_mod * util.unit_conversion("p", "eV", "GPa")))
    compound.bulk_mod = bulk_mod
    

# function for phonon calculation
from ase.phonons import Phonons
from ase.thermochemistry import CrystalThermo
from ase.io import read as ase_read
import numpy as np
import time
def perform_phonon(compound,settings):
    temperatures = settings.temperatures
    opt_struc = ase_read(compound.name+"_OPT.xyz",format="extxyz")
    phonon_start = time.time()
    if settings.f_max > 0.005:
        print("Fmax for the optimization is larger than 0.005 eV/Å. Consider a tighter convergence for phonon calculations.")
    
    N = [1,1,1]
    min_len = settings.min_len 

    for i in range(len(opt_struc.cell.lengths())):
        while N[i]*opt_struc.cell.lengths()[i] < min_len:
            N[i] += 1
        

    print(len(opt_struc.symbols))
    print(opt_struc.cell.lengths())
    displacements = len(opt_struc.symbols)*6*N[0]*N[1]*N[2]
    print(N,displacements)

    # Phonon analysis
    ph = Phonons(opt_struc, settings.calculator, supercell=(N[0], N[1], N[2]), delta=0.05)
    ph.run()
    print("Phonon calculation done.")

    displacement_end = time.time()
    
    ph.read(acoustic=True)
    phonon_energies, phonon_DOS = ph.dos(kpts=(20, 20, 20), npts=3000,
                                         delta=5e-4)
    t_check_start = time.time()
    path = []
    spread = 7
    for i in np.linspace(-0.5,0.5,spread):
        for j in np.linspace(-0.5,0.5,spread):
            for k in np.linspace(-0.5,0.5,spread):
                path.append((i,j,k))
                
    bs = ph.band_structure(path,verbose=False)
    flip = False
    for i in range(len(bs)):
       for j in range(len(bs[0])): 
           if bs[i][j] < 0:
               if flip == False:
                   print("Imaginary frequencies")
                   print("{:>19} {:>9} {:>7}".format("k-point","eV","1/cm"))
                   flip = True
               print("({:5.2f},{:5.2f},{:5.2f}) {:9.6f} {:7.2f}"
                     .format(path[i][0],path[i][1],path[i][2],bs[i][j],bs[i][j]*8065.54))


    t_check_end = time.time()
    print()

    # Calculate the Helmholtz free energy
    thermo = CrystalThermo(phonon_energies=phonon_energies,
                           phonon_DOS=phonon_DOS,
                           potentialenergy=compound.e_opt,
                           formula_units=1)
    for T in temperatures:
        print("===============================")
        print("Thermodynamics at {:9.2f} K".format(T))
        print("===============================")
        F = thermo.get_helmholtz_energy(temperature=T)


    print("Displacements/s: {:9.3f}".format(displacement_end-phonon_start))
    print("Phonons/s: {:9.3f}".format(t_check_start-displacement_end))
    print("Time for imaginary check: {:9.3f} s".format(t_check_end-t_check_start))


from ase.build import make_supercell
from ase.md.langevin import Langevin
from ase.md.nptberendsen import NPTBerendsen
from ase import units
from ase.md import MDLogger
import sys
def perform_md(compound,settings):
    
    N = [1,1,1]
    min_len = 3

    for i in range(len(compound.structure.cell.lengths())):
        while N[i]*compound.structure.cell.lengths()[i] < min_len:
            N[i] += 1

    # sys.stdout = open(name+'.out', 'a')

    print(len(compound.structure.symbols),len(compound.structure.symbols)*N[0]*N[1]*N[2])
    print(compound.structure.cell.lengths())


    sc = make_supercell(compound.structure,((N[0],0,0),(0,N[1],0),(0,0,N[2])))
    print(sc.cell.lengths())
    sc.set_calculator(settings.calculator)

    MD_init = time.time()

    # dyn = Langevin(sc, settings.step_size*units.fs, temperature_K=settings.md_T, friction=5e-3)
    dyn = NPTBerendsen(sc, timestep=settings.step_size*units.fs, temperature_K=settings.md_T,
                   taut=100 * units.fs, pressure_au=1.01325 * units.bar,
                   taup=1000 * units.fs, compressibility_au=4.57e-5 / units.bar)
    
    def write_frame():
        MD_step = time.time()
        dyn.atoms.write(compound.name+"_MD.xyz",format="extxyz",append=True)
        print(dyn.nsteps,dyn.max_steps,dyn.nsteps/(MD_step-MD_init))
        with open(compound.name+"_MD.out",mode="a") as f:
            f.write("{} {} {}\n".format(dyn.nsteps,dyn.max_steps,dyn.nsteps/(MD_step-MD_init)))
    dyn.attach(write_frame, interval=30)
    dyn.attach(MDLogger(dyn, sc, 'md.log', header=False, stress=False,
               peratom=True, mode="a"), interval=10)
    dyn.run(1000)
    print("MD finished!")
    
    sys.exit()
