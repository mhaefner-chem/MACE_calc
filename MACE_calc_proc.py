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
    

    struc = util.get_symmetry(compound.structure)

    struc.calc = settings.calculator
    if settings.sym == True:
        ucf = UnitCellFilter(struc)
    else:
        ucf = compound.structure
    print("Optimization with F_max = {} eV/Å".format(settings.opt_f_max))
    print("")
    relax = BFGS_LS(atoms=ucf)
    relax.run(fmax=settings.opt_f_max,steps=settings.opt_max_steps)
    
    print("Remaining pressure: {:9.3f} kbar".format(-np.mean(struc.get_stress()[0:3])* util.unit_conversion("p", "eV", "kbar")))
    
    # relax.run(fmax=0.1)
    
    print("")
    util.print_results(settings,compound,struc,settings.verbosity)
    
    # compound.structure = struc_sym
    if settings.write_geometries == True:
        struc.write(compound.name+"_OPT.xyz",format="extxyz")
    

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
        relax.run(fmax=settings.opt_f_max,steps=settings.opt_max_steps)
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
    if settings.opt_f_max > 0.005:
        print("Fmax for the optimization is larger than 0.005 eV/Å. Consider a tighter convergence for phonon calculations.")
    
    N = [1,1,1]
    min_len = settings.phon_min_len 

    for i in range(len(opt_struc.cell.lengths())):
        while N[i]*opt_struc.cell.lengths()[i] < min_len:
            N[i] += 1
    
    util.print_separator("-")
    print("Supercell size:")
    print("{}x{}x{}".format(N[0],N[1],N[2]))
    
    abc = opt_struc.cell.lengths()
    for i in range(3):
        abc[i] = abc[i] * N[i]
    displacements = len(opt_struc.symbols)*6*N[0]*N[1]*N[2]
    
    print("{:6.3f} {:6.3f} {:6.3f}".format(abc[0],abc[1],abc[2]))
    print()
    print("Atoms unit cell: {:4}".format(len(opt_struc.symbols)))
    print("Atoms supercell: {:4}".format(len(opt_struc.symbols)*N[0]*N[1]*N[2]))
    print("Displacements: {}".format(displacements))
    print()

    util.print_separator("-")
    

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
           if bs[i][j] < -1E-4:
               if flip == False:
                   print("Imaginary frequencies")
                   print("{:>19} {:>9} {:>7}".format("k-point","eV","1/cm"))
                   flip = True
               print("({:5.2f},{:5.2f},{:5.2f}) {:9.6f} {:7.2f}"
                     .format(path[i][0],path[i][1],path[i][2],bs[i][j],bs[i][j]*8065.54))
           elif bs[i][j] < 0:
               print("No significant imaginary frequencies found!")


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
from ase.io import read
from ase import units
import datetime
from ase.md import MDLogger
import sys, math, os
def perform_md(compound,settings):
    
    N = [1,1,1]
    
    if "x" in settings.md_sc:
        tmp = settings.md_sc.split("x")
        for i in range(3):
            N[i] = int(tmp[i])
    else:
        for i in range(len(compound.structure.cell.lengths())):
            while N[i]*compound.structure.cell.lengths()[i] < settings.md_min_len:
                N[i] += 1
    
    

    # sys.stdout = open(name+'.out', 'a')

    util.print_separator("-")
    abc = compound.structure.cell.lengths()
    print("Original unit cell:")
    print("{} atoms".format(len(compound.structure.symbols)))
    print("a {:6.3f} Å, b {:6.3f} Å, c {:6.3f} Å".format(abc[0],abc[1],abc[2]))
    print()
    


    sc = make_supercell(compound.structure,((N[0],0,0),(0,N[1],0),(0,0,N[2])))
    
    
    abc_sc = sc.cell.lengths()
    print("{}x{}x{} supercell:".format(N[0],N[1],N[2]))
    print("{} atoms".format(len(sc.symbols)))
    print("a {:6.3f} Å, b {:6.3f} Å, c {:6.3f} Å".format(abc_sc[0],abc_sc[1],abc_sc[2]))
    util.print_separator("-")
    
    if settings.md_interval_write_s != settings.md_interval_write_e:
        print('''
INFO! MD times for energies/temperatures and structures will not match.
Consider changing MD_interval_s (for structures) and MD_interval_e (for energies) to equal values for more simple plotting.''')

    sc.set_calculator(settings.calculator)

    
    
    # adapt settings to ramping temperature or constant temperature
    if settings.md_T_min > 0 or settings.md_T_max > 0:
        # subdivide by md_T_steps
        temperatures = np.linspace(settings.md_T_min, settings.md_T_max, num=settings.md_T_steps)
        steps = math.floor(settings.md_step_n/settings.md_T_steps)
        print("Using temperature ramp with {} steps per MD and following temperatures:".format(steps))
        for temperature in temperatures:
            print("{:6.2f} K".format(temperature))
    else:
        temperatures = [settings.md_T]
        steps = settings.md_step_n
        
    md_i = 0
    for temperature in temperatures: 
        MD_init = time.time()
        md_i += 1
        # check whether previous MD run exists, and use last step as initial structure for MD
        if os.path.isfile(compound.name+"_MD.xyz"):
            read(compound.name+"_MD.xyz",index=-1,format="extxyz")
            print("Read last structure from previous MD run.")
        
        # selects the thermostat for the dynamics
        if settings.md_algo == "nptb":
            from ase.md.nptberendsen import NPTBerendsen
            dyn = NPTBerendsen(sc, timestep=settings.md_step_size*units.fs, temperature_K=temperature,
                           taut=settings.md_taut * units.fs, pressure_au=settings.md_p * units.bar,
                           taup=settings.md_taup * units.fs, compressibility_au=1/compound.bulk_mod)
        elif settings.md_algo == "npt":
            from ase.md.npt import NPT
            dyn = NPT(sc, timestep=settings.md_step_size*units.fs,temperature_K=temperature,
                      externalstress=settings.md_p * units.bar,ttime=settings.md_ttime * units.fs,
                      pfactor=settings.md_ptime**2 * compound.bulk_mod * units.fs**2)
            print("Pfactor:",settings.md_ptime**2 * compound.bulk_mod/units.GPa,"GPa fs²")
        elif settings.md_algo == "nve":
            from ase.md.verlet import VelocityVerlet
            dyn = VelocityVerlet(sc, settings.md_step_size * units.fs)
        elif settings.md_algo == "nvt":
            from ase.md.langevin import Langevin
            dyn = Langevin(sc, settings.md_step_size*units.fs, temperature_K=temperature, friction=5e-3)
        
        # creates output with timings
        with open(compound.name+"_MD.out",mode="a") as f:
            header = "{:12} {:12} {:12} {:12} {:12} {:12} {:12} {} T: {} K\n".format("sim. time","curr. step","tot. steps","time/step","time ela.","time rem.","time tot.",datetime.datetime.now(),temperature)
            f.write(header)
            header_units = "{:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}\n".format("ps","","","s","s","s","s")
            f.write(header_units)
        
        # writes relevant timing info and structures
        def write_frame():
            MD_step = time.time()
            dyn.atoms.write(compound.name+"_MD.xyz",format="extxyz",append=True)
            
            # if dyn.nsteps == 0:
                # print(header,end="")
                # continue
            if dyn.nsteps != 0:
                time_step = (MD_step-MD_init)/dyn.nsteps
                rem_steps = settings.md_step_n - dyn.nsteps
                
                message = "{:12.3f} {:12} {:12} {:12.3f} {:12.3f} {:12.3f} {:12.3f}\n".format(dyn.nsteps*settings.md_step_size/1000,dyn.nsteps,steps,time_step,MD_step-MD_init,rem_steps * time_step, steps * time_step)
                # print(message,end="")
                with open(compound.name+"_MD.out",mode="a") as f:
                    f.write(message)
                    
        dyn.attach(write_frame, interval=settings.md_interval_write_s)
        dyn.attach(MDLogger(dyn, sc, compound.name+"_MD.log", header=True, stress=False,
                    peratom=True, mode="a"), interval=settings.md_interval_write_e)
        dyn.run(steps)
        if len(temperatures) == 1:
            print("MD finished!")
        else:
            print("MD {}/{} finished.".format(md_i,len(temperatures)))
        
    sys.exit()
