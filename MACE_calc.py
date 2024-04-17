#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 13:23:50 2024

@author: bt308570
"""
  
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
   
# function for single point calculation
def perform_singlepoint(compound):
    E_0 = compound.structure.get_potential_energy()
    n_atoms = len(compound.structure.symbols)
    
    print("Energy: {:9.3f} eV, {:6.3f} eV/atom".format(E_0,E_0/n_atoms))
    # write_results(file, , mode = "a")
    write_results(file_results, "{:16} {:8.3f} {:5}\n".format(compound.name,E_0,compound.multiplicity), "a")
    print("")
    
    
from ase.constraints import UnitCellFilter
from ase.optimize import BFGSLineSearch as BFGS_LS

def perform_optimization(compound):
    
    # relaxation
    E_0 = compound.structure.get_potential_energy()
    n_atoms = len(compound.structure.symbols)
    print("Initial Energy: {:9.3f} eV, {:6.3f} eV/atom".format(E_0,E_0/n_atoms))
    
    struc_sym = get_symmetry(compound.structure)
    struc_sym.calc = settings.set_calculator() 
    ucf = UnitCellFilter(struc_sym)
    relax = BFGS_LS(atoms=ucf)
    relax.run(fmax=settings.fmax,steps=settings.max_steps)
    
    E_fin = struc_sym.get_potential_energy()
    print("Final Energy: {:9.3f} eV, {:6.3f} eV/atom".format(E_fin,E_fin/n_atoms))
    
    write_results(file_results, "{:16} {:8.3f} {:5}\n".format(compound.name,E_0,compound.multiplicity), "a")

    # struc_sym.write(name+"_MACEOPT.xyz",format="extxyz")

        
        
def write_results(file,parameters,mode):
    with open(file, mode=mode) as f:
        f.write(parameters)


def print_separator(symbol,n=32):
    print("")
    for i in range(n):
        print(symbol,end="")
    print("")
    print("")
    
class time_tracker:
    def __init__(self):
        self.list_times = [("init",time.time())]
           
    def time_evaluation(self,label,mode):
        
        self.list_times.append((label,time.time()))
        
        if self.list_times[-1][0] == "setup":
            clean_label = "Setup Time:"
        elif "compound" in self.list_times[-1][0]:
            number = int(self.list_times[-1][0].split("_")[1])
            
            if number == 1:
                self.start_calculation = self.list_times[-2][1]
    
            clean_label = "Compound:"
        else:
            clean_label = "MISSING"
        
        if mode == "step":
            print("{:<18} {:8.3f} s".format(clean_label,self.list_times[-1][1]-self.list_times[-2][1]))
        elif mode == "total":
            print("Total Calculation: {:8.3f} s".format(self.list_times[-1][1]-self.list_times[0][1]))
  
    def time_estimator(self,n_compounds):
       number = int(time_tracker.list_times[-1][0].split("_")[1])
       time_current = time_tracker.list_times[-1][1]
       elapsed_time = time_current-time_tracker.start_calculation
       average = (time_current-time_tracker.start_calculation)/number
       estimated_time = average*n_compounds
        
       return elapsed_time, estimated_time-elapsed_time, average
   
   
import sys, os, time
# beginning of main program
if __name__ == "__main__":
    
    
    time_tracker = time_tracker()
    
    # fetch a link to the modules
    if not os.environ["PYTHON_MODULES"] in sys.path:
        sys.path.append(os.environ["PYTHON_MODULES"])
    # print(sys.path)
    
    import MACE_calc_setup as setup
    
    # handle SLURM job - separate file
    # DONE      read in arguments with basic data
    # PENDING   read in settings file with more complex data
    # DONE      set up model for calculation
    # DONE      handle basic single-point calculation
    # DONE      handle optimizations
    # DONE      track time of calculation steps
    # PENDING   implement graceful exit
    # PENDING   implement actually helpful help message    
    
    
    # PENDING   handle phonon calculations
    # PENDING   handle MD simulations
    
    
    # read in arguments with basic data
    
    compounds, settings, procedure = setup.read_arguments(sys.argv)
    calculator = settings.set_calculator()       
    
    print_separator("-")
    
    time_tracker.time_evaluation(label="setup",mode="step")
    
    print_separator("=")    
    # create output file
    file_results = "MACE_"+procedure+".dat"
    write_results(file_results, "{:16} {:8} {:5}\n".format("ID","E","Multi"), "w") 
   
    i = 0
    for compound in compounds:
        i += 1
        compound.structure.calc = calculator
        
        print("Compound: {} {}".format(compound.name,i))
        
        # basic structure analysis module
        if procedure == "analysis":
            print("Test")
            sys.exit()
        
        if settings.verbosity > 0:
            get_symmetry(compound.structure,settings.verbosity)
        
        # single-point module
        if procedure == "sp":
            perform_singlepoint(compound)
            
        # optimization module
        if procedure == "opt":
            perform_optimization(compound)
            
        time_tracker.time_evaluation("compound_"+str(i),"step")
        
               
        elapsed, left, average = time_tracker.time_estimator(len(compounds))
        print("Time elapsed {:8.3f} s, Time left: {:8.3f} s".format(elapsed, left))
        
        if elapsed + left > settings.max_time:
            print("Estimated runtime likely exceeds time limit of {} s.".format(settings.max_time))
        if elapsed > settings.max_time - 600:
            print("RUNNING OUT OF TIME!")
            # trigger graceful exit instead !!!
        
        print_separator("-")
     
        
print_separator("=")
elapsed, left, average = time_tracker.time_estimator(len(compounds))
print("Average time per calculation: {:8.3f} s".format(average))
time_tracker.time_evaluation("final","total")
    
    
    
    














