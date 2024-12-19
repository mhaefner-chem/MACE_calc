#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 13:23:50 2024

@author: bt308570
"""
  
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
            print("Total Time Calculation: {:8.3f} s".format(self.list_times[-1][1]-self.list_times[0][1]))
  
    def time_estimator(self,n_compounds):
       number = int(time_tracker.list_times[-1][0].split("_")[1])
       time_current = time_tracker.list_times[-1][1]
       elapsed_time = time_current-time_tracker.start_calculation
       average = (time_current-time_tracker.start_calculation)/number
       estimated_time = average*n_compounds
        
       return elapsed_time, estimated_time-elapsed_time, average
   
   
import sys, os, time, datetime
# beginning of main program
if __name__ == "__main__":
    
    
    time_tracker = time_tracker()
    
    # fetch a link to the modules
    if not os.environ["PYTHON_MODULES"] in sys.path:
        sys.path.append(os.environ["PYTHON_MODULES"])
    # print(sys.path)
    
    import MACE_calc_setup as setup
    import MACE_util as util
    import MACE_calc_proc as proc
    
    # handle SLURM job - separate file
    # DONE      read in arguments with basic data
    # WIP       read in settings file with more complex data
    # DONE      set up model for calculation
    # DONE      handle basic single-point calculation
    # DONE      handle optimizations
    # DONE      track time of calculation steps
    # DONE      implement graceful exit
    # DONE      implement actually helpful help message    
    
    
    # WIP       handle phonon calculations
    # WIP       handle MD simulations
    
    
    # read in arguments with basic data
    compounds, settings, procedure = setup.read_arguments(sys.argv)
    calculator = settings.set_calculator()       
    
    util.print_separator("-")
    
    settings.print_all()
    
    util.print_separator("-")
    
    time_tracker.time_evaluation(label="setup",mode="step")
    
    util.print_separator("=")    
    # create output file
    settings.file_results = "MACE_"+procedure+".dat"
    if settings.write_energies == True and settings.verbosity > -1:
        util.write_results(settings.file_results, "{:16} {:8} {:5} {}\n".format("ID","E","Multi",datetime.datetime.now()), "a") 
   
    pbc_numerical = []
    for letter in settings.pbc:
       if letter == "0" or letter == "1":
           pbc_numerical.append(int(letter))
    if len(pbc_numerical) != 3:
        print("ERROR in PBC settings!")
   
    i = 0
    for compound in compounds:
        
        
        N = [1,1,1]
        if "x" in settings.sc:
            from ase.build import make_supercell
            tmp = settings.sc.split("x")
            for j in range(3):
                N[j] = int(tmp[j])
            compound.structure = make_supercell(compound.structure,((N[0],0,0),(0,N[1],0),(0,0,N[2])))
        
        i += 1
        compound.structure.set_pbc(pbc_numerical)
        
        compound.structure.calc = calculator
        
        
        print("Compound: {} {}/{}".format(compound.name,i,len(compounds)))
        
        # basic structure analysis module
        if procedure == "analysis":
            print("Test")
            sys.exit()
        
        if settings.verbosity > 0:
            util.get_symmetry(compound.structure,settings.verbosity)
        
        # single-point module
        if procedure == "sp":
            proc.perform_singlepoint(compound,settings)
            
        # optimization module
        if procedure in ["opt"]:
            if not os.path.isfile(compound.name+"_OPT.xyz"):
                proc.perform_optimization(compound,settings)
            else:
                print("Calculation of "+compound.name+"_OPT.xyz already performed.")
            
        # optimization module
        if procedure == "bulk":
            if not os.path.isfile(compound.name+"_OPT.xyz"):
                proc.perform_optimization(compound,settings)
            else:
                print("Using structure from previous optimiziation {}_OPT.xyz".format(compound.name))
            proc.perform_bulkmod(compound,settings)
            
        if procedure == "phon":
            if not os.path.isfile(compound.name+"_OPT.xyz"):
                proc.perform_optimization(compound,settings)
            else:
                print("Using structure from previous optimiziation {}_OPT.xyz".format(compound.name))
            proc.perform_phonon(compound,settings)
        
        if procedure == "md":
            if settings.md_algo in ["nptb","npt"]:
                if settings.set_bulk_mod < 0:
                    proc.perform_optimization(compound,settings)
                    proc.perform_bulkmod(compound,settings)
                else:
                    compound.bulk_mod = settings.set_bulk_mod * util.unit_conversion("p", "GPa", "eV")
            if len(compounds) > 1:
                print("More than one compound specified for MD. Only one will be performed.")
            proc.perform_md(compound,settings)
        
        if procedure == "neb":
            proc.perform_neb(compound,settings)
            
        time_tracker.time_evaluation("compound_"+str(i),"step")
        
        if len(compounds) > 1:
            elapsed, left, average = time_tracker.time_estimator(len(compounds))
            print("Time elapsed {:8.3f} s, Time left: {:8.3f} s".format(elapsed, left))
            
            if elapsed + left > settings.max_time:
                print("Estimated runtime likely exceeds time limit of {} s.".format(settings.max_time))
            if elapsed > settings.max_time - 60:
                print("Approaching runtime limit, initiating graceful exit.")
                with open("compounds_remaining.dat", mode="w") as f:
                    for j in range(i,len(compounds)):
                        f.write(compounds[j].file+"\n")
                        
                    util.print_separator("=")
                    elapsed, left, average = time_tracker.time_estimator(len(compounds))
                    print("Average time per calculation: {:8.3f} s".format(average))
                    time_tracker.time_evaluation("final","total")
                sys.exit()
                    
            if i != len(compounds):
                util.print_separator("-")
     
        
util.print_separator("=")
elapsed, left, average = time_tracker.time_estimator(len(compounds))
print("Average time per calculation: {:8.3f} s".format(average))
time_tracker.time_evaluation("final","total")
    
    
    
    














