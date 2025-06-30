#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 13:22:35 2024

@author: bt308570

# input processing for MACE_calc

"""

import sys,os
from ase.io import read as ase_read
import MACE_calc_compound as compound
import numpy as np

class settings:
    def __init__(self):
        # I/O
        self.filename = "MACE_calc.inp"
        self.name = os.path.splitext(os.path.basename(self.filename))
        self.procedure = "analysis"
        self.file_results = "out.out"
        
        self.verbosity = 1
        self.write_geometries = True
        self.print_geometries = False
        self.write_energies = True
        self.print_forces = False
        self.print_stresses = False
        
        self.pbc = "111"
        
        
        # calculator parameters
        self.model = "macemp0"
        self.acc = "float64"
        self.dispersion = 0
        self.use_gpu = False
        
        # optimization
        self.opt_f_max = 0.01
        self.opt_max_steps = 250
    
        self.sym = True
        self.opt_cell = True
        self.opt_mask = [[1,1,1],[1,1,1],[1,1,1]]
        
        # bulk modulus calculation
        self.bulk_mod_delta = 0.003322
        self.set_bulk_mod = -1
        
        # phonons
        self.phon_min_len = 12.0
        self.phon_t_min = 300
        self.phon_t_max = 300
        self.phon_t_steps = 1
        self.temperatures = np.linspace(self.phon_t_min,self.phon_t_max,self.phon_t_steps)
        
        # MD
        self.md_T = 300
        self.md_T_start = -1
        self.md_T_end = -1
        self.md_T_steps = 3
        
        self.md_p = 1.01325
        
        self.md_step_size = 1.5
        self.md_step_n = 100
        
        self.md_min_len = 12.0
        self.sc = "" # explicitly set supercell size
        self.md_taut = 100
        self.md_taup = 1000
        self.md_algo = "nve"
        
        # NVT Nose-Hoover
        self.md_tdamp = 100
        self.md_tchain = 3
        self.md_tloop = 1

        # NPT Nose-Hoover
        self.md_ttime = 25.0
        self.md_ptime = 75.0
        
        self.md_interval_write_e = 10
        self.md_interval_write_s = 10
        
        # neb management
        self.neb_init_file = "01.xyz"
        self.neb_end_file = "02.xyz"
        self.neb_n_images = 5
        self.neb_f_max = 0.05
        self.neb_climb = True
        
        # time management
        self.max_time = 72000 # 72000 s/20 h as time limit
        
        
        
    def read_file(self,filename):
        # handle the settings file
        with open(filename,mode="r") as f:
            lines = f.readlines()
            for line in lines:
                print(line)
                if "#" in line.lower():
                    continue
                elif "dispersion" in line.lower():
                    self.dispersion = int(line.split()[-1])
                elif "verbosity" in line.lower():
                    self.verbosity = int(line.split()[-1])
                elif "model" in line.lower():
                    self.model = line.split()[-1]
                elif "acc" in line.lower():
                    self.acc = line.split()[-1]
                elif "opt_f_max" in line.lower():
                    self.opt_f_max = float(line.split()[-1])
                elif "max_steps" in line.lower():
                    self.opt_max_steps = int(line.split()[-1])
                elif "max_time" in line.lower():
                    self.max_time = float(line.split()[-1])
                elif "proc" in line.lower():
                    self.procedure = line.split()[-1]
                elif "pbc" in line.lower():
                    self.pbc = line.split()[-1]
                # bulk modulus settings
                elif "delta_bulk" in line.lower():
                    self.bulk_mod_delta = float(line.split()[-1])
                elif "set_bulk_mod" in line.lower():
                    self.set_bulk_mod = float(line.split()[-1])
                    
                
                
                
                # settings for phonon calculations
                elif "phon_min_len" in line.lower():
                    self.phon_min_len = float(line.split()[-1])
                elif "phon_t_min" in line.lower():
                    self.phon_t_min = float(line.split()[-1])
                elif "phon_t_max" in line.lower():
                    self.phon_t_max = float(line.split()[-1])
                elif "phon_t_steps" in line.lower():
                    self.phon_t_steps = int(line.split()[-1])
                    
                # settings for MD calculations
                elif "md_ttime" in line.lower():
                    self.md_ttime = float(line.split()[-1])
                elif "md_ptime" in line.lower():
                    self.md_ptime = float(line.split()[-1])
                elif "md_t_start" in line.lower():
                    self.md_T_start = float(line.split()[-1])
                elif "md_t_end" in line.lower():
                    self.md_T_end = float(line.split()[-1])
                elif "md_t_steps" in line.lower():
                    self.md_T_steps = int(line.split()[-1])  
                elif "md_p" in line.lower():
                    self.md_p = float(line.split()[-1])
                elif "md_step_size" in line.lower():
                    self.md_step_size = float(line.split()[-1])
                elif "md_min_len" in line.lower():
                    self.md_min_len = float(line.split()[-1])
                elif "sc" in line.lower():
                    self.sc = line.split()[-1]
                elif "md_step_n" in line.lower():
                    self.md_step_n = int(line.split()[-1])
                elif "md_interval_write_e" in line.lower():
                    self.md_interval_write_e = int(line.split()[-1])
                elif "md_interval_write_s" in line.lower():
                    self.md_interval_write_s = int(line.split()[-1])
                elif "md_taut" in line.lower():
                    self.md_taut = float(line.split()[-1])
                elif "md_taup" in line.lower():
                    self.md_taup = float(line.split()[-1])
                elif "md_tdamp" in line.lower():
                    self.md_tdamp = float(line.split()[-1])
                elif "md_tchain" in line.lower():
                    self.md_tchain = int(line.split()[-1])
                elif "md_tloop" in line.lower():
                    self.md_tloop = int(line.split()[-1])
                elif "md_algo" in line.lower():
                    self.md_algo = line.split()[-1]
                elif "md_t" in line.lower():
                    self.md_T = float(line.split()[-1])
                    
                    
                # settings for NEB
                elif "neb_init_file" in line.lower():
                    self.neb_init_file = line.split()[-1]
                elif "neb_end_file" in line.lower():
                    self.neb_end_file = line.split()[-1]
                elif "neb_climb" in line.lower():
                    if line.split()[-1].lower() == "true":
                        self.neb_climb = True
                    else:
                        self.neb_climb = False
                elif "neb_n_im" in line.lower():
                    self.neb_n_images = int(line.split()[-1])
                elif "neb_f" in line.lower():
                     print("TEST")
                     self.neb_f_max = float(line.split()[-1])
                
                    
                
                    
                
                # writes and prints
                elif "write_geom" in line.lower():
                    if line.split()[-1].lower() == "true":
                        self.write_geometries = True
                    else:
                        self.write_geometries = False
                elif "write_e" in line.lower():
                    if line.split()[-1].lower() == "true":
                        self.write_energies = True
                    else:
                        self.write_energies = False
                elif "print_g" in line.lower():
                    if line.split()[-1].lower() == "true":
                        self.print_geometries = True
                    else:
                        self.print_geometries = False
                elif "print_f" in line.lower():
                    if line.split()[-1].lower() == "true":
                        self.print_forces = True
                    else:
                        self.print_forces = False
                elif "print_s" in line.lower():
                    if line.split()[-1].lower() == "true":
                        self.print_stresses = True
                    else:
                        self.print_stresses = False
                elif "use_gpu" in line.lower():
                    if line.split()[-1].lower() == "true":
                        self.use_gpu = True
                    else:
                        self.use_gpu = False
                
                elif "sym" in line.lower():
                    if line.split()[-1].lower() == "true":
                        self.sym = True
                    else:
                        self.sym  = False
                elif "opt_cell" in line.lower():
                    if line.split()[-1].lower() == "true":
                        self.opt_cell = True
                    else:
                        self.opt_cell  = False
                elif "opt_angle" in line.lower():
                    if line.split()[-1].lower() == "true":
                        self.opt_mask = [1,1,1]
                    else:
                        self.opt_mask = [[1,0,0],[0,1,0],[0,0,1]]
            
        try:
            self.temperatures = np.linspace(self.phon_t_min,self.phon_t_max,self.phon_t_steps)
        except:
            self.temperatures = [300]
    
    # constructs the calculator
    def set_calculator(self):
        self.calculator = setup_model(self.model,self.dispersion,self.use_gpu,self.acc)
        return self.calculator

    def print_all(self):
        
        
        def print_item(text,item):
            format = "  {:20} {}"
            print(format.format(text,item))
        
        if self.verbosity > 0:
            print_item("Procedure:",self.procedure)
            print("")
            
            print("I/O")
            print("")
            print_item("File_settings:",self.filename)
            print_item("File_results:",self.file_results)

            print("")
            print_item("Verbosity:",self.verbosity)
            print_item("Write_geom:",self.write_geometries)
            print_item("Write_e:",self.write_energies)
            print_item("Print_f:",self.print_forces)
            print_item("Print_s:",self.print_stresses)
            print_item("PBC:",self.pbc)
            print("")
            
            # calculator parameters
            print("Calculator")
            print("")
            print_item("Model:",self.model)
            print_item("Accuracy:",self.acc)
            print_item("Dispersion (2=D3BJ):",self.dispersion)
            print_item("Use_GPU:",self.use_gpu)
            print_item("Sym:",self.sym)
            if self.set_bulk_mod > 0:
                print_item("set_bulk_mod [GPa]:",self.set_bulk_mod)
            print("")
            
            # optimization
            if self.procedure in ["opt","bulk","phon","analysis"]:
                print("Optimization")
                print("")
                print_item("OPT_F_max [eV/Å]:",self.opt_f_max)
                print_item("OPT_Max_steps:",self.opt_max_steps)
                print_item("OPT_cell:",self.opt_cell)
                print_item("OPT_mask:",self.opt_mask)
                print("")
                
            # bulk modulus calculation
            if self.procedure in ["bulk","analysis"]:
                print("Bulk Modulus")
                print("")
                print_item("BULK_MOD_Delta [%]:",self.bulk_mod_delta)
                print("")
            
            # phonons
            
            if self.procedure in ["phon","analysis"]:
                print("Phonons")
                print("")
                print_item("PHON_Min_len [Å]:",self.phon_min_len)
                print_item("PHON_T_min [K]:",self.phon_t_min)
                print_item("PHON_T_max [K]:",self.phon_t_max)
                print_item("PHON_T_steps:",self.phon_t_steps)
                print("")
            
            # MD
            if self.procedure in ["md","analysis"]:
                print("MD")
                print("")
                print_item("MD_T [K]:",self.md_T)
                if self.md_T_start > 0:
                    print_item("md_T_start [K]:",self.md_T_start)
                if self.md_T_end > 0:
                    print_item("md_T_end [K]:",self.md_T_end)
                    print_item("MD_T_steps [K]:",self.md_T_steps)
                print_item("MD_p [bar]:",self.md_p)
                print_item("MD_Step_size [fs]:",self.md_step_size)
                print_item("MD_step_n:",self.md_step_n)
                print_item("MD_Min_len [Å]:",self.md_min_len)
                print_item("MD_interval_write_E [steps]:",self.md_interval_write_e)
                print_item("MD_interval_write_S [steps]:",self.md_interval_write_s)
                if self.md_algo == "nptb" or self.md_algo == "nvt_bussi":
                    print_item("MD_tauT [1/fs]:",self.md_taut)
                if self.md_algo == 'nvt_nose':
                    print_item("MD_tdamp [fs]:",self.md_tdamp)
                    print_item("MD_tchain:",self.md_tchain)
                    print_item("MD_tloop:",self.md_tloop)
                if self.md_algo == "nptb":
                    print_item("MD_taup [1/fs]:",self.md_taup)
                elif self.md_algo == "npt":
                    print_item("MD_ttime [fs]:",self.md_ttime)
                    print_item("MD_ptime [fs]:",self.md_ptime)
                    print_item("OPT_mask:",self.opt_mask)
                print_item("MD_algo:",self.md_algo)
                print("")
                
            if self.procedure in ["neb","analysis"]:
                print("NEB")
                print("")
                print_item("NEB_init_file:",self.neb_init_file)
                print_item("NEB_end_file:",self.neb_end_file)
                print_item("NEB_climb:",self.neb_climb)
                print_item("NEB_n_images:",self.neb_n_images)
                print_item("NEB_f_max [eV/Å]:",self.neb_f_max)
                print("")
            
            # time management
            print("Time Management")
            print_item("Max_time [s]:",self.max_time)



# reads the arguments and extracts the information
def read_arguments(arguments):
    settings_calculation = settings()    
    compounds = [] 
    
    options = {}
    options["help"] = ["-h","--help"]
    options["compound"] = ["-c","--cell","--compound"]
    options["compound_list"] = ["-l","--list"]
    options["dispersion"] = ["-d","--d3","--dispersion","--d3bj"]
    options["model"] = ["-m","--model"]
    options["procedure"] = ["-p","--proc","--procedure"]
    options["settings"] = ["-s","--setting","--settings","-i","--input"]
    options["verbosity"] = ["-v","--verbosity"]
    options["symmetry"] = ["--nosym"]   
    options["gpu"] = ["-g","--gpu"]
    options["molecule"] = ["--mol","--molecule"]
 
    n = len(sys.argv)
    
    compounds = []

    for i in range(n):
        argument = sys.argv[i].lower()
        
        if argument in options["help"]:
            help_output(options)
        elif argument in options["compound"]:
            for item in read_compound(sys.argv[i+1]):
                compounds.append(item)
        elif argument in options["compound_list"]:
            for item in read_compound_list(sys.argv[i+1]):
                compounds.append(item)
        elif argument in options["dispersion"]:
            settings_calculation.dispersion = 3
        elif argument in options["procedure"]:
            settings_calculation.procedure = sys.argv[i+1]
        elif argument in options["settings"]:
            filename = sys.argv[i+1]
            settings_calculation.name = os.path.splitext(os.path.basename(filename))
        elif argument in options["model"]:
            settings_calculation.model = sys.argv[i+1]
        elif argument in options["verbosity"]:
            settings_calculation.verbosity = int(sys.argv[i+1])
        elif argument in options["symmetry"]:
            settings_calculation.sym = False
        elif argument in options["gpu"]:
            settings_calculation.use_gpu = True
        elif argument in options["molecule"]:
            settings_calculation.pbc = "000"
        
    
    try:
        settings_calculation.read_file(filename)
        print("Processing the input file {}.".format(filename))
    except:
        print("No valid settings/input file found. Proceeding with default settings.")
        
    
    procedures = {}
    procedures["analysis"] = ["analysis","an","a"]
    procedures["sp"] = ["sp","single","singlepoint","s"]
    procedures["opt"] = ["opt","optimization","o"]
    procedures["phon"] = ["p","phon","phonon"]
    procedures["md"] = ["m","md"]
    procedures["bulk"] = ["b","bulk","bulkmodulus"]
    procedures["neb"] = ["n","neb"]
    
    # print(settings_calculation.procedure)
    
    for key in procedures:
        if settings_calculation.procedure in procedures[key]:
            procedure = key
            settings_calculation.procedure = key
            
    return compounds, settings_calculation, procedure


# prints the help dialogue
import MACE_util as util
def help_output(options):
    util.print_separator("-",64,skip=False)
    print('''
This program coordinates single-point calculations, optimizations, phonon calculations, and MD simulations of a structure or set of structures with a MACE model.
''',end="")
    util.print_separator("-",64,skip=False)
    print("Basic usage:")
    print('''MACE_calc.py -c <COMPOUND> -s <SETTINGS_FILE> <OPTIONS>
MACE_calc.py -l <LIST_OF_COMPOUNDS -s <SETTINGS_FILE> <OPTIONS>''')
    util.print_separator("-",64,skip=False)
    print("Options:")
    formatting = "  {:<28} {:<32}"
    print(formatting.format("-h [ --help ]","Prints this help message."))
    print(formatting.format("-c [ --compound ] <ARG>","Input structure <ARG> formatted as CIF, VASP, or (EXT)XYZ."))
    print(formatting.format("-l [ --list ] <ARG>","<ARG> is a list of valid input structures."))
    print(formatting.format("-m [ --model ] <ARG>","<ARG> specifies the MACE model used for the calculation. Default = macemp0."))
    print("")
    print(formatting.format("-p [ --proc ] <ARG>","<ARG> specifies the mode that will be used."))
    print("  Modes:")
    print("    sp: single point")
    print("    opt: optimization")
    print("    phon: phonon calculation")
    print("    md: MD simulation")
    print("")
    print(formatting.format("-s [ --settings ] <ARG>","<ARG> contains more detailed settings for the calculation."))
    print(formatting.format("-d [ --d3bj ]","Activates Grimme D3(BJ) dispersion correction."))
    print(formatting.format("-v [ --verbosity ] <ARG>","Changes the verbosity of the output. Default = 1."))
    print(formatting.format("-g [ --gpu ] <ARG>","Use a GPU for calculations instead of CPU. Default = False."))

    sys.exit()

# reads a structure from a file
def read_compound(file):
    
    atoms = util.read_structure(file)  
        
    if len(atoms) > 1:
        print("{} structures found in input.".format(len(atoms)))
        
    compounds = []
    for i in range(len(atoms)):
        compounds.append(compound.compound(file,atoms[i]))
    return compounds

# processes a given list of files
def read_compound_list(file):
    
    with open(file,mode="r") as f:
        compounds = []
        for line in f.readlines():
            for compound in read_compound(line.rstrip("\n")):
                compounds.append(compound)
    return compounds

# initializes the calculator
from mace.calculators import mace_mp, MACECalculator
from ase.calculators.mixing import SumCalculator as Calc_Sum

def setup_model(file,dispersion,gpu,default_dtype="float64"):
        
    # check whether it defaults to MACE-MP-0
    if file == "macemp0":
        macemp = mace_mp(default_dtype=default_dtype)
    else:
        print("Using "+os.path.basename(file)+" as MACE model.")
        if gpu == True:
            import torch
            print("CUDA found:",torch.cuda.is_available())
            macemp = MACECalculator(model_paths=file,default_dtype=default_dtype,device="cuda")
        else:
            macemp = MACECalculator(model_paths=file,default_dtype=default_dtype,device="cpu")

    # activate dispersion, if requested
    if dispersion == 1:
        from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator as DFTD3
        macemp = DFTD3(dft=macemp)
    elif dispersion == 2:
        from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator as DFTD3
        macemp = DFTD3(dft=macemp,damping="bj")
    elif dispersion == 3:
        from dftd3.ase import DFTD3 as Simple_DFTD3
        print("Using Simple DFTD3.")
        macemp = Calc_Sum([macemp,Simple_DFTD3(method="PBE",damping="d3bj")])
    elif dispersion == 4:
        print("Using internal MACE-MP-0 dispersion.")
        macemp = mace_mp(default_dtype=default_dtype,dispersion=True)
    return macemp


















