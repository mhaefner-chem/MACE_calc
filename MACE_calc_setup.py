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
        
        
        # calculator parameters
        self.model = "macemp0"
        self.acc = "float64"
        self.dispersion = 0
        
        # optimization
        self.f_max = 0.01
        self.max_steps = 250
        self.sym = True
        
        # bulk modulus calculation
        self.bulk_mod_delta = 0.01
        
        # phonons
        self.min_len = 12.0
        self.t_min = 300
        self.t_max = 300
        self.t_steps = 1
        self.temperatures = np.linspace(self.t_min,self.t_max,self.t_steps)
        
        # MD
        self.md_T = 300
        self.step_size = 1.5
        
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
                elif "f_max" in line.lower():
                    self.f_max = float(line.split()[-1])
                elif "max_steps" in line.lower():
                    self.max_steps = int(line.split()[-1])
                elif "max_time" in line.lower():
                    self.max_time = float(line.split()[-1])
                elif "proc" in line.lower():
                    self.procedure = line.split()[-1]
                    
                # bulk modulus settings
                elif "delta_bulk" in line.lower():
                    self.bulk_mod_delta = float(line.split()[-1])
                    
                elif "min_len" in line.lower():
                    self.min_len = float(line.split()[-1])
                elif "md_t" in line.lower():
                    self.md_t = float(line.split()[-1])
                elif "step_size" in line.lower():
                    self.step_size = float(line.split()[-1])
                elif "t_min" in line.lower():
                    self.t_min = float(line.split()[-1])
                elif "t_max" in line.lower():
                    self.t_max = float(line.split()[-1])
                elif "t_steps" in line.lower():
                    self.t_steps = int(line.split()[-1])
                
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
            
        try:
            self.temperatures = np.linspace(self.t_min,self.t_max,self.t_steps)
        except:
            self.temperatures = [300]
    
    # constructs the calculator
    def set_calculator(self):
        self.calculator = setup_model(self.model,self.dispersion,self.acc)
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
            print("")
            
            # calculator parameters
            print("Calculator")
            print("")
            print_item("Model:",self.model)
            print_item("Accuracy:",self.acc)
            print_item("Dispersion (2=D3BJ):",self.dispersion)
            print("")
            
            # optimization
            if self.procedure in ["opt","bulk","phon","analysis"]:
                print("Optimization")
                print("")
                print_item("F_max [eV/Å]:",self.f_max)
                print_item("Max_steps:",self.max_steps)
                print_item("Sym:",self.sym)
                print("")
                
            # bulk modulus calculation
            if self.procedure in ["bulk","analysis"]:
                print("Bulk Modulus")
                print("")
                print_item("Delta [%]:",self.bulk_mod_delta)
                print("")
            
            # phonons
            
            if self.procedure in ["phon","analysis"]:
                print("Phonons")
                print("")
                print_item("Min_len [Å]:",self.min_len)
                print_item("T_min [K]:",self.t_min)
                print_item("T_max [K]:",self.t_max)
                print_item("T_steps:",self.t_steps)
                print("")
            
            # MD
            if self.procedure in ["md","analysis"]:
                print("MD")
                print("")
                print_item("T_MD [K]:",self.md_T)
                print_item("Step_size [fs]:",self.step_size)
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
            settings_calculation.dispersion = 2
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
    

    sys.exit()

# reads a structure from a file
def read_compound(file):
    file_ending = file.split(".")[-1].lower()
    
    if file_ending == "vasp" or file in ["POSCAR","CONTCAR"]:
        atoms = ase_read(file,index=":",format="vasp")
    elif file_ending == "cif":
        atoms = ase_read(file,index=":",format='cif')
    elif file_ending == "xyz":
        atoms = ase_read(file,index=":",format='extxyz')
    else:
        atoms = ase_read(file,index=":")
        
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
from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator as DFTD3
from dftd3.ase import DFTD3 as Simple_DFTD3
from ase.calculators.mixing import SumCalculator as Calc_Sum

def setup_model(file,dispersion,default_dtype="float64"):
        
    # check whether it defaults to MACE-MP-0
    if file == "macemp0":
        macemp = mace_mp(default_dtype=default_dtype)
    else:
        print("Using "+os.path.basename(file)+" as MACE model.")
        macemp = MACECalculator(model_paths=file,default_dtype=default_dtype,device="cpu")

    # activate dispersion, if requested
    if dispersion == 1:
        macemp = DFTD3(dft=macemp)
    elif dispersion == 2:
        macemp = DFTD3(dft=macemp,damping="bj")
    elif dispersion == 3:
        print("Using Simple DFTD3.")
        macemp = Calc_Sum([macemp,Simple_DFTD3(method="PBE",damping="d3bj")])
    elif dispersion == 4:
        print("Using internal MACE-MP-0 dispersion.")
        macemp = mace_mp(default_dtype=default_dtype,dispersion=True)
    return macemp


















