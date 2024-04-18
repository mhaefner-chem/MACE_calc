#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 13:22:35 2024

@author: bt308570

# input processing for MACE_calc

TODO:
    implement way to read temperatures as TMIN TMAX N
"""

import sys,os
from ase.io import read as ase_read
import MACE_calc_compound as compound

class settings:
    def __init__(self):
        # set all defaults
        self.filename = "MACE_calc.inp"
        self.verbosity = 1
        self.procedure = "analysis"
        
        self.model = "macemp0"
        self.acc = "float64"
        self.dispersion = 0
        
        # optimization
        self.fmax = 0.01
        self.max_steps = 250
        
        # phonons
        self.min_len = 12.0
        self.temperatures = [300]
        
        self.max_time = 72000 # 72000 s/20 h as time limit
        
        
        
    def read_file(self,filename):
        # handle the settings file
        flip = False
        try:
            open(filename,mode="r")
            print("Processing the input file {}.".format(filename))
            flip = True
        except:
            print("No valid settings/input file found. Proceeding with default settings.")
                
        if flip == True:
            with open(filename,mode="r") as f:
                lines = f.readlines()
                for line in lines:
                    # print(line)
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
                    elif "fmax" in line.lower():
                        self.fmax = float(line.split()[-1])
                    elif "max_steps" in line.lower():
                        self.max_steps = int(line.split()[-1])
                    elif "max_time" in line.lower():
                        self.max_time = float(line.split()[-1])
                    elif "proc" in line.lower():
                        self.procedure = line.split()[-1]
                    elif "min_len" in line.lower():
                        self.min_len = float(line.split()[-1])
                    
    
    # constructs the calculator
    def set_calculator(self):
        self.calculator = setup_model(self.model,self.dispersion,self.acc)
        return self.calculator





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
        elif argument in options["model"]:
            settings_calculation.model = sys.argv[i+1]
        elif argument in options["verbosity"]:
            settings_calculation.verbosity = int(sys.argv[i+1])
    
    settings_calculation.read_file(filename)
    
    procedures = {}
    procedures["analysis"] = ["analysis","an","a"]
    procedures["sp"] = ["sp","single","singlepoint","s"]
    procedures["opt"] = ["opt","optimization","o"]
    procedures["phon"] = ["p","phon","phonon"]
    procedures["md"] = ["m","md"]
    
    print(settings_calculation.procedure)
    
    for key in procedures:
        if settings_calculation.procedure in procedures[key]:
            procedure = key
            
    return compounds, settings_calculation, procedure


# prints the help dialogue
import MACE_calc_util as util
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
    return macemp


















