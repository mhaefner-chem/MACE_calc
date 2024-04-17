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

class settings:
    def __init__(self):
        # set all defaults
        self.filename = "MACE_calc.inp"
        self.verbosity = 1
        
        self.model = "macemp0"
        self.acc = "float64"
        self.dispersion = 0
        
        self.fmax = 0.01
        self.max_steps = 250
        
        self.max_time = 72000 # 72000 s/20 h as time limit
        
        self.read_file(self.filename)
        
    def read_file(self,filename):
        # handle the settings file
        try:
            with open(filename,mode="r") as f:
                lines = f.readlines()
        except:
            print("No settings/input file found. Proceeding with default settings.")
    
    # constructs the calculator
    def set_calculator(self):
        self.calculator = setup_model(self.model,self.dispersion,self.acc)
        return self.calculator



# reads the arguments and extracts the information
def read_arguments(arguments):
    settings_calculation = settings()
    options = {}
    options["help"] = ["-h","--help"]
    options["compound"] = ["-c","--cell","--compound"]
    options["compound_list"] = ["-l","--list"]
    options["dispersion"] = ["-d","--d3","--dispersion","--d3bj"]
    options["model"] = ["-m","--model"]
    options["procedure"] = ["-p","--proc","--procedure"]
    options["settings"] = ["-s","--setting","--settings"]
    options["verbosity"] = ["-v","--verbosity"]
    
    compounds = []
    procedure_input = "analysis"    
    
    n = len(sys.argv)

    for i in range(n):
        argument = sys.argv[i].lower()
        
        if argument in options["help"]:
            help_output(options)
        elif argument in options["compound"]:
            for compound in read_compound(sys.argv[i+1]):
                compounds.append(compound)
        elif argument in options["compound_list"]:
            compounds.append(read_compound_list(sys.argv[i+1]))
        elif argument in options["dispersion"]:
            settings_calculation.dispersion = 2
        elif argument in options["procedure"]:
            procedure_input = sys.argv[i+1]
        elif argument in options["settings"]:
            settings_calculation.filename = sys.argv[i+1]
        elif argument in options["model"]:
            settings_calculation.model = sys.argv[i+1]
        elif argument in options["verbosity"]:
            settings_calculation.verbosity = int(sys.argv[i+1])
                  
    procedures = {}
    procedures["analysis"] = ["analysis","an","a"]
    procedures["sp"] = ["sp","single","singlepoint","s"]
    procedures["opt"] = ["opt","optimization","o"]
    
    for key in procedures:
        if procedure_input in procedures[key]:
            procedure = key
            
    return compounds, settings_calculation, procedure

# prints the help dialogue
def help_output(options):
    print(123)
    for option in options.items():
        print(option)
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


















