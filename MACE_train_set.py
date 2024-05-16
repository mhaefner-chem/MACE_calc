#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:48:01 2024

@author: bt308570

generates the training set and input for MACE training
"""

def help_output(options):
    util.print_separator("-",64,skip=False)
    print('''
This program processes a list of vasprun.xml files from structural optimizations with VASP and condenses the trajectories into a set of sufficiently distinct structures. It also generates an executable snippet for subsequent MACE training.
''',end="")
    util.print_separator("-",64,skip=False)
    print("Basic usage:")
    print('''MACE_train_set.py -l <XML_LIST> -s <SETTINGS_FILE> <OPTIONS>''')
    util.print_separator("-",64,skip=False)
    print("Options:")
    formatting = "  {:<24} {:<32}"
    print(formatting.format("-h [ --help ]","Prints this help message."))
    print(formatting.format("-l [ --list ] <ARG>","<ARG> is a list of valid vasprun.xml files."))
    print(formatting.format("-s [ --settings ] <ARG>","<ARG> contains all settings for the MACE training."))
    print(formatting.format("-b [ --base ] <ARG>","<ARG> sets the cutoff for structure filtering."))
    print("    Default: 1.0 for ΔE = 1.0 eV, RMSD F = 0.33 eV/Å, RMSD S = 0.01 eV/Å², and RMSD coord = 0.1 Å")
    print()
    print(formatting.format("-d [ --dry]","Does a dry run without creating any files or directories."))
    

    sys.exit()
    

def evaluate_xml(base,xml="vasprun.xml"):
    steps = read(xml, format="vasp-xml",index=":")

    # collector for results of the calculation
    results = []
    for step in steps:
        results.append(step.calc.results)

    # determine important steps
    crucial_steps = []

    # ensure that the first and last step are included
    crucial_steps.append(steps[0])
    crucial_steps.append(steps[-1])

    # filter for energy difference, force difference, stress difference, and positional changes
    determiners = {}
    determiners["e_delta"] = []
    determiners["f_rmsd"] = []
    determiners["s_rmsd"] = []
    determiners["p_rmsd"] = []
    indices = []
    

    selected = 0
    # determine deltas to previously selected step
    for i in range(1,len(steps)):
        # energy difference
        delta_E = results[i]["energy"] - results[selected]["energy"]
        
        # RMSD stress
        parameter = []
        estimator = []
        for j in range(6):
            parameter.append(results[selected]["stress"][j])
            estimator.append(results[i]["stress"][j])
        
        rmsd_S = util.get_rmsd(parameter, estimator)
        
        # RMSD forces
        parameter = []
        estimator = []
        for j in range(len(results[i]["forces"])):
            
            for k in range(3):
                parameter.append(results[selected]["forces"][j][k])
                estimator.append(results[i]["forces"][j][k])
        
        rmsd_F = util.get_rmsd(parameter, estimator)
        

        # RMSD positions
        curr_pos = steps[i].get_positions()
        prev_pos = steps[i-1].get_positions()
        delta = []
        
        curr_list = []
        prev_list = []
        
        for j in range(len(curr_pos)):
            delta.append(util.get_distance(curr_pos[j], prev_pos[j]))
            for k in range(3):
                curr_list.append(curr_pos[j][k])
                prev_list.append(prev_pos[j][k])
        
        rmsd = util.get_rmsd(prev_list, curr_list)
        
        determiners["e_delta"].append(abs(delta_E))
        determiners["f_rmsd"].append(rmsd_F)
        determiners["s_rmsd"].append(rmsd_S)
        determiners["p_rmsd"].append(rmsd)
        indices.append(i)
        
        # cutoff = [0.4,0.4,0.04,0.4]
        cutoff = [base,base/3,base/100,base/10]
        
        if abs(delta_E) > cutoff[0] or rmsd_F > cutoff[1] or rmsd_S > cutoff[2] or rmsd > cutoff[3]: # or f_max > 0.3 or s_max > 0.3:
            selected = i
            
            if not i == len(steps)-1:
                crucial_steps.append(steps[i])


    print("Reduced number of structures from {} to {}.".format(len(steps),len(crucial_steps)))


    # write("for_train.xyz", images = crucial_steps, format="extxyz", append=True)
    return steps,crucial_steps


def find_target_dir(target):
    origin = os.getcwd()
    
    while True:
        if not os.path.isdir(target):
            os.chdir("..")
        else:
            print("Directory {} found.".format(target))
            target_dir = os.path.join(os.getcwd(),target)
            break
        if os.getcwd() == "/":
            print("Directory {} not found. Using current directory {} instead.".format(target,origin))
            target_dir = origin
            break
    os.chdir(origin)
    return target_dir
    

class settings:
    def __init__(self,xml_list,MACE_settings):
        # data set
        self.xml_list = xml_list
        self.MACE_settings = MACE_settings
        
        # MACE settings
        
        self.parameters = {}
        self.parameters["name"] = "MACE"
        self.parameters["mode"] = "scratch"
        
        # from scratch
        self.parameters["model"] = "MACE"
        self.parameters["hidden_irreps"] = "128x0e+128x1o"
        self.parameters["r_max"] = 5.0
        
        # using foundation model
        self.parameters["foundation_model"] = "medium"
        self.parameters["lr"] = 0.01
        self.parameters["scaling"] = "rms_forces_scaling" 
        
        # training set
        self.parameters["train_file"] = "train.xyz"
        self.parameters["valid_fraction"] = 0.05
        
        self.parameters["batch_size"] = 5
        self.parameters["max_num_epochs"] = 128
        
        # swa
        self.parameters["swa"] = True
        self.parameters["start_swa"] = 100
        
        # output and processing
        self.parameters["default_dtype"] = "float64"
        self.parameters["device"] = "cuda"
        self.parameters["seed"] = 3
        self.parameters["save_cpu"] = True
        
        # advanced
        self.parameters["energy_weight"] = 1.0
        self.parameters["forces_weight"] = 1.0
        self.parameters["e0s"] = "average"
        
        self.parameters["ema"] = True
        self.parameters["ema_decay"] = 0.99
        self.parameters["amsgrad"] = True
        
        self.parameters["restart_latest"] = False
        
        self.read_settings()
        
        
    def read_settings(self):

        ints = ["batch_size","max_num_epochs","start_swa","seed"]
        floats = ["r_max","lr","valid_fraction","energy_weight","forces_weight","ema_decay"]
        bools = ["swa","save_cpu","ema","amsgrad","restart_latest"]            

        # handle the settings file            
        with open(self.MACE_settings,mode="r") as f:
            lines = f.readlines()
            for line in lines:
                for key in self.parameters:
                    if "#" in line.lower():
                        continue
                    elif key in line.lower():
                        if key in ints:
                            self.parameters[key] = int(line.split()[-1])
                        elif key in floats:
                            self.parameters[key] = float(line.split()[-1])
                        elif key in bools:
                            self.parameters[key] = bool(line.split()[-1])
                        else:
                            self.parameters[key] = line.split()[-1]
              
        # set name that includes all relevant info
        dir_name = self.parameters["name"]
        dir_name += "_"+self.parameters["mode"]
        dir_name += "_"+str(self.parameters["batch_size"])
        dir_name += "_"+str(self.parameters["start_swa"])
        dir_name += "_"+str(self.parameters["max_num_epochs"])     
        self.parameters["name"] = dir_name          
             
        util.print_separator("=")
        print("Parameters for MACE training")
        for parameter in self.parameters:
            def print_item(text,item):
                format = "  {:16}: {}"
                print(format.format(text,item))
            print_item(parameter,self.parameters[parameter])
        util.print_separator("=")
                            

                        
def write_exec(program, settings):
    with open(settings.parameters["name"]+".inp", mode="w") as f:
        f.write(program)
        valid_parameters = settings.parameters
        if valid_parameters["mode"] == "scratch":
            exclude = ["mode","foundation_model","lr","scaling"]
            for item in exclude:
                valid_parameters.pop(item)
        elif valid_parameters["mode"] == "finetune":
            exclude = ["mode","model","hidden_irreps","r_max"]
            for item in exclude:
                valid_parameters.pop(item)
        else:
            print("Training mode not recognized. Please use 'scratch' or 'finetune'.")
            sys.exit()
        
        for key in valid_parameters:
            f.write(" ")
            if settings.parameters[key] == False:
                continue
            elif isinstance(settings.parameters[key],bool) == True:
                f.write("--"+key)
            elif key == "e0s":
                f.write("--E0s='"+settings.parameters[key]+"'")
            elif isinstance(settings.parameters[key],str) == True:
                f.write("--"+key+"='"+settings.parameters[key]+"'")
            else:
                f.write("--"+key+"="+str(settings.parameters[key]))
                
        f.write("\n")
        

    # exec_line += '--foundation_model=medium --train_file="train.xyz" --valid_fraction=0.05 --energy_weight=1.0 --forces_weight=1.0 --E0s="average" --lr=0.01 --scaling="rms_forces_scaling" --batch_size=5 --max_num_epochs=400 --swa --start_swa=320 --ema --ema_decay=0.99 --amsgrad --default_dtype=float64 --device=cuda --seed=3 --save_cpu\n'
    
    



from ase.io import read, write
import sys, os, math


if __name__ == "__main__":

    # fetch a link to the modules
    if not os.environ["PYTHON_MODULES"] in sys.path:
        sys.path.append(os.environ["PYTHON_MODULES"])
    
    
    import MACE_util as util
    
    # evaluate the command line parameters
    
    options = {}
    options["help"] = ["-h","--help"]
    options["xml_list"] = ["-l","--list"]
    options["settings"] = ["-s","--sett","--settings"]
    options["dryrun"] = ["-d","--dry-run","--dry_run"]
    options["base"] = ["-b","--base"]
    
    base = 1.0
    
    n = len(sys.argv)

    dryrun = False
    for i in range(n):
        argument = sys.argv[i].lower()
        
        if argument in options["help"]:
            help_output(options)
        elif argument in options["xml_list"]:
            xml_list = sys.argv[i+1]
        elif argument in options["settings"]:
            settings_list = sys.argv[i+1]
        elif argument in options["base"]:
            base = float(sys.argv[i+1])
        elif argument in options["dryrun"]:
            dryrun = True
    
    setup_settings = settings(xml_list,settings_list)
    
    
    # evaluate the given list of xml-files
    structures = []
    
    with open(setup_settings.xml_list, mode="r") as f:
        lines = f.readlines()
        n_xmls = len(lines)
        i = 0
        n_tot = 0
        n_cond = 0
        
        for line in lines:
            i += 1
            print("Processing {}. {}/{}".format(line.rstrip("\n"),i,n_xmls))
            xml_structures, xml_structures_condensed = evaluate_xml(base,line.rstrip("\n"))
            n_tot += len(xml_structures)
            n_cond += len(xml_structures_condensed)
            
            print("Current amount of structures: {}/{} ({:5.2f} %).".format(n_cond,n_tot,n_cond/n_tot*100.0))
            print("Predicted amount: {}".format(int(n_xmls/i*n_cond)))
            
            for structure in xml_structures_condensed:
                structures.append(structure)
            
            if not line == lines[-1]:
                util.print_separator("-")
                
    util.print_separator("=")
    # return suggestion for MACE settings
    
    batches = [1,5,10]
    adjusted_batch = math.floor(n_cond/n_xmls)
    batches.append(adjusted_batch)
    
    print("Amount of structures for training: {}".format(n_cond))
    print()
    print("Batch size | max_epochs | swa_start at 80%")
    
    for i in batches:
        max_epoch = math.ceil(200000 * i/n_cond)
        if i == batches[-1]:
            print()
            print("Adjusted to batch size = structures/compound:")
        print("{:3} {:6} {:6}".format(i,max_epoch,math.ceil(max_epoch*0.8)))
        
    util.print_separator("=")
    
    
    
    # create directory for MACE training, add file with structures
    super_dir = find_target_dir("SUPERCELL")
    
    if os.path.basename(super_dir) == "SUPERCELL":
        if not os.path.isdir(os.path.join(super_dir,"PREOPT")):
            os.mkdir(os.path.join(super_dir,"PREOPT"))
        os.chdir(os.path.join(super_dir,"PREOPT"))
    
    if dryrun == False:
        try:
            os.mkdir(setup_settings.parameters["name"])
        except:
            print("Folder already exists so creation was omitted.")
        os.chdir(setup_settings.parameters["name"])
        
        
        write("train.xyz", images = structures, format="extxyz")
        
        
        # create an input snippet for a submit script/direct submit
        program = os.environ["MACE_TRAIN"]
        
        write_exec(program,setup_settings)
        
        