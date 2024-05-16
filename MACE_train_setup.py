#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 12:03:32 2024

@author: bt308570

This script sets up supercell configurations for optimization with VASP to get data for MACE training
"""

class settings:
    def __init__(self):
        # I/O
        self.filename = "MACE_train_setup.inp"
        # self.procedure = "analysis"
        # self.file_results = "out.out"
        
        self.verbosity = 1        
        
        # calculator parameters
        self.model = "macemp0"
        self.acc = "float64"
        self.dispersion = 0
        
        # optimization
        self.f_max = 0.01
        self.max_steps = 250
        
        # phonons
        self.min_len = 12.0
        self.t_min = 300
        self.t_max = 300
        self.t_steps = 1
        
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
                elif "f_max" in line.lower():
                    self.f_max = float(line.split()[-1])
                elif "max_steps" in line.lower():
                    self.max_steps = int(line.split()[-1])
                elif "max_time" in line.lower():
                    self.max_time = float(line.split()[-1])
                elif "proc" in line.lower():
                    self.procedure = line.split()[-1]
                elif "min_len" in line.lower():
                    self.min_len = float(line.split()[-1])
                elif "md_t" in line.lower():
                    self.md_t = float(line.split()[-1])
                elif "step_size" in line.lower():
                    self.step_size = float(line.split()[-1])
                elif "t_min" in line.lower():
                    t_min = float(line.split()[-1])
                elif "t_max" in line.lower():
                    t_max = float(line.split()[-1])
                elif "t_steps" in line.lower():
                    t_steps = float(line.split()[-1])
                    
                # writes and prints
                elif "write_geom" in line.lower():
                    self.write_geometries = bool(line.split()[-1])
                elif "write_e" in line.lower():
                    self.write_energies = bool(line.split()[-1])
                elif "print_g" in line.lower():
                    self.print_geometries = bool(line.split()[-1])
                elif "print_f" in line.lower():
                    self.print_forces = bool(line.split()[-1])
                elif "print_s" in line.lower():
                    self.print_stresses = bool(line.split()[-1])
        
        
# reads the arguments and extracts the information
def read_arguments(arguments):
    settings_setup = settings()    
    
    options = {}
    options["help"] = ["-h","--help"]
    options["compound"] = ["-c","--cell","--compound"]
    options["dispersion"] = ["-d","--d3","--dispersion","--d3bj"]
    options["model"] = ["-m","--model"]
    options["procedure"] = ["-p","--proc","--procedure"]
    options["settings"] = ["-s","--setting","--settings","-i","--input"]
    options["verbosity"] = ["-v","--verbosity"]
    
    n = len(sys.argv)
    
    compounds = []

    for i in range(n):
        argument = sys.argv[i].lower()
        
        for i in range(0,n):
            # specify the cells and supercell size from cif
            if sys.argv[i] == "-c":
                print("Reading in the file "+sys.argv[i+1]+".")
                compound = Structure.from_file(sys.argv[i+1])
                cif_file = sys.argv[i+1]
                super = 1
            elif sys.argv[i] == "-r" or sys.argv[i] == "--reference":
                print("Reading in the file "+sys.argv[i+1]+".")
                reference = Structure.from_file(sys.argv[i+1])
                ref_file = sys.argv[i+1]
            elif sys.argv[i] == "-d" or sys.argv[i] == "--dopant":
                print("Reading in the file "+sys.argv[i+1]+".")
                dopant = Structure.from_file(sys.argv[i+1])
                dope_file = sys.argv[i+1]
            elif sys.argv[i] == "-s" or sys.argv[i] == "--supercell":
                sc_size = [int(x) for x in sys.argv[i+1].split("x")]
                sc_str = sys.argv[i+1]
            elif sys.argv[i] == "-n" or sys.argv[i] == "--name":
                id_tag = sys.argv[i+1]
            elif sys.argv[i] == "-m" or sys.argv[i] == "--mode":
                mode = sys.argv[i+1]
            elif sys.argv[i] == "-q" or sys.argv[i] == "--coulomb":
                selec = sys.argv[i+1]
            elif sys.argv[i] == "-p" or sys.argv[i] == "--potential":
                mlp_file = sys.argv[i+1]
                
                
        
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
    
    try:
        settings_calculation.read_file(filename)
        print("Processing the input file {}.".format(filename))
    except:
        print("No valid settings/input file found. Proceeding with default settings.")
        

    return settings_setup


import sys, os, time, datetime
# beginning of main program
if __name__ == "__main__":
    
    
    # time_tracker = time_tracker()
    
    # fetch a link to the modules
    if not os.environ["PYTHON_MODULES"] in sys.path:
        sys.path.append(os.environ["PYTHON_MODULES"])
    # print(sys.path)
    
    sys.stdout.flush()
    
    file_dir = os.getcwd()

    settings = read_arguments(sys.argv)
    
    
    
    selec = "100"
    super = 0
    mode = "vasp"
    mlp_file = "SHOULD_NOT_EXIST"
    n = len(sys.argv)
    
        
    
    warnings.filterwarnings("ignore", message="No POTCAR file with matching TITEL") # suppress the POTCAR warning, works regardless
    warnings.filterwarnings("ignore", message="POTCAR data with symbol") # suppress the POTCAR warning, works regardless
    
    
    
    
    print("Running script to carry out supercell calculations.")
    # locate head directory
    workdir = find_topdir()
    os.chdir(workdir)
    
    # import the directory tools 
    spec = importlib.util.spec_from_file_location("directory_tools",workdir + "/bin/directory_tools.py")
    directory_tools = importlib.util.module_from_spec(spec)
    sys.modules["name"] = directory_tools
    spec.loader.exec_module(directory_tools)   
    
    directory_tools.make_directory("SUPERCELL")
    os.chdir("SUPERCELL")
    directory_tools.make_directory("PREOPT")
    directory_tools.make_directory("OPT")
    os.chdir("PREOPT")
    
    
    # set up the calculation for the supercell-generated structures
    if super == 1:
        name_prefix = directory_tools.directory_name_constructor(structure)
        name_suffix = directory_tools.directory_name_constructor(reference)
        
        name = name_prefix+"_"+sc_str+"_"+id_tag+"_base_"+name_suffix
        directory_tools.make_directory(name)
        os.chdir(name)
        
        workdir = os.getcwd()
        
        shutil.copyfile(os.path.join(file_dir,cif_file),cif_file)
        if os.path.isfile(mlp_file):
            shutil.copyfile(os.path.join(file_dir,mlp_file),mlp_file)
        
        # run supercell to generate all relevant structures
        if not os.path.isfile("done"):
            process = subprocess.Popen(['supercell',"-i",cif_file,"-s",sc_str,"-n","l"+selec,"-q","-m"],
                         stdin=subprocess.DEVNULL,
                         start_new_session=False)
            process.wait()
            shutil.copyfile(cif_file,"done")
    
        
        files = glob("supercell*.cif")
        
        # shutil.copyfile(os.path.join(file_dir,ref_file),"ref.cif")
        # shutil.copyfile(os.path.join(file_dir,dope_file),"dope.cif")
        
        counter = 0
        for i in files:
            base_name = i.split(".")[0]
            if mode == "szp":
                base_name = base_name + "_gpaw_szp"
            elif mode == "dzp":
                base_name = base_name + "_gpaw_dzp"
            elif mode == "fast":
                  base_name = base_name + "_fast"
            elif mode == "slow":
                  base_name = base_name + "_slow"
            elif mode == "ml":
                  base_name = base_name + "_ml"
            elif mode == "cifs":
                sys.exit()
                  
                   
            directory_tools.make_directory(base_name)
            shutil.copyfile(i,os.path.join(base_name,i))
            if mode == "ml":
                shutil.copyfile(mlp_file,os.path.join(base_name,mlp_file))
            os.chdir(base_name)
            
            if not os.path.isfile("INCAR"):
                print("Setting up input!")
                counter += 1
                print(counter)
                if mode == "fast":
                    process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",i,"-n",name_prefix+"_"+base_name,"--fast","--ionly"],
                                 stdin=subprocess.DEVNULL,
                                 stdout=open(name_prefix+"_"+base_name+'.out', 'w'),
                                 stderr=subprocess.STDOUT,
                                 start_new_session=True,
                                 preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
                    process.wait()
                elif mode == "slow":
                    process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",i,"-n",name_prefix+"_"+base_name,"--ionly"],
                                 stdin=subprocess.DEVNULL,
                                 stdout=open(name_prefix+"_"+base_name+'.out', 'w'),
                                 stderr=subprocess.STDOUT,
                                 start_new_session=True,
                                 preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
                    process.wait()
                elif mode == "ml":
                    process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",i,"-n",name_prefix+"_"+base_name,"--ionly","--ml"],
                                 stdin=subprocess.DEVNULL,
                                 stdout=open(name_prefix+"_"+base_name+'.out', 'w'),
                                 stderr=subprocess.STDOUT,
                                 start_new_session=True,
                                 preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
                    process.wait()
                    shutil.copyfile("INCAR","INCAR_x")
                    with open("INCAR_x", "rt") as fin:
                        with open("INCAR", "wt") as fout:
                            for line in fin:
                                fout.write(line.replace('Run', 'run'))
                else:
                    
                    process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",i,"-n",name_prefix+"_"+base_name,"-p","--ionly"],
                                 stdin=subprocess.DEVNULL,
                                 stdout=open(name_prefix+"_"+base_name+'.out', 'w'),
                                 stderr=subprocess.STDOUT,
                                 start_new_session=True,
                                 preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
                    process.wait()
            else:
                print("Already set up.")
            
            os.chdir(workdir)
        
    # set up the calculations for the reference structure
    elif super == 0:
        os.chdir(os.path.join(find_topdir(),"SUPERCELL","PREOPT"))
        name = directory_tools.directory_name_constructor(reference)
        directory_tools.make_directory(name)
        os.chdir(name)
        shutil.copyfile(os.path.join(file_dir,ref_file),ref_file)
        process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",ref_file,"-n",name,"-p"],
                     stdin=subprocess.DEVNULL,
                     stdout=open(name+'.out', 'w'),
                     stderr=subprocess.STDOUT,
                     start_new_session=True,
                     preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
    
    # set up the calculations for the dopant structure
    
    # os.chdir(os.path.join(find_topdir(),"SUPERCELL","PREOPT"))
    # name = directory_tools.directory_name_constructor(dopant)
    # directory_tools.make_directory(name)
    # os.chdir(name)
    # shutil.copyfile(os.path.join(file_dir,dope_file),dope_file)
    # process = subprocess.Popen(["python",find_topdir()+"/bin/optimizer.py","-c",dope_file,"-n",name,"-p"],
    #              stdin=subprocess.DEVNULL,
    #              stdout=open(name+'.out', 'w'),
    #              stderr=subprocess.STDOUT,
    #              start_new_session=True,
    #              preexec_fn=(lambda: signal.signal(signal.SIGHUP, signal.SIG_IGN)))
    
    
    sys.exit()

