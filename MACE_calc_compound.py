#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 13:29:50 2024

@author: bt308570

contains everything involving the compound
"""

import os

class compound:
    def __init__(self,file,structure):
        self.structure = structure
        self.file = file
        
        self.multiplicity = 1
        
        self.get_name()
        
    def get_name(self):
        if "supercell" in self.file:
            multi = self.file.split("_")[2]
            multi = multi.split(".")[0]
            multi = int(multi.split("w")[1])
            self.multiplicity = multi
            
        name = os.path.basename(self.file)
        self.name = name.rsplit(".")[0]

