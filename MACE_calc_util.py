#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:41:06 2024

@author: bt308570

utilities for MACE_calc
"""

def print_separator(symbol,n=32,skip=True):
    if skip == True:
        print("")
    for i in range(n):
        print(symbol,end="")
    if skip == True:
        print("")