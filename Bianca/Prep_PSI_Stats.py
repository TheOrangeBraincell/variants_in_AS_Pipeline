# -*- coding: utf-8 -*-
"""
Date: Wed Mar 15 13:36:31 2023
File Name: Prep_PSI_Stats.py
Author: Mirjam Karlsson-Müller

Description: Prepares PSI output tables for statistical testing, by removing rows with <2 numerical PSI values.
As a minimum of 2 PSI values for one event would be needed to make a correlation test possible.
    

"""

#%% Imports

import sys

#%% Input handling

input_file= sys.argv[1]
output_file=sys.argv[2]

#%% Filter and output

#The PSI tables are the results of several PSI tables concatenated, so the headers are available several times.
#Only put a header once!
header=False

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            outfile.write(line)
            continue
        if line.startswith("Location"):
            if header==False:
                outfile.write(line)
                header=True
            continue
        
        #Count numerical PSI values.
        PSI_values=line.strip("\n").split("\t")[1:]
        
        c=0
        for psi in PSI_values:
            if type(psi)==int:
                c+=1
        #If there is minimum two numerical PSI values in that row, keep it.
        if c>1:
            outfile.write(line)
