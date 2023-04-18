# -*- coding: utf-8 -*-
"""
Date: Wed Mar 15 08:51:08 2023
File Name: Prep_Genotype_Stats.py
Author: Mirjam Karlsson-Müller


Description: Filters genotype entries for statistical testing. Only keeps locations that have minimum 2
distinct genotypes out of 0/0, 0/1, 1/1. 
"""

#%% Imports

import sys
import os

#%% Input handling

input_file=sys.argv[1]
output_file=sys.argv[2]

#Server option: to not overwrite output file.
if os.path.isfile(output_file):
    print("The output file already exists, we assume it is complete.")
    quit()

#%% Write output

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith("Location"):
            outfile.write(line)
            continue
        elif line.startswith("#"):
            outfile.write(line)
        #Flags to check for genotypes
        counts=[line.count("0/0"), line.count("0/1"), line.count("1/1"), line.count("0/2"), line.count("1/2"), line.count("2/2"), line.count("0/3"), line.count("1/3"),
                line.count("2/3"), line.count("3/3"),]
        #print(counts)
        #j needs to be minimum 2 for it to have min 2 distinct genotypes
        j=0
        for i in counts:
            if i !=0:
                j+=1
        
        if j>=2:
            outfile.write(line)
                
        
        
        
        


