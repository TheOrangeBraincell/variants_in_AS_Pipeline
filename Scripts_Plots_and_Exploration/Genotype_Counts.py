# -*- coding: utf-8 -*-
"""
Date: Wed Mar 15 08:51:08 2023
File Name: Prep_Genotype_Stats.py
Author: Mirjam Karlsson-Müller


Description: Filters genotype entries for statistical testing. Only keeps locations that have minimum 2
distinct genotypes out of 0/0, 0/1, 1/1 (incl 0/2, 2/2, 1/2 etc)

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
    #Header for counts
    outfile.write("Location\t0/0\t0/X\tX/X\n")
    for line in infile:
        if line.startswith("Location"):
            continue
        elif line.startswith("#"):
            continue
        #Flags to check for genotypes
        counts=[line.count("0/0"), line.count("0/1"), line.count("1/1"), line.count("0/2"), line.count("1/2"), line.count("2/2"), line.count("0/3"), line.count("1/3"),
                line.count("2/3"), line.count("3/3"),]
        
        #We only want counts for variants with at least 1% alternative allele.
        total=sum(counts)
        hmzr=counts[0]
        hetz=counts[1]+counts[3]+counts[4]+counts[6]+counts[7]+counts[8]
        hmza=counts[2]+counts[5]+counts[9]
        if (total-hmzr)/total > 0.01:
            outfile.write("{}\t{}\t{}\t{}\n".format(line.split("\t")[0], hmzr, hetz, hmza))
                
        
        
        
        


