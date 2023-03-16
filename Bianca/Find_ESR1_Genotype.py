# -*- coding: utf-8 -*-
"""
Date: Wed Mar 15 15:01:51 2023
File Name: Find_ESR1_Genotype.py
Author: Mirjam Karlsson-Müller

Description: Find Variant Locations in Genotype tables that are within the coordinates for ESR1
    
"""


#%% imports

import sys

#%% Input handling

infile=sys.argv[1]

#%% Code

#ESR1 coordinates
chrom="chr6"
start= 151656691
stop= 152129619

header=False
with open(infile, "r") as geno, open("Genotypes_ESR1.tsv", "w") as out:
    for line in geno:
        if line.startswith("Location"):
            if header==False:
                out.write(line)
                header=True
            continue
        if line.startswith("#"):
            continue
        
        #First column contains location identifier string. "chr_coord_refbase_(altbase)
        info=line.split("\t")[0]
        
        #Are we on chromosome 6?
        if info.split("_")[0]!= chrom:
            continue
        
        #Are coordinates in gene range of ESR1?
        if int(info.split("_")[1]) <start or int(info.split("_")[1])> stop:
            continue
        
        #If it passes write into output.
        out.write(line)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

