# -*- coding: utf-8 -*-
"""
Date: Mon Apr  3 11:00:06 2023
File Name: PSI_Counts.py
Author: Mirjam Karlsson-Müller

Description: In an attempt to not overwhelm R when plotting, I made count tables for each PSI table containing counts per bin.
This is to merge the count tables.

    

Useage:
    
    
Possible Bugs:
"""
#%% Imports

import sys

#%% Input


input_files=sys.argv[1:]


#%% Process

#initiate bin counts.
bin_counts={"[0,0.1)":0, "[0.1, 0.2)":0, "[0.2, 0.3)":0, "[0.3, 0.4)":0, "[0.4, 0.5)":0,
            "[0.5, 0.6)":0, "[0.6, 0.7)":0, "[0.8, 0.9)":0, "[0.9, 1]":0}

#Overall dict
Scores={"AA":bin_counts, "AD":bin_counts, "CE":bin_counts, "IR":bin_counts}

for file in input_files:
    with open(file, "r") as psi_counts:
        for line in psi_counts:
            if line.startswith("event"):
                #thats the header
                header=line
                continue
            
            #now we only got entry lines left!
            event= line.split("\t")[0]
            binn=line.split("\t")[1]
            count=line.strip("\n").split("\t")[2]
            Scores[event][binn]+=count

with open("PSI_Tables_Counts.tsv", "w") as out:
    out.write(header)
    for AS in list(Scores.keys()):
        line=AS
        for b in Scores[AS]:
            line+="\t"+str(Scores[AS][b])
        line+="\n"
        out.write(line)
            
            
            
            
            
            
            
            
            
            
            