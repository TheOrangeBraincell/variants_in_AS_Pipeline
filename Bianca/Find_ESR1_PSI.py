# -*- coding: utf-8 -*-
"""
Date: Wed Mar 15 13:56:37 2023
File Name: Find_ESR1_PSI.py
Author: Mirjam Karlsson-Müller

Description: Finds PSI entries related to ESR1 in PSI run outputs. 
    
    
"""

#%% Imports

import sys

#%% Handling inputs

in_files= sys.argv[1:] #This will be all outputs, so roughly 80...
out_file="ESR1_PSI.tsv"

#%% Process

#ESR1 coordinates
chrom="chr6"
start= 151656691
stop= 152129619

#Go through files
#only one header
header=False

out=open(out_file, "w")
for file in in_files:
    with open(file, "r") as infile:
        for line in infile:
            if line.startswith("Location"):
                if header==False:
                    out.write(line)
                continue
            if line.startswith("#"):
                continue
            infostring=line.split("\t")[0]
            
            if infostring.startswith("A"):
                #AA or AD
                #Check chromosome
                if infostring.split("_")[2]!= chrom:
                    continue
                
                #if the coordinate is smaller than start or bigger than stop, then its outside the range.
                elif int(infostring.split("_")[4])<start or int(infostring.split("_")[4]) > stop:
                    continue
            else:
                #CE or IR
                #Check chromosome
                if infostring.split("_")[1]!= chrom:
                    continue
                #If the smaller coordinate is bigger than the stop or the bigger coordinate is smaller thna start, then its outside of range..
                if int(infostring.split("_")[4])<start or int(infostring.split("_")[3])> stop:
                    continue
            
            #If it makes it here, print it.
            out.write(line)         

out.close()
                
                
            