# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 10:53:01 2022

name: format_interesting_samples.py
author: Mirjam

Description: Just a quick reformatting of how the samples are given in the table of the TripNeg Dataset.

The sample/folder names were extracted from "\\wsl.localhost\Ubuntu\home\mirjam\Trip_Neg_data\PD_ID_External_IDs_WGS_GSE96058"
by cutting out the right column in bash and saving it as "interesting_samples.txt"
"""

with open("interesting_samples.txt", "r") as file, open("formatted_samples_of_interest.txt","w") as out:
    for line in file:
        Sample_name=line.split(".")[0]
        folder_name=(".").join(line.split(".")[1:])
        out.write(Sample_name+"/"+folder_name)
        
        
        
        
