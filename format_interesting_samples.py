# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 10:53:01 2022

@author: mirja
"""

with open("interesting_samples.txt", "r") as file, open("formatted_samples_of_interest.txt","w") as out:
    for line in file:
        Sample_name=line.split(".")[0]
        folder_name=(".").join(line.split(".")[1:])
        out.write(Sample_name+"/"+folder_name)
        
        
        
        