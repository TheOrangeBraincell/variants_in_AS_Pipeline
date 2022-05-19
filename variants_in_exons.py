# -*- coding: utf-8 -*-
"""
date: 13-05-22
title: variants_in_exons.py
author: Mirjam Karlsson-Müller

Description:
    Takes a genotype table from VCF_parser.py and fills the gap using gene
    expression data from each samples corresponding gene.tsv file.


Instructions:
    Run in command line.
    
    python variants_in_exons.py -g genotype_out.txt -p CE_ESR1_all.txt -o genotypes_in_exons.txt


Possible Bugs:
    -
    

"""

#%% Import

import argparse
import time

#%% Time

start_time=time.time()


#%% argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Variants in Exons',
                                 description="""Filters out variants not in 
                                 exons from genotype table.""")

#
parser.add_argument('--out', '-o', required=True,
                    help="""Output txt file, containing only varians in exons,
                    includes a column for exon coordinates.""")
parser.add_argument('--genotype', '-g', required=True,
                    help="""txt file containing genotype per sample/location.""")
parser.add_argument('--PSI_table', '-p', required=True,
                    help="""output from events_gencode. Table containing PSI 
                    scores per sample and exon.""")

args = parser.parse_args()


#%% 1. Read in exon coordinates from PSI table

exon_list=[]
with open(args.PSI_table, "r") as psi:
    for line in psi:
        if line.startswith("Location"): #new version has #Location. change when PSI table is updated.
            continue
        if line.startswith("#"):
            continue
        coordinates=line.split("\t")[0]
        chrom, start, stop=coordinates.split("_")
        exon_list.append(list(coordinates.split("_")))
        


#%% 2. Read genotype file, write new one as we go along.

"""Only include variants in exons, add a second location column with the 
exon the variant is in """



with open(args.genotype, "r") as genotype, open(args.out, "w") as out:
    for line in genotype:
        if line.startswith("#Location"): 
            sample_names=line.split("\t")[2:]
            out.write(line.split("\t")[0]+"\t"+"exon"+"\t"+"\t".join(line.split("\t")[2:]))
            continue
        elif line.startswith("#"):
            out.write(line)
            continue
        
        location=line.split("\t")[0]
        reference=line.split("\t")[1]
        genotypes=line.strip().split("\t")[2:]
        
        #replace HMZ with reference genotype
        for item in genotypes:
            if item=="HMZ":
                genotypes[genotypes.index(item)]=reference
        
        #split location into chromosome and position
        chromosome=location.split("_")[0]
        position=location.split("_")[1]
            
        for exon in exon_list:
            chrom, start, stop=exon
            if chrom==chromosome and int(start)<=int(position)<=int(stop):
                out.write(location+"\t"+ chrom+"_"+start+"_"+stop+"\t"+"\t".join(genotypes)+"\n")
                    
                        


















#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))


