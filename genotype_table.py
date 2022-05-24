# -*- coding: utf-8 -*-
"""
date: 12-05-22
title: genotype_table.py
author: Mirjam Karlsson-Müller

Description:
    Takes a genotype table from VCF_parser.py and fills the NAN gaps, depending
    on whether a variant is within a gene or not. 
    This is enough crude filtering, as later on, we test the relation between
    variants in exons and their PSI score, which will only be a number if there
    has been more than 10 reads found of that exon. Hence we will only look
    at variants who are in exons which are expressed.


Instructions:
    Run in command line.

    python genotype_table.py -s . -o genotype_out.txt -t vcf_parser_ESR1.txt

Possible Bugs:
    -
    

"""




#%% Import

import glob
import argparse
import time

#%% Time

start_time=time.time()


#%% argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Genotype Table',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     -t VCFPARSER TABLE',
                                 description="""Fills gaps in VCF_parser.py's
                                 output table.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing among \
                        others, their gene.tsv files.')
parser.add_argument('--out', '-o', required=True,
                    help="""Output txt file, containing updated
                    genotype table, no more NAN.""")
parser.add_argument('--table', '-t', required=True,
                    help="""output from vcf parser. Table containing variants,
                    empty spots held by NAN.""")

args = parser.parse_args()



#%% Functions

def is_in_gene(sample, location):
    is_expressed=False
    chrom=location.split("_")[0]
    position=location.split("_")[1]
    for gene in tsv_info[sample]:
        start=gene[1]
        stop=gene[2]
        c=gene[0]
        if start<position<stop and chrom==c:
            is_expressed=True
    
    if is_expressed==False:
        return "NE"
    if is_expressed==True:
        return "HMZR"
    
#%% 1. Read in TSV gene file for each sample.    

#Find all the tsv files in the data folder.
argument_glob=args.samples+"/**/gene.tsv"
tsv_file_list=glob.glob(argument_glob, recursive=True)

#progress updates
total_files=len(tsv_file_list)
current_file=0
percentage=100*current_file/total_files
print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

tsv_info=dict()
for file in tsv_file_list:
    tsv_sample_name=file.split("/")[1]
    tsv_info[tsv_sample_name]=[]
    with open(file, "r") as tsv:
        for line in tsv:
            #if line starts with E then its a gene id
            if line.startswith("E"):
                chrom=line.split("\t")[2]
                gene_id=line.split("\t")[1]
                start=line.split("\t")[4]
                stop=line.split("\t")[5]
                
                tsv_info[tsv_sample_name].append([chrom, start, stop, gene_id])
    current_file+=1
    percentage=100*current_file/total_files
    print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

print("Reading in gene.tsv: Done!            \n", end="\r")

#%% 3. Read in VCF_Parser's output table, update NAN values

#progress updates
total_lines=sum(1 for line in open(args.table))
current_line=0
percentage=100*current_line/total_lines
print("Replacing NAN values: {:.2f}%".format(percentage),end="\r")

#read in current genotype table, write new one simultaneously.
with open(args.table, "r") as table, open(args.out, "w") as out:
    out.write("# NE=Not Expressed\n# HMZR=Homozygote Reference, HETZ=Heterozygote, HMZA=Homozygote Alternative\n")
    for line in table:
        if line.startswith("#"):
            sample_names=line.strip().split("\t")[1:]
            out.write(line)
            continue
        location=line.split("\t")[0]
        genotypes=line.strip().split("\t")[1:]
        
        for value in genotypes:
            if value=="NAN":
                new_value=is_in_gene(sample_names[genotypes.index(value)],location)
            elif value=="0/1":
                new_value="HETZ"
            elif value=="1/1":
                new_value="HMZA"
            genotypes[genotypes.index(value)]=new_value
            
        out.write(location+"\t"+"\t".join(genotypes)+"\n")
        current_line+=1
        percentage=100*current_line/total_lines
        print("Replacing NAN values: {:.2f}%".format(percentage),end="\r")

print("Replacing NAN values: Done!            \n", end="\r")
#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))
