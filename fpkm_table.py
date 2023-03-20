# -*- coding: utf-8 -*-
"""
Date: Mon Mar 20 15:11:39 2023
File Name: fpkm_table.py
Author: Mirjam Karlsson-Müller

Description:
    Reads in all gene.tsv files for all samples and saves them in a gene x Sample table,
    which contains a 1 if the fpkm is bigger equal 10, and a 0 if its below.
    
Useage:
    python ../fpkm_table.py -s ../../Sample_Data -o fpkm_table.tsv 
    
Possible Bugs:
"""
#%% Imports

import argparse
import glob
import time
import re


#%% Time

start_time=time.time()

#%% 0. argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='VCF Parser',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     [-c] "chrX:XXXXXX-XXXXXX"',
                                 description="""Creates a genotype table out of
                                 several samples. Containing location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing among \
                        others, the vcf files.')
parser.add_argument('--out', '-o', required=True,
                    help="""Output tsv file, containing genotype table.""")

args = parser.parse_args()


#%% Parse 

#Find all the tsv files in the data folder.
argument_glob=args.samples+"/**/gene.tsv"
tsv_file_list=glob.glob(argument_glob, recursive=True)

#progress updates
total_files=len(tsv_file_list)
current_file=0
percentage=100*current_file/total_files
print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

gene_info=dict()
tsv_info=dict()
sample_names=[]
for file in tsv_file_list:
    tsv_sample_name=file.split("/")[-5]
    with open(file, "r") as tsv:
        for line in tsv:
            #if line starts with E then its a gene id
            if line.startswith("E"):
                chrom=line.split("\t")[2]
                gene_id=line.split("\t")[1]
                strand=line.split("\t")[3]
                start=line.split("\t")[4]
                stop=line.split("\t")[5]
                fpkm=line.split("\t")[7]
                
                sample_names.append(tsv_sample_name)
                gene_info[gene_id]=[chrom, strand, start, stop]
                tsv_info[gene_id]=dict()
                tsv_info[gene_id][tsv_sample_name]=fpkm
                
    current_file+=1
    percentage=100*current_file/total_files
    print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

print("Reading in gene.tsv: Done!            \n", end="\r")

sample_names=sorted(sample_names)

#%% Write output file


current_gene=0
total_genes=len(list(gene_info.keys()))
percentage=100*current_gene/total_genes
print("Printing fpkm table:: {:.2f}%".format(percentage),end="\r")

with open(args.out, "w") as out:
    out.write("Gene\tChrom\tStrand\tStart\tStop\t"+"\t".join(sample_names)+"\n")
    for gene in gene_info:
        #if the fpkm is not available for this sample and gene, we assume assume fpkm <10.
        newline=gene+"\t"+"\t".join(gene_info[gene])
        for sample in sample_names:
            if sample not in tsv_info[gene]:
                newline+="\t0"
                continue
            if float(tsv_info[gene][sample])>=10:
                newline+="\t1"
            else:
                newline+="\t0"
        #finish line.
        out.write(newline+"\n")

        current_gene+=1
        percentage=100*current_gene/total_genes
        print("Printing fpkm table: {:.2f}%".format(percentage),end="\r")
print("Printing fpkm table: Done!                           ",end="\r")            

#%% Time
print("Run time: {:.2f} seconds.".format(time.time()-start_time))   
