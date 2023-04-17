# -*- coding: utf-8 -*-
"""
Date: Tue Apr 11 11:26:55 2023
File Name: fpkm_table.py
Author: Mirjam Karlsson-Müller

Description:
    Creates a table containing all fpkm scores in the 3455 samples per gene name.
    
    
List of Functions:
    none

Useage:
    python fpkm_table.py -s ../../Sample_Data/ -o fpkm_table.tsv    
    
Possible Bugs:

"""


#%% Imports

import argparse
import time
import glob
import re

#%% Initialize time


start_time=time.time()

#%% argparse

parser = argparse.ArgumentParser(prog='creating fpkm table',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT',
                                 description="""Creates a fpkm table out of
                                 gene expression information files of samples.""")

parser.add_argument('--samples', '-s', required=True,
                    help='file containing the path to the vcf files.')
parser.add_argument('--out', '-o', required=True,
                    help="""Output file containing fpkm values""")

                    

args = parser.parse_args()


#%% Parsing fpkm tables: gene.tsv files.

#no checking for multiple files per sample, as cohort is already prepped like that.
#Make fpkm file list.

argument_glob=args.samples+"/**/gene.tsv"
tsv_file_list=glob.glob(argument_glob, recursive=True)

#If no vcf files are found, quit the program.
if len(tsv_file_list)==0:
    print("""There were no gene.tsv files found in the input folder. Please make 
          sure to use the right input folder. The files can be in any 
          subfolder of the input folder.""")
    quit()


#Make sample names list.
sample_names=[]
for file in tsv_file_list:
    info=re.search(r'(S\d{6})', file)
    name=info.group(1)
    sample_names.append(name)

#Initialize dictionary
fpkm=dict()

for file in tsv_file_list:
    info=re.search(r'(S\d{6})', file)
    sample=info.group(1)
    with open(file, "r") as infile:
        #go through lines.
        for line in infile:
            #Skip header
            if line.startswith("Gene"):
                continue
            #Otherwise its entries! fpkm[gene][sample]=score
            if line.split("\t")[1] not in fpkm:
                fpkm[line.split("\t")[1]]=dict()
            fpkm[line.split("\t")[1]][sample]=line.split("\t")[7]

#%% Write into output file


with open(args.out, "w") as out:
    #write header
    out.write("Location\t"+"\t".join(sample_names)+"\n")
    for gene in fpkm:
        new_line=gene
        for sample in sample_names:
            #Check if sample has an fpkm
            if sample not in fpkm[gene]:
                new_line+="\t0"
            else:
                new_line+="\t"+fpkm[gene][sample]
        out.write(new_line+"\n")
            

#%% Stop timer


print("Run time: {:.2f} seconds.".format(time.time()-start_time))   

