# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 09:52:19 2022

@author: mirja
"""


#%% Imports

import argparse
import glob
import re
import gzip
from collections import OrderedDict
import time

#%% Time

start_time=time.time()

#%% 0. argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='VCF Parser',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     [-c] "chrX:XXXXXX-XXXXXX"',
                                 description="""Creates a variant table out of
                                 several samples. Containing location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing among \
                        others, the vcf files.')


args = parser.parse_args()

        
#%% Do we have 1 vcf per bam?

argument_glob=args.samples+"/**/*.vcf.gz"
vcf_file_list=glob.glob(argument_glob, recursive=True)

argument_glob=args.samples+"/**/*.bam"
bam_file_list=glob.glob(argument_glob, recursive=True)

argument_glob=args.samples+"/**/gene.tsv"
tsv_file_list=glob.glob(argument_glob, recursive=True)


print(len(vcf_file_list), len(bam_file_list), len(tsv_file_list))
quit()

#%% How big is the percentage of genes having fpkm <1 per sample?

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
    tsv_sample_name=file.split("/")[-5]
    tsv_info[tsv_sample_name]=[]
    with open(file, "r") as tsv:
        for line in tsv:
            #if line starts with E then its a gene id
            if line.startswith("E"):
                chrom=line.split("\t")[2]
                gene_id=line.split("\t")[1]
                start=line.split("\t")[4]
                stop=line.split("\t")[5]
                fpkm=line.split("\t")[7]
                
                tsv_info[tsv_sample_name].append([chrom, start, stop, gene_id, fpkm])
    current_file+=1
    percentage=100*current_file/total_files
    print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

print("Reading in gene.tsv: Done!            \n", end="\r")