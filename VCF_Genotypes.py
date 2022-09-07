# -*- coding: utf-8 -*-
"""
date: 07-09-22
title: VCF_Genotypes.py
author: Mirjam Karlsson-Müller

Description:
    Filters vcf files for one or more samples, and summarizes the genotypes
    in a matrix (location x sample). For gaps in the table, uses gene expression
    data to predict genotype HMZR or NOEX.
    
Abbreviations:
    HMZA= Genotype Homozygous Alternative Allele
    HETZ= Heterozygous
    HMZR= Homozygous Reference Allele
    NEDA= Not enough data to be sure. Means there is evidence of variants, or
            alternative genotypes, but they got filtered out.
    NOEX= no expression at this point.
    


Instructions:
    Run in command line.
    
    #with coordinates for Estrogen Receptor
    python VCF_Genotypes.py -s . -o genotype_ESR1.tsv -c "chr6:151690496-152103274"

    #no coordinates
    python VCF_Genotypes.py -s . -o genoytpe_wholegenome.tsv

Possible Bugs:
    -
    

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
parser.add_argument('--out', '-o', required=True,
                    help="""Output tsv file, containing genotype table.""")
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")

args = parser.parse_args()

# Extract input coordinates, check their format.
if args.coordinates:
    if re.search(r'[a-z]{3}[MXY]?\d*:\d+-\d+', args.coordinates):
        coord = re.search(r'([a-z]{3}[MXY]?\d*):(\d+)-(\d+)', args.coordinates)
        coord_chrom = coord.group(1)
        coord_start = int(coord.group(2))
        coord_stop = int(coord.group(3))

    else:
        raise argparse.ArgumentTypeError("""The given coordinates are in the 
                                         wrong format. Input as 
                                         chrX:XXXX-XXXX.""")
        quit()
        

#%% 1. Parsing VCFs. 

#Find all the vcf files in the data folder.
argument_glob=args.samples+"/**/*.vcf.gz"
vcf_file_list=glob.glob(argument_glob, recursive=True)

#If no vcf files are found, quit the program.
if len(vcf_file_list)==0:
    print("""There were no vcf files found in the input folder. Please make 
          sure to use the right input folder. The vcf files can be in any 
          subfolder of the input folder.""")
    quit()

#Files are zipped though.
variants=dict()
sample_names=[]
#Progress update
total_files=len(vcf_file_list)
current_number_files=0
percentage=100*(current_number_files/total_files)
print("Filtering vcf {:.2f}%".format(percentage), end="\r")
for file in vcf_file_list:
    vcf = gzip.open(file, "rt")
    #read through entries
    for line in vcf:
        #Skip headers
        if line.startswith("#"):
            continue
        chrom, position, ID, ref, alt, qual, filt, info, form, sample =line.split("\t")
        genotype=sample.split(":")[-7].strip(" ")
        
        #If coordinates are given, no need to save outside of coordinates.
        if args.coordinates:
            if chrom != coord_chrom:
                continue
            else:
                if int(position)< int(coord_start) or int(position)> int(coord_stop):
                    continue
        
        "Before filtering, add all variants to genotype table."
        sample_name=re.search(r"SAMPLE=(S\d+)\.", info).group(1)
        if sample_name not in sample_names:    
            sample_names.append(sample_name)
        variant_ID=chrom+"_"+position
        if sample_name in variants:
            variants[variant_ID][sample_name]=[chrom, position, "NEDA"]
        else:
            variants[variant_ID]={sample_name:[chrom, position, "NEDA"]}
        
        "Now filter the entries."
        #Keep values with MSI<7
        if re.search(r"MSI=(\d+);", info):
            if int(re.search(r"MSI=(\d+);", info).group(1))>=7:
                continue
        else:
            continue
        
        #print("line 104")
        #Keep values with HMPOL<6
        if re.search(r"HMPOL=(\d+);", info):
            if int(re.search(r"HMPOL=(\d+);", info).group(1))>=6:
                continue
        else:
            continue
        
        #print("line 111")
        #Keep entries with GC_cont < 78%
        if re.search(r"GC_CONT=(0\.\d+);", info):
            if float(re.search(r"GC_CONT=(0\.\d+);", info).group(1))>=0.78:
                continue
        else:
            continue
        
        #print("line 118")
        #Keep variant depth >=5
        if re.search(r"VD=(\d+);", info):
            if int(re.search(r"VD=(\d+);", info).group(1))<5:
                continue
        else:
            continue
        
        #print("line 125")
        #Is a flag, so it will only be there if it applies.
        if re.search(r"low_complexity_region", info):
            continue
        
        #print("line 129")
        #6th column, filter bad quality reads.
        if qual=="." or float(qual)<55:
            continue
        #print("line 133")
        
        if re.search(r"ucsc_rep=([a-z]+);", info):
            #if re.search(r"ucsc_rep=([a-z]+);", info).group(1)=="segdup":
            continue
        
        variants[variant_ID][sample_name][2]= genotype
    
    "Progress updates on number of vcf files."
    current_number_files+=1
    percentage=100*(current_number_files/total_files)
    print("Filtering vcf: {:.2f}%".format(percentage),end="\r")
        
        
print("Filtering vcf: Done!            \n",end="\r")

sample_names=sorted(sample_names)

#%% 2. Read in gene.tsv for each sample for gene expression data.

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
                fpkm=line.split("\t")[7]
                
                tsv_info[tsv_sample_name].append([chrom, start, stop, gene_id, fpkm])
    current_file+=1
    percentage=100*current_file/total_files
    print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

print("Reading in gene.tsv: Done!            \n", end="\r")



#%% 3. Fill holes in just created table. 

with open(args.out, "w") as out:
    for variant in variants:
        new_line=[variant]
        #iterating through sorted sample names, to always have same order.
        for sample in sample_names:
            if sample in variant:
                new_line.append(variant[sample][2])
            else:
                #figure out what genotype is. HMZR or NOEX.
                
     
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
