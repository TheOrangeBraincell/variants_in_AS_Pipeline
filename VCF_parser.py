# -*- coding: utf-8 -*-
"""
date: 10-05-22
title: VCF_parser.py
author: Mirjam Karlsson-Müller

Description:
    Filters vcf files for one or more samples, and summarizes the genotypes
    in a matrix (location x sample).


Instructions:
    Run in command line.
    
    #with coordinates for Estrogen Receptor
    python VCF_parser.py -s . -o vcf_parser_ESR1.txt -c "chr6:151690496-152103274"

    #no coordinates
    python VCF_parser.py -s . -o vcf_parser_all.txt

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
                    help="""Output txt file, containing genotype table, positions
                    with no variant in a sample but variants in other sample(s)
                    are marked with NAN.""")
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
        
#%% 1. Finding all VCF files, Reading them in + Filtering

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
sample_vcf=dict()
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
        
        #print(info)
        "Filter entries based on info string"
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
        
        #If coordinates are given:
        if args.coordinates:
            if chrom != coord_chrom:
                continue
            else:
                if int(position)< int(coord_start) or int(position)> int(coord_stop):
                    continue
        
        "Create dictionary entry, for those who pass filtering"
        coord=chrom+"_"+position
        sample_name=re.search(r"SAMPLE=(S\d+)\.", info).group(1)
        sample_names.append(sample_name)
        sample_vcf[coord]={sample_name:genotype}
        
    current_number_files+=1
    percentage=100*(current_number_files/total_files)
    print("Filtering vcf: {:.2f}%".format(percentage),end="\r")

print("Filtering vcf: Done!            \n",end="\r")

#%% 2. Order the dictionary by positions

location_dict=dict()
chromosomes_numbers=[]
chromosomes_letters=[]
for coord in sample_vcf:
    #Only samples on known chromosomes. 
    if not coord.startswith("chr"):
        continue
    chrom=coord.split("_")[0].strip("chr")
    if chrom in ["X", "Y", "M"]:
        chromosomes_letters.append(chrom)
    else:
        chromosomes_numbers.append(int(chrom))
    if coord.split("_")[0] in location_dict:
        location_dict[coord.split("_")[0]].append(int(coord.split("_")[1]))
    else:
        location_dict[coord.split("_")[0]]=[int(coord.split("_")[1])]
    

#sort chromosomes
chromosomes_letters=sorted(list(set(chromosomes_letters)))
chromosomes_numbers=sorted(list(set(chromosomes_numbers)))
chromosomes_numbers.extend(chromosomes_letters)
chromosomes=chromosomes_numbers

#sort positions per chromosome
for coord in location_dict:
    location_dict[coord]=sorted(location_dict[coord])


#Assemble ordered dict:
ordered_genotype_dict=OrderedDict()
for chromosome in chromosomes:
    for position in location_dict["chr"+str(chromosome)]:
        ordered_genotype_dict["chr"+str(chromosome)+"_"+str(position)]=sample_vcf["chr"+str(chromosome)+"_"+str(position)]

sample_names=sorted(list(set(sample_names)))

print("Sorting Variants: Done!                     \n",end="\r")
    
#%% 3. Write output file

#initiate progress update
total_genotypes=len(ordered_genotype_dict)
current_genotype=0
percentage=100*current_genotype/total_genotypes
print("Writing Output table: {:.2f}%".format(percentage),end="\r")

with open(args.out, "w") as out:
    title="#Location\t"+"\t".join(sample_names)
    out.write(title+"\n")
    for position in ordered_genotype_dict:
        #make sample list, insert spaceholders
        entry=""

        for sample in sample_names:
            variant_found=False
            for key, value in ordered_genotype_dict[position].items():
                if key==sample:
                    entry+="\t"+value
                    variant_found=True
            if variant_found==False:
                entry+="\tNAN"

        
        out.write("{}{}\n".format(position, entry ))
        current_genotype+=1
        percentage=100*current_genotype/total_genotypes
        print("Writing Output table: {:.2f}%".format(percentage),end="\r")

print("Writing Output table: Done!               \n",end="\r")
#%% End time

print("Run time: {:.2f} seconds.           ".format(time.time()-start_time))

























