# -*- coding: utf-8 -*-
"""
date: 27-09-22
title: Variants_AS_altered.py
author: Mirjam Karlsson-Müller

Description:
    Creates a table containing sample, chrom, pos and genotype for each variant entry.
    Returns a filtered and an unfiltered version of the table.
    
    Made to be compared to tnbc variants. To be run on ca. 250 samples that are
    shared between cohorts.
    
Usage:
    python Variants_AS_altered.py -s ../Sample_Data/
"""

#%% Imports

import argparse
import glob
import re
import gzip
import time

#%% Time

start_time=time.time()

#%% 0. argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Variants_AS_altered',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     [-c] "chrX:XXXXXX-XXXXXX"',
                                 description="""Creates a genotype table out of
                                 several samples. Containing location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='file containing sample folder names containing among \
                        others, the vcf, bam and gene.tsv files.')



args = parser.parse_args()


#%% 2. Read in variants and filter.
#Find all the vcf files in the data folder.
argument_glob=args.samples+"/**/*.vcf.gz"
vcf_file_list=glob.glob(argument_glob, recursive=True)

#If no vcf files are found, quit the program.
if len(vcf_file_list)==0:
    print("""There were no vcf files found in the input folder. Please make 
          sure to use the right input folder. The vcf files can be in any 
          subfolder of the input folder.""")
    quit()

#Progress update
total_files=len(vcf_file_list)
current_number_files=0
percentage=100*(current_number_files/total_files)
print("Filtering vcf {:.2f}%".format(percentage), end="\r")
with open("rna_variants_filtered.txt", "w") as filtered, open("rna_variants.txt", "w") as unfiltered:
    filtered.write("#Sample\tchrom\tpos\tgenotype\n")
    unfiltered.write("#Sample\tchrom\tpos\tgenotype\n")
    for file in vcf_file_list:
        vcf = gzip.open(file, "rt")
        sample_name=file.split("/")[-5]
        #read through entries
        for line in vcf:
            #Skip headers
            if line.startswith("#"):
                continue
            chrom, position, ID, ref, alt, qual, filt, info, form, sample =line.split("\t")
            genotype=sample.split(":")[-7].strip(" ")
            
            unfiltered.write("{}\t{}\t{}\t{}\n".format(sample_name, chrom, position, genotype))
            
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
            
            filtered.write("{}\t{}\t{}\t{}\n".format(sample_name, chrom, position, genotype))
        
        "Progress updates on number of vcf files."
        current_number_files+=1
        percentage=100*(current_number_files/total_files)
        print("Filtering vcf: {:.2f}%".format(percentage),end="\r")
            
print("Filtering vcf: Done!            \n",end="\r")


#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))      

































