# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:40:14 2022

python variants_vs_FPKM_one_sample.py -s ../Sample_Data/ -c "chr6:151650000-152130000" 

@author: mirja
"""

"""The goal is to get an idea how high the gene expression/how many reads there
need to be for us to be able to make a reliable genotype prediction. 
(Though we will eventually probably consult the DNA data)"""

#%% Imports

import argparse
import glob
import re
import gzip
import matplotlib.pylab as plt
import numpy as np
import math

#%% Argparse
"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Genotype Plots',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     [-c] "chrX:XXXXXX-XXXXXX"',
                                 description="""Creates a variant table out of
                                 several samples. Containing location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing among \
                        others, the vcf files.')
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


#%% Idea 1: Var/kb over FPKM

"For 1 sample, but all genes"

#need fpkms and variant files.
"Read in fpkms"

#Find all the tsv files in the data folder.
argument_glob=args.samples+"/**/gene.tsv"
tsv_file_list=glob.glob(argument_glob, recursive=True)

#progress updates
total_files=len(tsv_file_list)
current_file=0
percentage=100*current_file/total_files
print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

fpkm=dict()
for file in tsv_file_list:
    with open(file, "r") as tsv:
        for line in tsv:
            if line.startswith("Gene"):
                continue
            gene_name=line.split("\t")[1]
            fpkm[gene_name]=line.split("\t")[4:7]
            fpkm[gene_name].append(int(fpkm[gene_name][1])-int(fpkm[gene_name][0]))
    current_file+=1
    percentage=100*current_file/total_files
    print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

print("Reading in gene.tsv: Done!            \n", end="\r")


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
#Progress update
total_files=len(vcf_file_list)
current_number_files=0
percentage=100*(current_number_files/total_files)
print("Filtering vcf {:.2f}%".format(percentage), end="\r")
variant_coord=[]
for file in vcf_file_list:
    vcf = gzip.open(file, "rt")
    #read through entries
    for line in vcf:
        #Skip headers
        if line.startswith("#"):
            continue
        chrom, position, ID, ref, alt, qual, filt, info, form, sample =line.split("\t")
        genotype=sample.split(":")[-7].strip(" ")
        
        #If coordinates are given:
        if args.coordinates:
            if chrom != coord_chrom:
                continue
            else:
                if int(position)< int(coord_start) or int(position)> int(coord_stop):
                    continue
    
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
        

        "Add coordinates, for those who pass filtering"
        variant_coord.append(position)
        
        
        
    current_number_files+=1
    percentage=100*(current_number_files/total_files)
    print("Filtering vcf: {:.2f}%".format(percentage),end="\r")

print("Filtering vcf: Done!            \n",end="\r")


print("Counting...", end="\r")
counts=dict()
for gene in fpkm:
    counts[gene]=0
    start=fpkm[gene][0]
    stop=fpkm[gene][1]
    for position in variant_coord:
        if start<=position<=stop:
            counts[gene]+=1

print("Counting: Done!       \n", end="\r")
#make plot
x=[]
y=[]

for gene in fpkm:
    if float(fpkm[gene][2])!=0:    
        x.append(math.log2(float(fpkm[gene][2])+1))
        y.append(float(counts[gene]*1000/fpkm[gene][3]))
"""
print(x)
print("\n\n")
print(y)
"""
plt.plot(x,y,"ro")
plt.xlabel("log2(fpkm)")
plt.ylabel("Variant Counts/kb")
plt.title("Variant Counts/FPKM for "+args.samples.split("/")[-2])
plt.ylim(0, 2000)
plt.savefig("Variants_vs_FPKM_Sample1_2000.png")



#%% Idea 2: Distributions over expression levels.





#%% Idea 3: allele fractions/total read depth.






