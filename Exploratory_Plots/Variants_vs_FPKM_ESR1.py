# -*- coding: utf-8 -*-
"""
Exploring Genotype Prediction/Gene Expression Coherence.

#Variants vs FPKM  for one gene!

python Variants_vs_FPKM_ESR1.py -s ../Sample_Data/ -c "chr6:151650000-152130000"
"""
#%% Imports

import argparse
import glob
import re
import gzip
import matplotlib.pylab as plt
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

"For 1 gene, but all samples. Just var over fpkm. no length correction needed."

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
    tsv_sample_name=file.split("/")[-5]
    with open(file, "r") as tsv:
        for line in tsv:
            #if line starts with E then its a gene id
            if line.split("\t")[1]=="ESR1":
                fpkm[tsv_sample_name]=float(line.strip("\n").split("\t")[6])
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
sample_names=[]
#Progress update
total_files=len(vcf_file_list)
current_number_files=0
percentage=100*(current_number_files/total_files)
print("Filtering vcf {:.2f}%".format(percentage), end="\r")
variant_counts=dict()
for file in vcf_file_list:
    sample_name=file.split("/")[-5]
    if sample_name not in sample_names:
        sample_names.append(sample_name)
    variant_counts[sample_name]=0
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
        

        "Count the entry, for those who pass filtering"
        variant_counts[sample_name]+=1
        
        
        
    current_number_files+=1
    percentage=100*(current_number_files/total_files)
    print("Filtering vcf: {:.2f}%".format(percentage),end="\r")

print("Filtering vcf: Done!            \n",end="\r")


x=[]
y=[]
coordinates=[]
#for sample in sample_names:
#    coordinates.append((sample, variant_counts[sample], fpkm[sample]))
#for x in coordinates: 
#    plt.annotate(x[0], (x[1], x[2]))
for sample in sample_names:    
    x.append(math.log2(fpkm[sample]+1))
    y.append(variant_counts[sample])

plt.plot(x,y,"ro")
plt.xlabel("FPKM")
plt.ylabel("Variant Counts")
plt.title("Variant Counts/FPKM for ESR1")
plt.savefig("Variants_vs_FPKM_ESR1.png")

