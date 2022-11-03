# -*- coding: utf-8 -*-
"""
Date: Thu Nov  3 13:01:32 2022
File Name: vcf_location_table.py
Author: Mirjam Karlsson-Müller

Description:
    Parses through all variant calling files from a cohort and makes a list of all
    variants found (filtered and unfiltered).
    Outputs a table sample x location.
    
List of Functions:
    
Procedure: 
    1. Argparse
    2. Read in variant calling files, generate table (incomplete) of all found locations/variants.
    3. Print said table into a file.
    
Useage: Run in command line
    
    With coordinates for f.e. estrogen receptor:
        python vcf_location_table.py -s ../Sample_Data/ -o location_table_ESR1.tsv -c "chr6:151690496-152103274"
    
    For whole genome:
        python vcf_location_table.py -s ../Sample_Data/ -o location_table_WG.tsv
    
    Server:
        python vcf_location_table.py -s /raidset/mi3258mu-se/Mirjam -o location_table_WG.tsv
    
Possible Bugs:
    
    
"""

#%% Imports
import argparse
import glob
import gzip
import os
import re
import time

#%% Start Timer

start_time=time.time()


#%% Argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='vcf location table',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     [-c] "chrX:XXXXXX-XXXXXX"',
                                 description="""Creates a location of variants table out of
                                 several samples. Containing location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing the vcf files.')
parser.add_argument('--out', '-o', required=True,
                    help="""Output file, containing vcf location table.""")
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


#Check if output file already exists.
while True:
    if os.path.isfile(args.out):
        #If it does, ask user if they want to overwrite the output.
        answer=str(input("The output file "+ args.out+ " already exists. Would you like to overwrite it? (Y/N) "))
        if answer.upper()=="Y":
            break
        elif answer.upper()=="N":
            new_output=str(input("Please enter a new output name: "))
            if all([e in [" ", "", "\t", "\n"] for e in new_output]):
                print("This is not a valid output name. Lets try this again.")
            else:
                output=new_output
        else:
            print("This is not a valid response. Please answer with Y (yes) or N (no).")
    else:
        #There is no problem with this output file name. We proceed with the code.
        break

#%% Read vcf file

#Find all the vcf files in the data folder.
argument_glob=args.samples+"/**/*.vcf.gz"
vcf_file_list=glob.glob(argument_glob, recursive=True)

#If no vcf files are found, quit the program.
if len(vcf_file_list)==0:
    print("""There were no vcf files found in the input folder. Please make 
          sure to use the right input folder. The vcf files can be in any 
          subfolder of the input folder.""")
    quit()

#Initialize dictionary
variants=dict()
#Sample list
sample_names=[]
#Progress update
total_files=len(vcf_file_list)
current_number_files=0
percentage=100*(current_number_files/total_files)
print("Reading vcf {:.2f}%".format(percentage), end="\r")
for file in vcf_file_list:
    vcf = gzip.open(file, "rt")
    sample_name=re.search(r"/(S\d+)/",file).group(1)
    if sample_name not in sample_names:    
        sample_names.append(sample_name)
    #read through entries
    for line in vcf:
        #Skip headers
        if line.startswith("#"):
            continue
        chrom, position, ID, ref, alt, qual, filt, info, form, sample =line.split("\t")
        genotype=sample.split(":")[-7].strip(" ")
        
        #If coordinates or name are given, restrict area of variants.
        if args.coordinates:
            if chrom != coord_chrom:
                continue
            else:
                if int(position)< int(coord_start) or int(position)> int(coord_stop):
                    continue
        "Before filtering, add all variants to genotype table."
        variant_ID=chrom+"_"+position
        if variant_ID not in variants:
            variants[variant_ID]=dict()
            
        variants[variant_ID][sample_name]=[chrom, position, "ND"]
        
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
        
        if genotype=="0/1" or genotype=="1/0":    
            variants[variant_ID][sample_name][2]= "0/1"
        elif genotype=="1/1":
            variants[variant_ID][sample_name][2]= "1/1"
        #although it should not happen....  But it clearly is.
        elif genotype=="0/0":
            variants[variant_ID][sample_name][2]= "0/0"
        else:
            print("Invalid genotype ", genotype)
    
    "Progress updates on number of vcf files."
    current_number_files+=1
    percentage=100*(current_number_files/total_files)
    print("Reading vcf: {:.2f}%".format(percentage),end="\r")
        
        
print("Reading vcf: Done!            \n",end="\r")

sample_names=sorted(sample_names)


#%% Write output file

counter=0
percentage=100*counter/len(list(variants.keys()))
print("Writing Location Table: {:.2f}%".format(percentage),end="\r")

with open(args.out, "w") as out:
    #Create Header
    if args.coordinates:
        out.write("#Variant Location Table for variants in range "+ args.coordinates + "\n")
    else:
        out.write("#Variant Location Table for variants in whole genome.\n")
    out.write("Location\t"+"\t".join(sample_names)+"\n")
    
    for variant in variants:
        #Starting string for file:
        new_line=[variant]
        #Go through sorted samples, to always have the same order.
        for sample in sample_names:
            if sample in variants[variant]:
                new_line.append(variants[variant][sample][2])
            else:
                #Put a spaceholder. Genotype will be determined using gene expression data by next script in pipeline.
                new_line.append("-")
        out.write("\t".join(new_line)+"\n")
        counter+=1
        percentage=100*counter/len(list(variants.keys()))
        print("Writing Location Table: {:.2f}%".format(percentage),end="\r")

print("Writing Location Table: Done          \n",end="\r")
        


#%% Stop Timer

print("Run time: {:.2f} seconds.".format(time.time()-start_time))      
