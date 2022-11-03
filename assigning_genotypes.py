# -*- coding: utf-8 -*-
"""
Date: Thu Nov  3 13:45:42 2022
File Name: assigning_genotypes.py
Author: Mirjam Karlsson-Müller

Description:
    Parses through all gene.tsv files and extracts FPKMs per gene & sample.
    Based on these FPKMs, assigns genotypes where they are not determined yet (marked with space holder.)
    If the FPKM>10, the genotype 0/0 is assigned. Otherwise NE (not expressed.)
    The genotypes ND (not enough data=filtered out), 0/1 and 1/1 stay as they are.
    
    
List of Functions:
    -HMZR_NOEX()
    
Procedure: 
    1. Read in gene.tsv files to collect FPKM values for different samples/genes.
    2. Go through location table. Where there is a spaceholder '-', assign genotype based on FPKM
    3. Write new completed information into output file.
    
Useage: console
    
        python assigning_genotypes.py -s ../Sample_Data/ -i location_table_WG.tsv -o genotype_table.tsv
    
    #For server:
        
        python assigning_genotypes.py -s /raidset/mi3258mu-se/Mirjam -i location_table_WG.tsv -o genotype_table.tsv
    
Possible Bugs:
"""

#%% Imports
import argparse
import glob
import os
import time

#%% Start Timer 

start_time=time.time()

#%% Argparse

parser = argparse.ArgumentParser(prog='assigning genotypes',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     [-c] "chrX:XXXXXX-XXXXXX"',
                                 description="""Creates a genotype table out of
                                 a location table. Containing genotypes for location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing gene.tsv files.')
parser.add_argument('--input', '-i', required=True,
                    help="""Input file containing location table.""")
parser.add_argument('--out', '-o', required=True,
                    help="""Output file containing genotypes.""")

args = parser.parse_args()

#Check if input file exists:
if not os.path.exists(args.samples):
    print("The input folder is invalid. Did you misstype the directory structure? Try again.")
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


#%% Read in gene.tsv's

#Find all the gene.tsv files in the data folder.
argument_glob=args.samples+"/**/gene.tsv"
tsv_list=glob.glob(argument_glob, recursive=True)

#If no vcf files are found, quit the program.
if len(tsv_list)==0:
    print("""There were no gene.tsv files found in the input folder. Please make 
          sure to use the right input folder. The gene.tsv files can be in any 
          subfolder of the input folder.""")
    quit()

#initialize dictionary
tsv_info=dict()
#Progress updates
total_files=len(tsv_list)
current_file=0
percentage=100*current_file/total_files
print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

for file in tsv_list:
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


#%% Read off location table, assign genotypes as we go along and print to output file.

def HMZR_NOEX(variant_ID, sample):
    """
    Assigns a genotype at a specific location for a specific sample, if no
    variant is found at that location.
    
    Genotype HMZR (1/1) is given if the variant lies in a gene, which has fpkm >= 10.
    
    Otherwise genotype NE is given.

    Parameters
    ----------
    variant_ID= string
                "chrom_position"
    sample= string
            sample name f.e. S00001

    Returns
    -------
    genotype

    """
    genotype="NE"
    for gene in tsv_info[sample]:
        chrom, position= variant_ID.split("_")
        #print(gene)
        #Find gene that the variant is in
        if chrom==gene[0]:
            if int(gene[1])<int(position)<int(gene[2]):
                #Check FPKM value:
                if float(gene[4])>=10:
                    #print(gene[4])
                    genotype="0/0"
                    return genotype
        
    return genotype


print("Assigning Genotypes...", end="\r")

with open(args.input, "r") as infile, open(args.out, "w") as outfile:
    variants=dict()
    for line in infile:
        #Skip header
        if line.startswith("#"):
            continue
        if line.startswith("Location"):
            sample_names=line.strip("\n").split("\t")[1::]
            #Write header for output file
            outfile.write(line)
            continue
        
        entry=line.strip("\n").split("\t")
        location=entry[0]
        for i in range(1, len(sample_names)+1):
            if entry[i]=="-":
                entry[i]=HMZR_NOEX(location, sample_names[i-1])
            
        outfile.write("\t".join(entry)+"\n")
        
print("Assigning Genotypes Done! \n", end="\r")


#%% Stop Timer

print("Run time: {:.2f} seconds.".format(time.time()-start_time))    
