# -*- coding: utf-8 -*-
"""
Date: Thu Jan 19 10:24:12 2023
File Name: genotype.py
Author: Mirjam Karlsson-Müller

Description:
    Assigns genotypes NE or 0/0 to location table entries marked with "-"
    
List of Functions:
    none
    
Procedure: 
    1. Reads in gene range file created by gene_ranges.py
    2. Iterates through genes and location table at the same time (as they are sorted by coordinates thats ok)
    3. If variant in gene, assign genotype based on fpkm (>=10 0/0, otherwise NE), otherwise assign genotype NE (not expressed)

Input: 
    - Sorted Location table as created by vcf_location_table.py and sorted by Sort_Locations.sh
    - sample folder containing gene.tsv


Useage:
    for ESR1 f.e.
    python genotype.py -s ../Sample_Data/ -i location_table_ESR1.tsv -o genotype_table_ESR1.tsv -r gene_ranges.tsv -c "chr6:151690496-152103274" 
    
    
Possible Bugs:
    location tables need to be created by vcf_location_table.py and sorted with Sort_Location.sh
    gene_ranges.tsv needs to be created with gene_ranges.py
"""

#%% Imports
import argparse
import glob
import os
import time
import re

#%% Start Timer 

start_time=time.time()

#%% Argparse

parser = argparse.ArgumentParser(prog='assigning genotypes',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     -i LOCATION-TABLE [-c] "chrX:XXXXXX-XXXXXX"\
                                         -r RANGE-TSV',
                                 description="""Creates a genotype table out of
                                 a location table. Containing genotypes for location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing gene.tsv files.')
parser.add_argument('--input', '-i', required=True,
                    help="""Input file containing sorted location table.""")
parser.add_argument('--out', '-o', required=True,
                    help="""Output file containing incomplete genotype table.""")
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")
parser.add_argument('--ranges','-r', required=True, help= "File containing \
                    gene ranges created with refseq and gencode.")
                    

args = parser.parse_args()

#Check if input file exists:
if not os.path.exists(args.samples):
    print("The input folder is invalid. Did you misstype the directory structure? Try again.")
    quit()

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
""" Manual option. Need server option.
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

#Server option
if os.path.isfile(args.out):
    print("The output file already exists, we assume it is complete.")
    quit()
"""    
    
#%% List of gene.tsv's

#Find all the gene.tsv files in the data folder.
argument_glob=args.samples+"/**/gene.tsv"
tsv_list=glob.glob(argument_glob, recursive=True)

#If no tsv files are found, quit the program.
if len(tsv_list)==0:
    print("""There were no gene.tsv files found in the input folder. Please make 
          sure to use the right input folder. The gene.tsv files can be in any 
          subfolder of the input folder.""")
    quit()

#%% Dictionary of gene ranges from gene range file (made with gencode and refseq) 

gene_ranges=dict()

with open(args.ranges, "r") as ranges:
    for line in ranges:
        #skip header line
        if line.startswith("Gene"):
            continue
        #Input file is of format Gene\tChrom\tStrand\tMin\tMax
        #Only save information within input coordinates
        if args.coordinates:
            if line.split("\t")[1]!=coord_chrom:
                continue
            elif int(line.split("\t")[4])<coord_start or int(line.split("\t")[3])>coord_stop:
                continue
        #Save information in dictionary
        gene_ranges[line.split("\t")[0]]=line.strip("\n").split("\t")[1:]

#%% Open location file and output genotype table

locations=open(args.input, "r")
genotypes=open(args.out, "r+")

#List of locations of genotype table
already_exists=[]
for line in genotypes:
    if line.startswith("#"):
        continue
    if line.startswith("Location"):
        continue
    
    already_exists.append(line.split("\t")[0])

"""
#Write header for genotypes table
genotypes.write("#Genotype table per sample for variants in range " + args.coordinates + "\n")
"""
#%% Assign Genotypes

genes=list(gene_ranges.keys())

#First new line
new_line=False
#Iterate through lines in location table as well as genes simultaneously..
for line in locations:
    if line.startswith("#"):
        continue
    elif line.startswith("Location"):
        #Thats the column names. We write those and extract sample names.
        sample_names=line.strip("\n").split("\t")[1::]
        continue
    
    #if line already in genotypes, skip
    elif line.split("\t")[0] in already_exists:
        continue
    else:
        #When the first new line appears, we need to initialize current gene, range and index.
        if new_line==False:
            variant_position=int(line.split("\t")[0].split("_")[1])
            #Then we need to find out where we are on the chromosome
            current_index=0
            #initialize variable previous_stop
            previous_stop=0
            for gene in gene_ranges:
                start_gene=int(gene_ranges[gene][2])
                stop_gene=int(gene_ranges[gene][3])
                
                #If variant is in a gene, then this will find it. Otherwise not. Crashes.
                if start_gene<variant_position and stop_gene>variant_position:
                    #current gene found
                    current_gene=gene
                    current_range=gene_ranges[current_gene]
                    #Initiate fpkm dict for this gene
                    #Make fpkm dict for this current gene:
                    fpkm=dict()
                    #go through samples
                    for sample in sample_names:
                        for tsv in tsv_list:
                            if sample in tsv: 
                                #Find corresponding gene, save for sample name
                                file=open(tsv, "r")
                                #If the gene name is not in the gene.tsv file, we assume not expressed due to no information
                                gene_found=False
                                for l in file:
                                    if current_gene in l:
                                        gene_found=True
                                        #save genotype dictionary if fpkm >=10 or not
                                        if float(l.split("\t")[7])>=10:
                                            fpkm[sample]="0/0"
                                        else:
                                            fpkm[sample]="NE"
                                        break
                                if gene_found==False:
                                    for sample in sample_names:
                                        fpkm[sample]="NE"
                                file.close()
                    #If for some sample, there was no match of the name found, notice here:
                    for sample in sample_names:
                        if sample not in fpkm:
                            fpkm[sample]=="NE"
                    break
                #Check if its between this gene and the previous one
                elif start_gene > variant_position and previous_stop< variant_position:
                    #In between genes. 
                    current_gene=gene
                    current_range=gene_ranges[current_gene]
                    #But we still need fpkms for this current gene, because the next variant could be in it.
                    #Initiate fpkm dict for this gene
                    #Make fpkm dict for this current gene:
                    fpkm=dict()
                    #go through samples
                    for sample in sample_names:
                        for tsv in tsv_list:
                            if sample in tsv: 
                                #Find corresponding gene, save for sample name
                                file=open(tsv, "r")
                                #If the gene name is not in the gene.tsv file, we assume not expressed due to no information
                                gene_found=False
                                for l in file:
                                    if current_gene in l:
                                        gene_found=True
                                        #save genotype dictionary if fpkm >=10 or not
                                        if float(l.split("\t")[7])>=10:
                                            fpkm[sample]="0/0"
                                        else:
                                            fpkm[sample]="NE"
                                        break
                                if gene_found==False:
                                    for sample in sample_names:
                                        fpkm[sample]="NE"
                                file.close()
                    #If for some sample, there was no match of the name found, notice here:
                    for sample in sample_names:
                        if sample not in fpkm:
                            fpkm[sample]=="NE"
                    break
                else:
                    previous_stop=stop_gene
                    current_index+=1
            
            new_line=True
        
    
        #Now all that is left is entries. shape: chr_position\t samples
        variant_chrom,variant_position = line.split("\t")[0].split("_")
        entries=line.strip("\n").split("\t")[1:]
        
        #print(current_gene)
        #But first dummy check if we are on the right chromosome:
        if variant_chrom!= current_range[0]:
            print("The range of interest when starting the genotype script does not match the range of interest of your location table. Please adjust.")
            quit()
        #We try and find where the position is regarding genes. Until we find it, we keep repeating it.
        variant_found=False
        while variant_found==False:
            #The variant locations are saved in increasing order. So if the variant is before current gene start
            #then it is either before the first gene or between genes. either way genotype is NE.
            if int(variant_position) < int(current_range[2]):
                #then we can replace all "-" entries with NE as they are outside of a gene.
                for i in range(0, len(entries)):
                    if entries[i]=="-":
                        entries[i]="NE"
                #We found the variant, so we can go to next line
                variant_found=True
            
            #if the variant position is in the current gene range, use fpkm dict to assign right genotype
            elif int(variant_position)>= int(current_range[2]) and int(variant_position)<=int(current_range[3]):
                #go through entries
                for i in range(0, len(entries)):
                    if entries[i]=="-":
                        sample= sample_names[i]
                        #assign genotype from fpkm dict for corresponding sample
                        entries[i]=fpkm[sample]
                #We found the variant, so we can go to next line
                variant_found=True
            
            #else the variant is after current gene. Then we need to update the gene and check again.
            else:          
                #update index:
                current_index+=1
                #If last gene, stop
                if current_index>=len(genes):
                    #Then the rest of the entries will be after and thus outside of a gene.
                    for i in range(0, len(entries)):
                        if entries[i]=="-":
                            entries[i]="NE"
                    variant_found=True
                else:
                    #update gene and range
                    current_gene=genes[current_index]
                    current_range=gene_ranges[current_gene]
                    #Make fpkm dict for this current gene:
                    fpkm=dict()
                    #go through samples
                    for sample in sample_names:
                        for tsv in tsv_list:
                            if sample in tsv: 
                                #Find corresponding gene, save for sample name
                                file=open(tsv, "r")
                                #If the gene name is not in the gene.tsv file, we assume not expressed due to no information
                                gene_found=False
                                for l in file:
                                    #wrong chrom skip
                                    if variant_chrom!=l.split("\t")[2]:
                                        continue
                                    
                                    if current_gene in l:
                                        gene_found=True
                                        #save genotype dictionary if fpkm >=10 or not
                                        if float(l.split("\t")[7])>=10:
                                            fpkm[sample]="0/0"
                                        else:
                                            fpkm[sample]="NE"
                                        break
                                if gene_found==False:
                                    for sample in sample_names:
                                        fpkm[sample]="NE"
                                file.close()
                    #If for some sample, there was no match of the name found, notice here:
                    for sample in sample_names:
                        if sample not in fpkm:
                            fpkm[sample]=="NE"
                #do not go to next line, as we did not find the variants location in regards to the gene.
        
        #Now that we have replaced all "-" in the line, we write the line out to the output file.
        genotypes.write(variant_chrom+"_"+variant_position+"\t"+"\t".join(entries)+"\n")


#%% Close location file and output file.

locations.close()
genotypes.close()

#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    