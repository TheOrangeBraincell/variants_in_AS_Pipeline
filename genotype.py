# -*- coding: utf-8 -*-
"""
Date: Thu Jan 19 10:24:12 2023
File Name: genotype.py
Author: Mirjam Karlsson-Müller

Description:
    Assigns genotypes NE (not expressed) or 0/0 (homozygous reference) to location table entries marked with "-"
    
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
    python variants_in_AS_Pipeline/genotype.py -f Database/fpkm_table.tsv -i location_table_ESR1.tsv -o genotype_table_ESR1.tsv -r Database/gene_ranges.tsv -c "chr6:151656691-152129619"
    
    
Possible Bugs:
    location tables need to be created by vcf_location_table.py
    gene_ranges.tsv needs to be created with gene_ranges.py
"""

#%% Imports
import argparse
import os
import time
import re

#%% Start Timer 

start_time=time.time()

#%% Argparse

parser = argparse.ArgumentParser(prog='assigning genotypes',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     -i LOCATION-TABLE [-c] "chrX:XXXXXX-XXXXXX"\
                                         -r RANGE-TSV -f FPKM-TABLE',
                                 description="""Creates a genotype table out of
                                 a location table. Containing genotypes for location x sample.""")

parser.add_argument('--fpkm', '-f', required=True,
                    help='table containing fpkm values for all genes and samples')
parser.add_argument('--input', '-i', required=True,
                    help="""Input file containing location table.""")
parser.add_argument('--out', '-o', required=True,
                    help="""Output file containing genotypes.""")
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")
parser.add_argument('--ranges','-r', required=True, help= "File containing \
                    gene ranges created with refseq and gencode.")
                    

args = parser.parse_args()


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
""" Manual option: Interactive at console
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
"""

#Server option: Make sure no outputs are accidentally deleted.
if os.path.isfile(args.out):
    print("The output file already exists, we assume it is complete.")
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
        
        
print("Gene Ranges dictionary made: {:.2f} seconds.".format(time.time()-start_time))
#%% Read in fpkm table

fpkms=dict()
with open(args.fpkm, "r") as fpkm:
    for line in fpkm:
        if line.startswith("Location"):
            #header
            sample_names=line.strip("\n").split("\t")[1:]
            continue
        #print(sample_names.index("S001298"))        
        #only gene names in range of input coordinates
        if line.strip("\n").split("\t")[0] not in gene_ranges:
            continue
        #That should cut down memory requirement and speed things up further down.
        #add fpkm values to dictionary.
        values=line.strip("\n").split("\t")[1:]
        #initialize dict
        fpkms[line.split("\t")[0]]=dict()
        for sample in sample_names:
            #if sample=="S001298":
                #print("The sample is in the header of the fpkm table")
            fpkms[line.split("\t")[0]][sample]=values[sample_names.index(sample)]
        
print("FPKM file is read in: {:.2f} seconds.".format(time.time()-start_time))

#%% Open location file and output genotype table

locations=open(args.input, "r")
genotypes=open(args.out, "w")

if args.coordinates:
    #Write header for genotypes table
    genotypes.write("#Genotype table per sample for variants in range " + args.coordinates + "\n")
else:
    genotypes.write("#Genotype table per sample for variants in whole genome range\n")

print("Opened input and output file: {:.2f} seconds.".format(time.time()-start_time))

#%% Assign Genotypes

genes=list(gene_ranges.keys())
total=len(genes)
count=0
percentage=100*(count/total)
print("Genotyping {:.2f}%".format(percentage), end="\r")

#Iterate through lines in location table as well as genes simultaneously..
for line in locations:
    if line.startswith("#"):
        genotypes.write(line)
        continue
    if line.startswith("Location"):
        #Thats the column names. We write those and extract sample names.
        sample_names=line.strip("\n").split("\t")[1:]
        #Write header for output file + column for gene name!
        genotypes.write(line.split("\t")[0]+ "\tGene\t"+"\t".join(sample_names)+"\n")
        
        #Initiate gene walk through, requires sample names.
        #initiate gene counter
        current_index=0
        current_gene=genes[current_index]
        current_range=gene_ranges[current_gene]
        #Make fpkm dict for this current gene:
        fpkm=dict()
        #Do we have expression values for this gene?
        if current_gene not in fpkms:
            #nope. All NE.
            for sample in sample_names:
                fpkm[sample]="NE"
        else:
            #go through samples
            for sample in sample_names:
                if sample in fpkms[current_gene]:
                    #if sample=="S001298":
                    #    if sample in fpkms[current_gene]:
                    #        print(fpkms[current_gene]["S001298"])
                    #save genotype dictionary if fpkm >=10 or not
                    if float(fpkms[current_gene][sample])>=10:
                        fpkm[sample]="0/0"
                    else:
                        fpkm[sample]="NE"
                else:
                    fpkm[sample]="NE"
        #print(current_gene, fpkm)
        continue
    #Now all that is left is entries. shape: chr_position_ref_(alt)\t samples
    #print(line.split("\t")[0].split("_")[0:2])
    variant_chrom,variant_position = line.split("\t")[0].split("_")[0:2]
    #infostring to use in first column of output.
    infostring=line.split("\t")[0]
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
            out_gene="intergenic"
            #We found the variant, so we can go to next line
            variant_found=True
        
        #if the variant position is in the current gene range, use fpkm dict to assign right genotype
        elif int(variant_position)>= int(current_range[2]) and int(variant_position)<=int(current_range[3]):
            #go through entries
            for i in range(0, len(entries)):
                if entries[i]=="-":
                    sample= sample_names[i]
                    #assign genotype from fpkm dict for corresponding sample
                    if sample in fpkm:
                        entries[i]=fpkm[sample]
                    else:
                        entries[i]="NE"
            out_gene=current_gene
            #We found the variant, so we can go to next line
            variant_found=True
        
        #else the variant is after current gene. Then we need to update the gene and check again.
        else:          
            #update index:
            current_index+=1
            #Update on progress
            count+=1
            percentage=100*(count/total)
            print("Genotyping {:.2f}%".format(percentage), end="\r")
            #If last gene, stop
            if current_index>=len(genes):
                #Then the rest of the entries will be after and thus outside of a gene.
                for i in range(0, len(entries)):
                    if entries[i]=="-":
                        entries[i]="NE"
                out_gene="intergenic"
                variant_found=True
            else:
                #update gene and range
                current_gene=genes[current_index]
                current_range=gene_ranges[current_gene]
                #Make fpkm dict for this current gene:
                fpkm=dict()
                #Do we have expression values for this gene?
                if current_gene not in fpkms:
                    #nope. All NE.
                    for sample in sample_names:
                        fpkm[sample]="NE"
                else:
                    #go through samples
                    for sample in sample_names:
                        if sample in fpkms[current_gene]:
                            #if sample=="S001298":
                            #    if sample in fpkms[current_gene]:
                            #        print(fpkms[current_gene]["S001298"])
                            #save genotype dictionary if fpkm >=10 or not
                            if float(fpkms[current_gene][sample])>=10:
                                fpkm[sample]="0/0"
                            else:
                                fpkm[sample]="NE"
                        else:
                            fpkm[sample]="NE"
            #do not go to next line, as we did not find the variants location in regards to the gene.
    
    #Now that we have replaced all "-" in the line, we write the line out to the output file.
    #print(current_gene, entries)
    if variant_found==True:
        genotypes.write(infostring+"\t"+out_gene+"\t"+"\t".join(entries)+"\n")


print("Genotyping Done!                 \n", end="\r")


#%% Close location file and output file.

locations.close()
genotypes.close()

print("Files closed!")

#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
