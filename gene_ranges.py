# -*- coding: utf-8 -*-
"""
Date: Thu Jan 19 10:28:18 2023
File Name: gene_ranges.py
Author: Mirjam Karlsson-Müller

Description:
    Produces tsv file containing gene ranges.
    
List of Functions:
    none

Procedure: 
    1. Reads in database GENCODE and REFSEQ annotation files
    2. Extracts exons
    3. Finds smallest start and biggest stop, saves as gene range
    
Useage:
    #with coordinates f.e. Estrogen Receptor
    python variants_in_AS_Pipeline/gene_ranges.py -o gene_ranges.tsv -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv
    
Possible Bugs:
"""
#%% Imports

import argparse
import glob
import re
import time
import os
import math


#%% Time

start_time=time.time()


#%% argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Create gene range tsv',
                                 usage='%(prog)s -o OUTPUT-FILE \
                                     -g GENCODE-FILE -r REFSEQ-FILE \
                                         [-c] "chrX:XXXXXX-XXXXXX" -as AS-TYPE \
                                             -is INSERT-SIZE',
                                 description="""Per AS event of interest, creates
                                 a table with the PSI scores supporting said event
                                 per sample in sample folder.""")

parser.add_argument('--out', '-o', required=True,
                    help="""Output file, containing gene ranges.""")
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")
parser.add_argument('--gencode', '-g', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from GENCODE39 as well as gene names.""")
parser.add_argument('--refseq', '-r', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from RefSeq as well as gene names.""")


args = parser.parse_args()

# Extract input coordinates, check their format.
if args.coordinates:
    if re.search(r'[a-z]{3}[MXY]?\d*:\d+-\d+', args.coordinates):
        coord = re.search(r'([a-z]{3}[MXY]?\d*):(\d+)-(\d+)', args.coordinates)
        coord_chrom = coord.group(1)
        #To adjust for different inclusivity of stop/start for plus and minus strand, expand coordinate range by 1.
        coord_start = int(coord.group(2))-1
        coord_stop = int(coord.group(3))+1

    else:
        raise argparse.ArgumentTypeError("""The given coordinates are in the 
                                         wrong format. Input as 
                                         chrX:XXXX-XXXX.""")
        quit()


#%% 1. Process Database input

#outer dictionary with gene names as key
gene_dict=dict()
print("Creating Database Dictionary...", end="\r")

for file in [args.gencode, args.refseq]:
    with open(file, "r") as infile:
        for line in infile:
            # To exclude potential title lines/empty lines, formatting mistakes
            # Only takes chr[] and chr[]_random lines, in accordance with bam.
            if re.search(r"(.+)\t(?:([a-z]{3}[X,M,Y]?\d*)|([a-z]{3}[X,M,Y]?\d*)"
                         r".+_random)\t(\-?\+?)\t(\d+)\t(\d+)\t(\d+)\t([\d,]+)"
                         r"\t([\d,]+)\t(.+)", line):
                # specify groups.
                entry = re.search(r"(.+)\t([a-z]{3}[X,M,Y]?\d*).*\t"
                             r"(\-?\+?)\t(\d+)\t(\d+)\t(\d+)\t([\d,]+)"
                             r"\t([\d,]+)\t(.+)", line)
                #Assign variables to groups
                trans_ID=entry.group(1)
                chrom=entry.group(2)
                strand=entry.group(3)
                
                """To not get caught in start/stop, -/+ strand complications,
                coordinates are referred to as bigger and smaller instead."""
                
                number_exons=int(entry.group(6))
                exon_smaller=entry.group(7).split(",")[0:-1]
                exon_bigger=entry.group(8).split(",")[0:-1]
                gene_name=entry.group(9)
                if file==args.gencode:
                    db="G"
                elif file==args.refseq:
                    db="R"
                else:
                    db="You got a bug."
                
                #If choordinates are given:
                if args.coordinates:
                    #exclude entries outside of coordinates
                    if int(exon_smaller[0])<coord_start or \
                        int(exon_bigger[-1])>coord_stop or chrom!=coord_chrom:
                        continue
                
                #if its a new gene symbol, initialize inner dictionary
                if gene_name not in list(gene_dict.keys()):
                    gene_dict[gene_name]={trans_ID:[]}
                    
                else:
                    #add transcript ID to gene_dict dictionary
                    gene_dict[gene_name][trans_ID]=[]
                
                
                #make entries for each exon.
                for i in range(0, number_exons):
                    if i==0:
                        if strand=="+":
                            position= "first"
                        else:
                            position="last"
                    elif i==number_exons-1:
                        if strand=="+":
                            position="last"
                        else:
                            position="first"
                    else:
                        position="middle"
                    gene_dict[gene_name][trans_ID].append([chrom, 
                                                           exon_smaller[i], 
                                                           exon_bigger[i], 
                                                           strand, position,db])

print("Creating Database Dictionary: Done! \n", end="\r")

#%%
#Make gene ranges dictionary to use in other parts of analysis
gene_ranges=dict()
for gene in gene_dict:
    starts=[]
    stops=[]
    for trans_ID in gene_dict[gene]:
        starts.append(int(gene_dict[gene][trans_ID][0][1]))
        stops.append(int(gene_dict[gene][trans_ID][-1][2]))
        chrom=gene_dict[gene][trans_ID][0][0]
        strand=gene_dict[gene][trans_ID][0][3]
    #Take smallest start and biggest stop as range of gene.
    gene_ranges[gene]=[chrom, strand, min(starts), max(stops)]
    
    
    
    
#%% Write gene ranges into output file.

with open(args.out, "w") as out:
    out.write("Gene\tStrand\tstart\tstop\n")
    for gene in gene_ranges:
        out.write(gene+"\t"+"\t".join([str(i) for i in gene_ranges[gene]])+"\n")

#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  

