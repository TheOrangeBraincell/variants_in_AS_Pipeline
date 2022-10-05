# -*- coding: utf-8 -*-
"""
name: Parse_Lifted_Bed.py
author: Mirjam Karlsson-Müller
date: 30.9.22

Description: Reformats the output bed from LiftOver of ucsc.

Usage:
    python Parse_Lifted_Bed.py -v hglft_genome_34d4d_aad970.bed -db ../Database/gencode.v38.annotation.gff3.gz -c ../Database/chromosome_sizes.txt -o tnbc_variants.tsv


"""

#%% Imports
import argparse
import re
import gzip

#%% argparse

parser = argparse.ArgumentParser(prog='Parse TNBC')

parser.add_argument('--variants', '-v', required=True,
                    help='variant calling file tnbc')
parser.add_argument('--database', '-db', type=str,
                    help="""GENCODE annotated exons hg38""")
parser.add_argument('--out', '-o', required=True,
                    help="Output file")
parser.add_argument('--chromosome', '-c', type=str,
                    help="""GENCODE chromosome sizes""")


args = parser.parse_args()

#%% 1. Read gtf database file



chrom_sizes=dict()
with open(args.chromosome, "r") as chromosome:
    for line in chromosome:
        chrom, size= line.strip("\n").split("\t")
        chrom_sizes[chrom]=int(size)
        
        
#Split the sizes into areas (easier to read.)
chromosomes=dict()
for chrom in chrom_sizes:
    i=0
    while i<chrom_sizes[chrom]:
        start=i
        stop=i+1000000
        if stop>=chrom_sizes[chrom]:
            stop=chrom_sizes[chrom]
        if chrom not in chromosomes:
            chromosomes[chrom]=dict()
        
        chromosomes[chrom][str(start)+"-"+str(stop)]=dict()
        i=stop


# ready to add genes and their coord.
#%%exon coordinates.

db = gzip.open(args.database, "rt")

annotated_exons=dict()
gene_coordinates=dict()
for line in db:
    if line.startswith("#"):
        continue
    chrom, typ, start, stop, strand, info=list(map(line.split("\t").__getitem__,[0, 2, 3, 4, 6, 8]))
    gene_name=re.search(r"gene_name=(.+?);", info).group(1)
    
    if typ=="gene":
        gene_found=False
        for area in chromosomes[chrom]:
            """To account for the case where a gene overlaps the border between
            regions, we sort the genes to regions based only on their "smaller"
            coordinate"""
            if int(start)>=int(area.split("-")[0]) and int(start)<=int(area.split("-")[1]):
                chromosomes[chrom][area][gene_name]=[]
                gene_found=True
                gene_coordinates[gene_name]=[start, stop]
                
        if gene_found==False:
            print("gene ", gene_name, " has not been found on its chromosome.")
    
    if typ=="exon":
        exon_found=False
        for area in chromosomes[chrom]:
            for gene in chromosomes[chrom][area]:
                if gene==gene_name:
                    chromosomes[chrom][area][gene].append([start, stop])
                    exon_found=True
                    break
            if exon_found==True:
                break



#%% 2. Read tnbc data and check wether they are in exon etc.


with open(args.variants,"r") as bed, open(args.out, "w") as out:
    out.write("Sample\tchrom\tposition\tref\talt\tgenotype\tlocation\n")
    for line in bed:
        chrom, start, stop, name= line.strip("\n").split("\t")[0:4]
        sample, genotype, ref, alt= name.split("_")
        
        #Check whether the variant is in an exon.
        location="Intergenic"
        #We sort out unlocated entries.
        if chrom not in chromosomes:
            continue
        for area in chromosomes[chrom]:
            if int(stop)>= int(area.split("-")[0]) and int(stop)<=int(area.split("-")[1]):
                for gene in chromosomes[chrom][area]:
                    if int(stop)>= int(gene_coordinates[gene][0]) and int(stop)<= int(gene_coordinates[gene][1]):
                        location="Intron"
                        for exon in chromosomes[chrom][area][gene]:
                            if int(stop)>=int(exon[0]) and int(stop)<=int(exon[1]):
                                location="Exon"
                                break
                        break
                break
        
        
        
        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, chrom, stop, 
                                                    ref, alt, genotype, location))
        
        