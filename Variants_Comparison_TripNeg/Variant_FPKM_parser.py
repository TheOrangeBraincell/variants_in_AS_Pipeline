# -*- coding: utf-8 -*-
"""
date: 05-10-22
title: Variant_FPKM_parser.py
author: Mirjam Karlsson-Müller

Description:
    Reads a table with variant information, checks what gene the variant lies in,
    if any, and then retrieves the fpkm value for that gene in that sample.
    Returns table with an additional FPKM column.
    
Usage:
    python -s ../Sample_Data/ -i rna_variants.txt -o rna_variants_FPKM.tsv -c ../Database/chromosome_sizes.txt
"""

import argparse
import time
#%% Time

start_time=time.time()

#%%

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Variants FPKM Parser',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     -i INPUT',
                                 description="""Creates a genotype table out of
                                 several samples. Containing location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='file containing sample folder names containing among \
                        others, the vcf, bam and gene.tsv files.')
parser.add_argument('--input', '-i', type=str,
                    help="""Input Table""")
parser.add_argument('--output', '-o', type=str,
                    help="""Output Table""")
parser.add_argument('--chromosome', '-c', type=str,
                    help="""GENCODE chromosome sizes""")


args = parser.parse_args()

#%% 1. Chromosome Dictionary

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

#%% 2. Read genes and FPKMs from gene.tsv.













#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time)) 