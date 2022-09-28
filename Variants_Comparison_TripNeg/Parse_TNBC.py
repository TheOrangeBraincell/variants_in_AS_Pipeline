# -*- coding: utf-8 -*-
"""
date: 28.09.22
author: Mirjam Karlsson-Müller
name: Parse_TNBC.py

Description:
    Creates a table containing sample, chrom, pos and genotype for each variant entry.
    
    To be run on tnbc cohort for comparison to our variant files.
Usage:
    python Parse_TNBC.py -v ../Trip_Neg_data/subs.bulk.1639 -i ../Trip_Neg_data/PD_ID_External_IDs_WGS_GSE96058.txt -o tnbc_variants.tsv 
    
"""

#%% Imports
import argparse

#%% argparse

parser = argparse.ArgumentParser(prog='Parse TNBC')

parser.add_argument('--variants', '-v', required=True,
                    help='variant calling file tnbc')
parser.add_argument('--ids', '-i', type=str,
                    help="""Pa_ids linking to the sample names we use""")
parser.add_argument('--out', '-o', required=True,
                    help="Output file")


args = parser.parse_args()


#%% Create ID dictionary

id_dict=dict()
with open(args.ids, "r") as ids:
    for line in ids:
        PD_ID=line.split("\t")[0]
        sample=line.split("\t")[2].split(".")[0]
        id_dict[PD_ID]=sample


#%% Parse

counter=0
unmatched_ids=set()
with open(args.variants, "r") as tnbc, open(args.out, "w") as out:
    out.write("#Sample\tchrom\tpos\tgenotype\n")
    for line in tnbc:
        #skip header
        if line.startswith("#"):
            continue
        
        if line.split("\t")[1] in id_dict:
            sample=id_dict[line.split("\t")[1]]
        else:
            counter+=1
            unmatched_ids.add(line.split("\t")[1])
            continue
        chrom, pos, ref, alt=line.split("\t")[4:8]
        genotype=line.split("\t")[41]
        out.write("{}\t{}\t{}\t{}\n".format(sample, chrom, pos, genotype))        
        
        


print(counter, unmatched_ids)