# -*- coding: utf-8 -*-
"""
date: 29.03.22
name: splicesites_in_bam.py
author: Mirjam Karlsson-Müller


"""




#input: BAM file, GFF3 file.
#output: file with common splicesites, and unique ones to each source.

#%%
"1. For now: Take sam file, find unique splice sites:"

import argparse
import re

parser = argparse.ArgumentParser(prog='Find Splicesites',
                                 usage='%(prog)s -s SAM-INPUT -o OUTPUT ',
                                 description='Returns coordinates of unique splicesites in SAM file.')

parser.add_argument('--samfile', '-s', 
                    help='SAM file, containing alignment information.')

parser.add_argument('--gff', '-g',
                    help='GFF3 file, containing annotation for genome in question.')
parser.add_argument('--out', '-o',
                    help='Output file, where comparison data should be printed to.')

args=parser.parse_args()


#Extract information from SAM:

format_splicesite= re.compile('\d+M\d+N\d+M')

splicesite=dict()
if args.samfile:
    with open(args.samfile, 'r') as sam:
        for line in sam:
            if line.startswith("@"):
                continue
            else:
                entry=line.strip().split("\t")
                if bool(re.search(format_splicesite, entry[5]))==True:
                    #key=chromosome name, value=[coordinates first exon, coordinates second exon, length intron]
                    coordinates=entry[5].split("M")
                    exon1=coordinates[0]
                    intron_length=coordinates[1].split("N")[0]
                    exon2=coordinates[1].split("N")[1]
                    splicesite[entry[5]]=[[entry[3], int(exon1)+int(entry[3])], intron_length, [int(entry[3])+int(exon1)+int(intron_length),int(exon2)+int(entry[3])+int(exon1)+int(intron_length)]]

print(splicesite)


#Find corresponding exons in GFF3:

exons=[]
if args.gff:
    with open(args.gff, 'r') as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            else:
                entry=line.strip().split("\t")
                if entry[2]=="exon":
                    if entry[8].split("gene_name=")[1].split(";")[0]=="DDX11L1":
                        exons.append([entry[3],entry[4]])
                    
#sort exons from gff3, to see intron length:
exons_s=sorted(exons)
for i in range(1, len(exons_s)-1):
    print(int(exons_s[i][0])-int(exons_s[i-1][1]))

print(exons_s)
                    
                
            
                                
