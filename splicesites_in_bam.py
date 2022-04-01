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
import subprocess

parser = argparse.ArgumentParser(prog='Find Splicesites',
                                 usage='%(prog)s -s SAM-INPUT -o OUTPUT ',
                                 description='Returns coordinates of unique splicesites in SAM file.')

parser.add_argument('--bamfile', '-b', 
                    help='BAM file, containing alignment information.')

parser.add_argument('--gff', '-g',
                    help='GFF3 file, containing annotation for genome in question.')
parser.add_argument('--out', '-o',
                    help='Output bed file, where comparison data should be printed to.')

args=parser.parse_args()


#Extract information from SAM:
format_splicesite= re.compile(r'\d+M\d+N\d+M')

#%%
"Functions"

def which_strand(flag):
    binary_flag=bin(int(flag))
    if str(binary_flag)[-5]=="1":
        return "-"
    else:
        return "+"


#%%

"1. Find splice sites in SAM file, save their information in dictionary"
splicesite=dict()
if args.bamfile:
    command=['samtools','view','-h','-F','0x100', args.bamfile,'chr6:151625000-152103274']
    p=subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output=p.communicate(input=args.bamfile.encode())[0]
    sam=output.decode('utf-8')
    sam=sam.split("\n")
    for line in sam:
        if line.startswith("@"):
            continue
        else:
            if bool(re.search(r"(\d+)\t([a-z]{3}\d+)\t(\d+)\t(\d+)\t(\d+M\d+N\d+M)\t.+YT:Z:([A-Z]{2})", line))==True:
                entries=re.search(r"(\d+)\t([a-z]{3}\d+)\t(\d+)\t(\d+)\t(\d+M\d+N\d+M)\t.+YT:Z:([A-Z]{2})",line)
                #reads=which_strand(entries.group(1)) #not needed, determines if read is reverse or not.
                chrom=entries.group(2)
                start=entries.group(3)
                cigar=entries.group(5)
                pair=entries.group(6)
                #Only include concordant pairs:
                if pair=="CP":
                    #Find coordinates for junction:
                    current_cigar=cigar
                    current_length=0
                    #To allow for several junctions in one cigar string.
                    while bool(re.search(format_splicesite,current_cigar))==True:
                        #Only process part of cigar string which matches the format (softclips removed)
                        junction=re.search(r'(\d+)M(\d+)N(\d+)M',current_cigar)
                        exon1=int(junction.group(1))
                        intron=int(junction.group(2))
                        exon2=int(junction.group(3))
                        #Exclude junctions which have less than 3 bp on either side:
                        if exon1>2 and exon2>2:
                            exon1_end=int(start)+exon1+current_length-1
                            exon2_start=exon1_end+intron+1
                            #If the whole junction is in range (only works for ESR1 atm)
                            if exon1_end>151625000 and exon2_start<152103274:
                                #Generate data for bed file
                                score=0
                                RGB="RGB"
                                exon_count="2,2"
                                block_sizes="tbd"
                                name=str(exon1_end)+"_"+str(exon2_start)
                                #The plus for strand is given rn due to it running only for ESR1. Adjust later with "reads", l68
                                splicesite[name]=[chrom, exon1_end-1, exon2_start, score, "+",
                                                  RGB, exon_count, block_sizes, "0,"+str(exon2_start-exon1_end)]
                        #remove used part of the string
                        current_cigar=re.sub(r'^.*?N','N', current_cigar).lstrip("N")
                        current_length+=exon2_start



"Write the results in BED file"
with open(args.out, 'w') as out:
    out.write("\n\n")
    counter=0
    for name in sorted(splicesite):
        #estrogen is on + strand, so we exclude - for now.
        #if splicesite[name][4]=="+":
        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            splicesite[name][0],splicesite[name][1], splicesite[name][2], name, 
            splicesite[name][3], splicesite[name][4], splicesite[name][1], splicesite[name][2],
            splicesite[name][5], splicesite[name][6], splicesite[name][7], splicesite[name][8]))
        counter+=1

        
                
            
                                
