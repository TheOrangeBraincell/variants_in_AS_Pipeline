# -*- coding: utf-8 -*-
"""
date: 29.03.22
name: splicesites_in_bam.py
author: Mirjam Karlsson-Müller


"""


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

parser.add_argument('--database', '-db',
                    help='BED file containing already known/annotated splice junctions from the ucsc.')
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
                                exon_count="1,1"
                                block_sizes="tbd"
                                #Flag for being in database, default =False
                                in_db=False
                                name=str(exon1_end)+"_"+str(exon2_start)
                                #The plus for strand is given rn due to it running only for ESR1. Adjust later with "reads", l68
                                splicesite[name]=[chrom, exon1_end-1, exon2_start, "+",
                                                  exon_count, block_sizes, "0,"+str(exon2_start-exon1_end), in_db]
                                
                        #remove used part of the string
                        current_cigar=re.sub(r'^.*?N','N', 
                                             current_cigar).lstrip("N")
                        current_length+=exon2_start
#%%
"2. Read in known splice junctions for corresponding region."

db_junctions=dict()
with open(args.database, 'r') as db:
    for line in db:
        #To exclude potential title lines/empty lines, formatting mistakes
        if bool(re.search(r"([a-z]{3}[X,M,Y]?\d*)(?:.+[random]{6})?\t(\d+)\t(\d+)\t([A-Z]+\d+\.*\d*)\t0\t(\-?\+?)\t\d+\t\d+\t0\t(\d+)\t([\d,]+)\t([\d,]+)", line))==True:
            entry=re.search(r"([a-z]{3}[X,M,Y]?\d*)(?:.+[random]{6})?\t(\d+)\t(\d+)\t([A-Z]+\d+\.*\d*)\t0\t(\-?\+?)\t\d+\t\d+\t0\t(\d+)\t([\d,]+)\t([\d,]+)", line)
            chrom=entry.group(1)
            start=entry.group(2)
            stop=entry.group(3)
            gene_name=entry.group(4)
            strand=entry.group(5)
            number_exons=entry.group(6)
            exon_size=entry.group(7)
            rel_start=entry.group(8)
            name=str(int(start)+1)+"_"+stop
            db_junctions[name]=[chrom, start, stop, gene_name, strand, number_exons, exon_size, rel_start]


  
#%%

"3. Write the results in BED file"
with open(args.out, 'w') as out:
    out.write("\n\n")
    #iterate through junctions in ucsc file
    for name in db_junctions:
        if name in splicesite:
            RGB=0,0,255
            splicesite[name][7]=True
        else:
            RGB=0,255,0
        #estrogen is on + strand, so we exclude - for now.
        if db_junctions[name][4]=="+" and db_junctions[name][0]=="chr6" and int(db_junctions[name][1])>151625000 and int(db_junctions[name][2])<152103274:
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    db_junctions[name][0], db_junctions[name][1], db_junctions[name][2],
                    db_junctions[name][3], 0, db_junctions[name][4], db_junctions[name][1],
                    db_junctions[name][2], RGB, db_junctions[name][5], db_junctions[name][6],
                    db_junctions[name][7]))
    #add new ones.
    for name in splicesite:
        if splicesite[name][7]==False:
            RGB=255,0,0
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                splicesite[name][0], splicesite[name][1], splicesite[name][2], name,
                0, splicesite[name][3], splicesite[name][1], splicesite[name][2], RGB,
                splicesite[name][4], splicesite[name][5], splicesite[name][6]))    

        
                

                                
