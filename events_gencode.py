# -*- coding: utf-8 -*-
"""
date: 21-04-22
title: events_gencode.py
author: Mirjam Karlsson-Müller
    
Abbreviations:
CE= casette exon
AA= alternative acceptor
AD= alternative donor

Finds alternative splicing events in gencode bed file, and calculates the PSI
for these events, based on the bam files from the SCAN-B database.

Currently only returns event table with CE events, no AA or AD or intron
retention.

python events_gencode.py -b alignment.bam -gc geneID_hg38_GENCODE39.tsv -db hg38_GENCODE39.bed -o test_out.txt -c "chr6:151690496-152103274"


"""

# %% Imports
import argparse
import re
import pysam

# %% 0. argparse

"1. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Find Events',
                                 usage='%(prog)s -b BAM-INPUT -o OUTPUT ',
                                 description="""Returns alt splice events
                                 and their PSI scores, based on samples.""")

# Bam file should be later replaced with bam folder. (to allow several samples)
parser.add_argument('--bamfile', '-b', required=True,
                    help='BAM file, containing alignment information.')
parser.add_argument('--gencode', '-gc', required=True,
                    help="""tsv file containing gene_ID and transcript IDs.""")
parser.add_argument('--database', '-db', required=True,
                    help="""BED file containing already known/annotated splice
                    junctions from the ucsc.""")
parser.add_argument('--out', '-o', required=True,
                    help="""Output bed file, where comparison data should be 
                    printed to.""")
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")

args = parser.parse_args()

# Extract input coordinates:
if args.coordinates:
    if bool(re.search(r'[a-z]{3}[MXY]?\d*:\d+-\d+', args.coordinates)) == False:
        raise argparse.ArgumentTypeError("""The given coordinates are in the 
                                         wrong format. Input as 
                                         chrX:XXXX-XXXX.""")
        quit()
    else:
        coord = re.search(r'([a-z]{3}[MXY]?\d*):(\d+)-(\d+)', args.coordinates)
        coord_chrom = coord.group(1)
        coord_start = int(coord.group(2))
        coord_stop = int(coord.group(3))

# %% 1. Dictionary linking gene_ID with transcript_IDs using gencode gtf.

"1. Make genename_transcript id dictionary based on gencode table."

gene_dict = dict()
with open(args.gencode, "r") as gene:
    for line in gene:
        if bool(re.search("(.+)\t([a-z]{3}\d*[X,M,Y]*)\t(.{1})\t(.+)", line)) == True:
            entry = re.search(
                "(.+)\t([a-z]{3}\d*[X,M,Y]*)\t(.{1})\t(.+)", line)
            transcript_id = entry.group(1)
            gene_id = entry.group(4)

            # If only looking for a region, only look at corresponding chrom
            if args.coordinates:
                if entry.group(2) != args.coordinates.split(":")[0]:
                    continue

                else:
                    if gene_id in gene_dict:
                        gene_dict[gene_id].append(transcript_id)
                    else:
                        gene_dict[gene_id] = [transcript_id]

            else:
                if gene_id in gene_dict:
                    gene_dict[gene_id].append(transcript_id)
                else:
                    gene_dict[gene_id] = [transcript_id]

print("1 done")

# %% 2. Extract exons from gencode bed file
"""  """
transID_exons = dict()

with open(args.database, "r") as database:
    for line in database:
        # To exclude potential title lines/empty lines, formatting mistakes
        # Only takes chr[] and chr[]_random lines, in accordance with bam.
        if bool(re.search(r"([a-z]{3}[X,M,Y]?\d*)\t(\d+)\t(\d+)\t([A-Z]+\d+"
                          r".*\d*)\t0\t(\-?\+?)\t\d+\t\d+\t0\t(\d+)\t([\d,]+)"
                          r"\t([\d,]+)", line)) == True or bool(re.search(
                              r"([a-z]{3}[X,M,Y]?\d*).+_random\t(\d+)\t(\d+)\t"
                              r"([A-Z]+\d+\.*\d*)\t0\t(\-?\+?)\t\d+\t\d+\t0\t"
                              r"(\d+)\t([\d,]+)\t([\d,]+)", line)) == True:
            # specify what the groups in the line correspond to.
            entry = re.search(r"([a-z]{3}[X,M,Y]?\d*).*\t(\d+)"
                              r"\t(\d+)\t([A-Z]+\d+\.*\d*)\t0\t(\-?\+?)\t\d+\t"
                              r"\d+\t0\t(\d+)\t([\d,]+)\t([\d,]+)", line)
            chrom = entry.group(1)
            start = int(entry.group(2))
            stop = entry.group(3)
            trans_name = entry.group(4)
            strand = entry.group(5)
            number_exons = entry.group(6)
            exon_size = entry.group(7)
            rel_start = entry.group(8)
            
            #If coordinates are given, we only want entries within.
            if args.coordinates:
                if int(start)<coord_start or int(stop) >coord_stop or chrom != coord_chrom:
                    continue
            # We are interested in the single exons of each entry.
            size_list = exon_size.split(",")
            start_list = rel_start.split(",")
            current_start = start

            for i in range(0, int(number_exons)-1):
                exon_start = current_start+1
                exon_end = current_start+int(size_list[i])-1

                if trans_name in transID_exons:
                    transID_exons[trans_name].append([chrom,exon_start, exon_end,0])
                else:
                    transID_exons[trans_name] = [[chrom,exon_start, exon_end,0]]

print("2 done")

# %% 3. Make gene/exon dictionary.
"""Every exon that is not the first or last exon in a gene, is treated as 
a potential CE"""

gene_exons = dict()

# Make a dictionary containing all exons belonging to a gene id.
for gene_id in gene_dict:
    for trans_id in gene_dict[gene_id]:
        if gene_id in gene_exons:
            if trans_id in transID_exons:
                for exon in transID_exons[trans_id]:
                    gene_exons[gene_id].append(exon)
        else:
            if trans_id in transID_exons:
                for exon in transID_exons[trans_id]:
                    if transID_exons[trans_id].index(exon) == 0:
                        gene_exons[gene_id] = [exon]
                    else:
                        gene_exons[gene_id].append(exon)

"""
# Sort the exons belonging to one gene id, to find first and last.
for gene_id in gene_exons:
    #sort after start positions of exons. remove first exon.
    gene_exons[gene_id].sort(key=lambda x: x[1])
    
    gene_exons[gene_id]=gene_exons[gene_id][1::]
    #sort after end positions. remove last exon.
    gene_exons[gene_id] = gene_exons[gene_id][0:-1]
    """
print("3 done")

# %% 4. Read BAM file
"""
exon_reads = []

# Splicejunction has alignment pattern:
pattern_junction = re.compile(r'\d+M\d+N\d+M')

if args.coordinates:
    command = ['samtools', 'view', '-h', '-F', '0x100', args.bamfile,
               args.coordinates]
else:
    command = ['samtools', 'view', '-h', '-F', '0x100', args.bamfile]


    for line in sam:
        if line.startswith("@"):
            continue
        # if the line follows the pattern containing one or more splice junctions:
        if bool(re.search(r"(\d+)\t([a-z]{3}\d+)\t(\d+)\t(\d+)\t"
                          r"(\d+M\d+N\d+M)*\t.+YT:Z:([A-Z]{2})",
                          line)) == True:
            entries = re.search(r"(\d+)\t([a-z]{3}\d+)\t(\d+)\t(\d+)\t"
                                r"(\d+M\d+N\d+M)\t.+YT:Z:([A-Z]{2})", line)
            # Assign variables to entries in line.
            flag = entries.group(1)
            chrom = entries.group(2)
            start = entries.group(3)
            cigar = entries.group(5)
            pair = entries.group(6)
    
            # We only include concordant pairs:
            if pair != "CP":
                continue
    
            #To allow for several junctions in one cigar string, we require a loop that keeps 
            #looking for a pattern.
            current_cigar = cigar
            current_length = 0
            while bool(re.search(pattern_junction, current_cigar)) == True:
                print(current_cigar)
                junction = re.search(r'(\d+)M(\d+)N(\d+)M', current_cigar)
                exon1 = int(junction.group(1))
                intron = int(junction.group(2))
                exon2 = int(junction.group(3))
                exon1_start = int(start)+current_length-1
                exon1_end = exon1_start+exon1
                exon2_start = exon1_end+intron + 1
                exon2_end = exon2_start+exon2
    
               # If junction has less than 3 bp on either side of the intron, 
               # remove it:
                if exon1 < 3 or exon2 < 3:
                    # update cigar string
                    current_cigar = re.sub(r'^.*?N', 'N',
                                           current_cigar).lstrip("N")
                    current_length += exon2_start
                    continue
    
                # if there's input coordinates, exons need to be in the region.
                if exon1_end < coord_start or exon2_start > coord_stop:
                    # update cigar string
                    current_cigar = re.sub(r'^.*?N', 'N',
                                           current_cigar).lstrip("N")
                    current_length += exon2_start
                    continue
    
                # Exclude non-primary alignments (flag 256)
                if len(bin(int(flag))) >= 11 and str(bin(int(flag)))[-9] == "1":
                    # update cigar string
                    current_cigar = re.sub(r'^.*?N', 'N',
                                           current_cigar).lstrip("N")
                    current_length += exon2_start
                    continue
    
                exon_reads.append([exon1_start, exon1_end],
                                  [exon2_start, exon2_end])
    
                # update cigar string
                current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                current_length += exon2_start

# count reads per exon.
total_reads = len(exon_reads)

print("before loop fine.")

for gene in gene_exons:
    for exons in gene_exons[gene]:
        #iterate through reads:
        for reads in exon_reads:
            if reads[0]>=exons[0] and reads[1]<=exons[1]:
                exons[2]+=1
"""

# This command would count all reads per exon, we only want spliced.
# exons[3]+=pysam.AlignmentFile.count(samfile, contig=exons[0], start=exons[1],
#                      stop=exons[2])
# "reads" not defined. 
# exons[3]+=pysam.AlignmentFile.count(samfile, contig=exons[0], 
#                                     start=exons[1], stop=exons[2], 
#                                     read_callback=check_read(read))

def check_read(read):
    if re.search(r'\d+M\d+N\d+M', read.split("\t")[4]):
        return True

"open output file"
out=open(args.out, "w")
"Parse through file manually."
samfile=pysam.AlignmentFile(args.bamfile, 'rb')

"""Trying out stuff with pysam
for read in samfile:
    #print(read.to_string()) #THATS THE LINE AS I KNOW IT. HORRAY.
    cigar=read.cigarstring
    #flag=read.flag, not needed.
    chrom=read.reference_id #this ones still fucked.
    #print(cigar)
    #exclude secondary alignments
    if read.is_secondary:
        continue
    #filter out not proper pairs. 
    
    #read1: read.is_read1, read2: read.is_read2
    "read off strand"
    if read.mate_is_reverse and read.is_read1:
        strand="-"
    elif read.mate_is_reverse and read.is_read2:
        strand="+"
    elif read.mate_is_forward and read.is_read1:
        strand="+"
    elif read.mate_is_forward and read.is_read2:
        strand="-"
    
    
"""
#pysam gives the line a new structure. unlike samtools. 
#chromosome entry complete nonsense. way too many different options.
for gene_id in gene_exons:
    #for some genes, there might only be one or two annotated exons. Thus no
    #CE. These we skip to avoid indexing problems.
    if len(gene_exons[gene_id])<3:
        continue
    #sort after start positions, specify start exon and remove it.
    gene_exons[gene_id].sort(key=lambda x: x[1])
    start_first_exon=gene_exons[gene_id][0][1]
    gene_exons[gene_id]=gene_exons[gene_id][1::]
    #sort after stop positions, specify stop exon and remove it
    gene_exons[gene_id].sort(key=lambda x: x[2])
    stop_last_exon=gene_exons[gene_id][-1][2]
    gene_exons[gene_id] = gene_exons[gene_id][0:-1]
    chrom=gene_exons[gene_id][0][0]
    
    
    
    #Open bam file at the coordinates of this gene. returns iterable.
    samfile_fetch=samfile.fetch(chrom, start_first_exon, stop_last_exon)
    #iterate through exons of gene
    for exons in gene_exons[gene_id]:
        print(exons)
        #initiate count values:
        junction3=0
        junction5=0
        splice_junction=0
        #iterate through reads at gene's location:
        for read in samfile_fetch:
            #turn into string object.
            read_str=read.to_string()
            if re.search(r'\d+M\d+N\d+M', read.cigarstring):
                cigar=read.cigarstring
                start=read.reference_start
                chrom=read_str.split("\t")[2]
                
                "Filtering pt 1"
                # Exclude non-primary alignments
                if read.is_secondary:
                    continue
                
                "get strand information"
                if read.mate_is_reverse and read.is_read1:
                    strand="-"
                elif read.mate_is_reverse and read.is_read2:
                    strand="+"
                elif read.mate_is_forward and read.is_read1:
                    strand="+"
                elif read.mate_is_forward and read.is_read2:
                    strand="-"
                
                #To allow for several junctions in one cigar string, we require a loop that keeps 
                #looking for a pattern.
                current_cigar = cigar
                current_length = 0
                while re.search(r'\d+M\d+N\d+M', current_cigar):
                    #assign splice junction variables
                    junction = re.search(r'(\d+)M(\d+)N(\d+)M', current_cigar)
                    exon1 = int(junction.group(1))
                    intron = int(junction.group(2))
                    exon2 = int(junction.group(3))
                    exon1_start = int(start)+current_length-1
                    exon1_end = exon1_start+exon1
                    exon2_start = exon1_end+intron + 1
                    exon2_end = exon2_start+exon2
                    
                    "Filtering pt 2"
                    # If junction has less than 3 bp on either side of the intron, 
                    # skip it:
                    if exon1 < 3 or exon2 < 3:
                        # update cigar string
                        current_cigar = re.sub(r'^.*?N', 'N',
                                               current_cigar).lstrip("N")
                        current_length += exon2_start
                        continue
                    
                    "Count: this is for forward string"
                    #print(exon1_end, exons[1], exon2_start, exons[2])
                    if exon1_end < exons[1] and exon2_start>exons[2]:
                        print("case1")
                        splice_junction+=1
                    elif exon1_end<exons[1] and exons[1]<exon2_end<exons[2]:
                        print("case2")
                        junction5+=1
                    elif exon2_end>exons[2] and exons[1]<exon1_start<exons[2]:
                        print("case3")
                        junction3+=1
                    
                    # update cigar string
                    current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                    current_length += exon2_start
        
        "Calculate PSI for exon, write it into file"
        #those that did not appear in reads, will be 0. They get a NAN value.
        if (junction5+junction3+splice_junction)==0:
            PSI="NAN"
        else:
            PSI=(junction5+junction3)/(junction5+junction3+splice_junction)
        out.write("{}\t{}\n".format(str(exons[1])+"_"+str(exons[2]), PSI))
        
out.close       
        

print("4. done")

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        