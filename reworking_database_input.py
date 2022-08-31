# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:15:01 2022


python reworking_database_input.py -s ../Sample_Data/ -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -o CE_ESR1.tsv -c "chr6:151600000-152150000"

@author: mirjam
"""
#%% Imports
import argparse
import re
import pysam
import glob
import time

#%% Time

start_time=time.time()


#%% argparse

parser = argparse.ArgumentParser(prog='Find Alternative Splicing events',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     -g GENCODE-TSV -n REFSEQ-TSV \
                                     [-c] "chrX:XXXXXX-XXXXXX" \
                                         [-n] "ABCDE" ', 
                                 description="""Returns potential
                                 Casette Exons based on GENCODE39 and RefSeq
                                 and their PSI scores, based on sample 
                                 bam files.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing among \
                        others, the bam files.')
parser.add_argument('--gencode', '-g', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from GENCODE39 as well as gene names.""")
parser.add_argument('--refseq', '-r', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from RefSeq as well as gene names.""")
parser.add_argument('--out', '-o', required=True,
                    help="Output tsv file, where PSI scores are printed to.")
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")
parser.add_argument('--name', '-n', type=str,
                    help="""Symbol of gene of interest. f.e. ESR1. 
                    (Symbol by HUGO Gene Nomenclature Committee).""")

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


#%% 1. Database Input.

"""NCBI and UCSC use two different systems of Nomenclauture. To be able to link 
the gene, both sets of gene IDs need to be linked to a common gene name."""

#Create a dictionary linking GENCODE ids to Gene Names
gencode_names=dict()
with open(args.gcnames, "r") as gcnames:
    for line in gcnames:
        if line.startswith("#"):
            continue
        line=line.strip("\n")
        name=line.split("\t")[5]
        if name in gencode_names:
            gencode_names[name].append(line.split("\t")[0])
        else:
            gencode_names[line.split("\t")[5]]=[line.split("\t")[0]]
        

#Create a dictionary linking Refseq ids to gene names

refseq_names=dict()
with open(args.rsnames, "r") as rsnames:
    for line in rsnames:
        if line.startswith("#"):
            continue
        line=line.strip("\n")
        name=line.split("\t")[5]
        if name in refseq_names:
            refseq_names[name].append(line.split("\t")[0])
        else:
            refseq_names[line.split("\t")[5]]=[line.split("\t")[0]]

#make a list of all gene_names
gene_names=list(gencode_names.keys())
for keys in refseq_names:
    if keys in gene_names:
        continue
    else:
        gene_names.append(keys)


" geneID_transcript id dictionary based on gencode table."

gene_dict = dict()
with open(args.gencode, "r") as gene:
    for line in gene:
        if re.search("(.+)\t([a-z]{3}\d*[X,M,Y]*)\t(.{1})\t(.+)", line):
            entry = re.search(
                "(.+)\t([a-z]{3}\d*[X,M,Y]*)\t(.{1})\t(.+)", line)
            transcript_id = entry.group(1)
            gene_id = entry.group(4)
            
            # If only looking for a region, only look at corresponding chrom
            if args.coordinates:
                if entry.group(2) != coord_chrom:
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


"This is not necessary for Refseq as transcriptIDs contain the geneIDs"


#%% 2. Read in annotated exons per transcript ID.

transID_exons=dict()

with open(args.gcbed, "r") as gc:
    for line in gc:
        # To exclude potential title lines/empty lines, formatting mistakes
        # Only takes chr[] and chr[]_random lines, in accordance with bam.
        if re.search(r"(?:([a-z]{3}[X,M,Y]?\d*)|([a-z]{3}[X,M,Y]?\d*).+_random)"
                     r"\t(\d+)\t(\d+)\t([A-Z]+\d+.*\d*)\t0\t(\-?\+?)\t\d+\t\d+"
                     r"\t0\t(\d+)\t([\d,]+)\t([\d,]+)", line):
            # specify what the groups in the line correspond to.
            entry = re.search(r"([a-z]{3}[X,M,Y]?\d*).*\t(\d+)"
                              r"\t(\d+)\t([A-Z]+\d+\.*\d*)\t0\t(\-?\+?)\t\d+\t"
                              r"\d+\t0\t(\d+)\t([\d,]+)\t([\d,]+)", line)
            
            strand = entry.group(5)
            chrom = entry.group(1)
            trans_name = entry.group(4)
            number_exons = entry.group(6)
            exon_size = entry.group(7)
            
            """Defining "smaller" and "bigger" coordinate for each exon,
            depending on strand, dictionary gets fixed 
            [chrom, smaller, bigger, strand] instead of confusing start and
            stop"""
            
            small=int(entry.group(2))
            big=int(entry.group(3))
            """
            #Exclude entries outside of coordinates, if given.
            if args.coordinates:
                if small<coord_start or big>coord_stop or chrom!=coord_chrom:
                    continue
            """
            # We are interested in the single exons of each entry.
            rel_small=entry.group(8)
            size_list = exon_size.split(",")
            small_list = rel_small.split(",")
            
            
            for i in range(0, int(number_exons)-1):
                if strand=="+":
                    current_start=small+int(small_list[i])
                    exon_start=current_start
                    exon_end=current_start+int(size_list[i])+1
                    exon_cor=[exon_start, exon_end]
                
                if strand=="-":
                    current_stop=small+int(small_list[i]) 
                    exon_end = current_stop-1 
                    exon_start = current_stop+int(size_list[i])+1 #exon start
                    exon_cor=[exon_end, exon_start] #so that smaller, bigger
                
                #append exon to dictionary
                if trans_name in transID_exons:
                    transID_exons[trans_name].append([chrom, exon_cor[0], 
                                                      exon_cor[1], strand, "G"])
                else:
                    transID_exons[trans_name] = [[chrom, exon_cor[0], 
                                                  exon_cor[1], strand, "G"]]

with open(args.rsbed, "r") as rs:
    for line in rs:
        # To exclude potential title lines/empty lines, formatting mistakes
        # Only takes chr[] and chr[]_random lines, in accordance with bam.
        if re.search(r"(?:([a-z]{3}[X,M,Y]?\d*)|([a-z]{3}[X,M,Y]?\d*).+_random)"
                     r"\t(\d+)\t(\d+)\t([A-Z]+\d+.*\d*)\t0\t(\-?\+?)\t\d+\t\d+"
                     r"\t0\t(\d+)\t([\d,]+)\t([\d,]+)", line):
            # specify what the groups in the line correspond to.
            entry = re.search(r"([a-z]{3}[X,M,Y]?\d*).*\t(\d+)"
                              r"\t(\d+)\t([A-Z]+\d+\.*\d*)\t0\t(\-?\+?)\t\d+\t"
                              r"\d+\t0\t(\d+)\t([\d,]+)\t([\d,]+)", line)
            
            strand = entry.group(5)
            chrom = entry.group(1)
            trans_name = entry.group(4)
            number_exons = entry.group(6)
            exon_size = entry.group(7)
            
            """Defining "smaller" and "bigger" coordinate for each exon,
            depending on strand, dictionary gets fixed 
            [chrom, smaller, bigger, strand] instead of confusing start and
            stop"""
            
            small=int(entry.group(2))
            big=int(entry.group(3))
            
            #Exclude entries outside of coordinates, if given.
            if args.coordinates:
                if small<coord_start or big>coord_stop or chrom!=coord_chrom:
                    continue
            
            # We are interested in the single exons of each entry.
            rel_small=entry.group(8)
            size_list = exon_size.split(",")
            small_list = rel_small.split(",")
            
            
            for i in range(0, int(number_exons)-1):
                if strand=="+":
                    current_start=small+int(small_list[i])
                    exon_start=current_start
                    exon_end=current_start+int(size_list[i])+1
                    exon_cor=[exon_start, exon_end]
                
                if strand=="-":
                    current_stop=small+int(small_list[i]) 
                    exon_end = current_stop-1 
                    exon_start = current_stop+int(size_list[i])+1 #exon start
                    exon_cor=[exon_end, exon_start] #so that smaller, bigger
                
                #append exon to dictionary
                if trans_name in transID_exons:
                    transID_exons[trans_name].append([chrom, exon_cor[0], 
                                                      exon_cor[1], strand, "R"])
                else:
                    transID_exons[trans_name] = [[chrom, exon_cor[0], 
                                                  exon_cor[1], strand, "R"]]


#%% 3. Assemble them into a gene_name:transID:exon nested dictionary

gene_exons=dict()

for gene in gene_names:
    gene_exons[gene]=dict()
    print(gencode_names[gene])
    "Add GENCODE ones into dictionary first"
    for gene_id in gencode_names[gene]:
        for key,value in gene_dict.items():
            print(key+"\t"+ ",".join(value)+"\n")
        quit()
        for trans_id in gene_dict[gene_id]:
            gene_exons[gene][trans_id]=transID_exons[trans_id]
    
    "Now add RefSeq into the dictionary"
    for gene_id in refseq_names[gene]:
        for trans_id in transID_exons:
            if trans_id.split(".")[0]==gene_id:
                gene_exons[gene][trans_id]=transID_exons[trans_id]

















