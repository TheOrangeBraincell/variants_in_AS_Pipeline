# -*- coding: utf-8 -*-
"""
date: 10-05-22
title: events_gencode.py
author: Mirjam Karlsson-Müller
    
Abbreviations:
    CE= casette exon
    AA= alternative acceptor
    AD= alternative donor

Description:
    Finds potential CE in gencode bed file, and calculates the PSI
    for these events, based on the bam files from the SCAN-B database.


Instructions:
    Run in command line. Requires database files.
    The tsv file contains information on what transcription ID's belong to what
    gene ID and the bed file contains annotated exons per transcription ID.
    Also requires a folder with alignment files (bam) and corresponding index
    files (bai) with same name, but different ending.
    
    #with coordinates for Estrogen Receptor
    python events_gencode.py -s . -gc geneID_hg38_GENCODE39.tsv -db hg38_GENCODE39.bed -o CE_ESR1.txt -c "chr6:151600000-152150000" 

    #no coordinates
    python events_gencode.py -s . -gc geneID_hg38_GENCODE39.tsv -db hg38_GENCODE39.bed -o CE_PSI_all_out.txt  

Possible Bugs:
    - If there is no .bai file for a .bam file, then the bam file cannot be
    opened and thus not processed. The same happens if the .bai file has a 
    different name from the bam file or is in a different folder.
    

"""

# %% Imports
import argparse
import re
import pysam
import glob
import time

#%% Time

start_time=time.time()

# %% 0. argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Find GENECODE events',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     [-c] "chrX:XXXXXX-XXXXXX" -db BEDFILE \
                                         -gc GENCODE_TSV',
                                 description="""Returns potential
                                 Casette Exons and their PSI scores, 
                                 based on sample bam files.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing among \
                        others, the bam files.')
parser.add_argument('--gencode', '-gc', required=True,
                    help="""tsv file containing gene_ID and transcript IDs.""")
parser.add_argument('--database', '-db', required=True,
                    help="""BED file containing already known/annotated splice
                    junctions from the ucsc.""")
parser.add_argument('--out', '-o', required=True,
                    help="""Output txt file, where PSI scores are printed to.
                    """)
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


# %% 2. Extract exons from gencode bed file
transID_exons = dict()

with open(args.database, "r") as database:
    for line in database:
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
                                                      exon_cor[1], strand])
                else:
                    transID_exons[trans_name] = [[chrom, exon_cor[0], 
                                                  exon_cor[1], strand]]



#Remove first and last exon for each transcript ID, as they can not be CE.

filtered = {k: v for k, v in transID_exons.items() if len(v) > 2}
transID_exons=filtered

for trans_name in transID_exons:
    #sort after smaller coordinate, specify "first" exon and remove it.
    transID_exons[trans_name].sort(key=lambda x: x[1])
    small_first_exon=transID_exons[trans_name][0][1] #smaller coord.
    transID_exons[trans_name]=transID_exons[trans_name][1::]
    #sort after bigger positions, specify last exon and remove it
    transID_exons[trans_name].sort(key=lambda x: x[2])
    big_last_exon=transID_exons[trans_name][-1][2]
    transID_exons[trans_name] = transID_exons[trans_name][0:-1]

#print("2 done")
#for key in transID_exons:
#    print(transID_exons[key])

# %% 3. Make gene/exon dictionary.
"""Every exon that is left in the transID dictionary is treated as 
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
                    #Avoid duplicate exons. 
                    if exon in transID_exons[trans_id]:
                        continue
                    else:
                        gene_exons[gene_id].append(exon)

print("Creating Database Dictionary: Done!")

# %% 4. Extract reads from bam files.

"""Iterating through all BAM files, extracting the splicejunction reads from 
and saving them by sample of origin. """
#Find all the bam files in the data folder.
argument_glob=args.samples+"/**/*.bam"
bam_file_list=glob.glob(argument_glob, recursive=True)

#If no bam files are found, quit the program.
if len(bam_file_list)==0:
    print("""There were no bam files found in the input folder. Please make 
          sure to use the right input folder. The bam files can be in any 
          subfolder of the input folder.""")
    quit()



"Initializing variables"
#Progress tracker
current_file=0
total_files=len(bam_file_list)
percentage=round(100*current_file/total_files,2)
print("Reading alignment files: ", "{:.2f}".format(percentage), "%", end="\r")
#to keep track of alignment files
previous_sample=""
sample_names=[]
#Initiate read dictionary
read_dict=dict()
#Initiate sample dictionary
sample_dict=dict()
#start loop
for file in bam_file_list:
    sample_name=file.split("/")[2]
    #If there is no entry for this sample in the dictionary, initiate new entry
    if sample_name!=previous_sample:
        sample_dict[sample_name]=[]
    #Index file has the same name, but bai ending instead of bam.
    index_file=file[0:-1]+"i"
    #open file
    samfile=pysam.AlignmentFile(file, 'rb', index_filename=index_file)
    total_reads=samfile.count()
    #make it iterable.
    #if coordinates are given, only fetch that part of file
    if args.coordinates:
        samfile=samfile.fetch(coord_chrom, coord_start, coord_stop)
    else:
        samfile=samfile.fetch()
    #Progress tracker
    current_read=0
    for read in samfile:
        #print("\nHERE\n")
        "Filtering the reads:"
        #if read is not from a splice junction
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        # Exclude non-primary alignments
        if read.is_secondary:
            continue
        #exclude second read of pair, if maps to overlapping region.
        name=read.query_name
        start=int(read.reference_start)
        chrom=read.reference_name
        read_length=int(read.infer_query_length())
        if name in read_dict:
            if read_dict[name][0]<=start<=read_dict[name][1] or \
            read_dict[name][0]<=start+read_length<=read_dict[name][1]:
                continue
        else:
            read_dict[name]=[start, start+read_length]
        
        "Get strand information"
        if read.mate_is_reverse and read.is_read1:
            strand="-"
        elif read.mate_is_reverse and read.is_read2:
            strand="+"
        elif read.mate_is_forward and read.is_read1:
            strand="+"
        elif read.mate_is_forward and read.is_read2:
            strand="-"
        
        "to allow for several splice junctions in one read"
        current_cigar = read.cigarstring
        current_start = int(start)
        while re.search(r'\d+M\d+N\d+M', current_cigar):
            #assign splice junction variables
            junction = re.search(r'(\d+)M(\d+)N(\d+)M', current_cigar)
            if strand=="+":
                exon1 = int(junction.group(1))
                intron = int(junction.group(2))
                exon2 = int(junction.group(3))
                exon1_start = current_start
                exon1_end = exon1_start+exon1+1
                exon2_start = exon1_end+intron
                exon2_end = exon2_start+exon2 +1
                smaller_ex=[exon1_start, exon1_end]
                bigger_ex=[exon2_start, exon2_end]
                
            if strand =="-":
                exon2 = int(junction.group(1))
                intron = int(junction.group(2))
                exon1 = int(junction.group(3))
                exon2_end=current_start-1
                exon2_start=exon2_end+exon2
                exon1_end=exon2_start+intron
                exon1_start=exon1_end+exon1
                smaller_ex=[exon2_end, exon2_start]
                bigger_ex=[exon1_end, exon1_start]
            
            
            "Filtering short reads"
            # If junction has less than 3 bp on either side of the 
            # intron, skip it:
            if exon1 > 3 and exon2 > 3:
                sample_dict[sample_name].append([chrom, smaller_ex, bigger_ex, 
                                                 strand])
            
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= bigger_ex[0]
        
        #Progress Update
        current_read+=1
        percentage_reads=percentage+round(100*current_read/(total_reads*3),2)
        print("Reading alignment files: ", "{:.2f}".format(percentage_reads), "%", end="\r")
        
    sample_names.append(sample_name)
    previous_sample=sample_name
    current_file+=1
    percentage=round(100*current_file/total_files,2)
    print("Reading alignment files: ", "{:.2f}".format(percentage), "%", end="\r")
    #Sort entries for each sample
    sample_dict[sample_name]=sorted(sample_dict[sample_name])

print("Reading alignment files: Done!         \n", end="\r")

#Ordering so the data is in the sample order
sample_names=sorted(list(set(sample_names)))
#%% 5. Counting reads for each exon in each gene, for each sample.

   
"open output file"
out=open(args.out, "w")
title="Location\t"+"\t".join(sample_names)

    
out.write(title+"\n")

#Progress tracker
total_iterations= len(gene_exons)*len(sample_dict)*len(list(gene_exons.values()))
current_iteration=0
percentage=round(current_iteration/total_iterations,2)
#print("Counting Reads: ",percentage, "%", end="\r")

for gene_id in gene_exons:
    #extract strand and chrom from first exon
    strand=gene_exons[gene_id][0][3]
    
    #remove duplicate exons.
    gene_exons[gene_id]=sorted(list(set(map(tuple,gene_exons[gene_id]))))
    #print(gene_exons[gene_id])
    
    #Line for gene id, so that entries in table are seperated.
    out.write("#"+gene_id+", "+strand+"\n")
    for exon in gene_exons[gene_id]:
        PSI_scores=[]
        chrom=exon[0]
        start=exon[1] #Note that these are start and stop coor
        stop=exon[2] #not start and stop of exon. Strand dependent.
        strand=exon[3]
                
        for sample in sample_names:
            current_iteration+=1
            percentage=round(current_iteration/total_iterations,2)
            #print("Counting Reads: ",percentage, "%", end="\r")
            #initiate count values for PSI scores
            junction3=0
            junction5=0
            splice_junction=0
            
            for read in sample_dict[sample]:
                smaller_exon=read[1]
                bigger_exon=read[2]
                
                #Exclude read if it is on the wrong strand.
                if read[3]!=strand:
                    continue
                
                "Counts:"
                if smaller_exon[1]<start and bigger_exon[0]>stop:
                    splice_junction+=1
                elif smaller_exon[1]<start and \
                    start<=bigger_exon[1] <= stop:
                    junction5+=1
                elif bigger_exon[0]>stop and \
                    start<=smaller_exon[0]<=stop:
                    junction3+=1
            
            
            "Calculate PSI for exon, write it into file"
            #those that did not appear in reads, will be 0. They get a NAN value.
            #Also filter out low counts, as they are not significant %.
            if int((junction5+junction3+splice_junction))<=10:
                PSI="NAN"
            else:
                PSI=round((junction5+junction3)/(junction5+junction3+splice_junction),3) 
               
            PSI_scores.append(str(PSI))    
            #print(PSI, sample)   
        out.write("{}\t{}\n".format(exon[0]+"_"+str(exon[1])+"_"+str(exon[2]), 
                             "\t".join(PSI_scores)))


out.close

print("Counting Reads: Done!         \n", end="\r")


#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
