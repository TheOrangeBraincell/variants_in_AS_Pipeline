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

go more generous on the coordinates.

python events_gencode.py -b alignment.bam -gc geneID_hg38_GENCODE39.tsv -db hg38_GENCODE39.bed -o test_out.txt -c "chr6:151690000-15210000"

but if we go too generous on the coordinates, we might find things on the wrong strand.
Though this is only based on gencode, and it would tell us if a transcript would 
belong to a different gene. So for now it isnt a problem.  
But it might become one when we include novel junctions. keep that in mind.  

"""

# %% Imports
import argparse
import re
import pysam
import glob

# %% 0. argparse

"1. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Find Events',
                                 usage='%(prog)s -b BAM-INPUT -o OUTPUT ',
                                 description="""Returns alternative potential
                                 Casette Exons and their PSI scores, 
                                 based on sample bam files.""")

# Bam file should be later replaced with bam folder. (to allow several samples)
parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing among \
                        others, the bam files.')
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

#print("1 done")
#print(gene_dict)
#here we still have all trans_ids per chromosome.



# %% 2. Extract exons from gencode bed file
"""  """
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

#print("2 done")
#for key in transID_exons:
#    print(transID_exons[key])

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

#print("3 done")


# %% 4. Count reads in BAM file for every CE and calculate PSI.

# This command would count all reads per exon, we only want spliced.
# exons[3]+=pysam.AlignmentFile.count(samfile, contig=exons[0], start=exons[1],
#                      stop=exons[2])




"""Run for bam file of each sample"""

table=dict()   
#Initiate counter to keep track of progress.
#counter=0
for gene_id in gene_exons:
    #prints updates on progress of script. 
    #print("sample:"+bam_file_list.index(file)+"/"+len(bam_file_list)+" ,gene id", counter, "/", len(gene_exons), end="\r")
    #counter+=1
    
    #for some genes, there might only be one or two annotated exons. Thus no
    #CE. These we skip to avoid indexing problems.
    if len(gene_exons[gene_id])<3:
        continue
    #Extract strand and chromosome of gene.
    gene_strand=gene_exons[gene_id][0][3]
    chrom=gene_exons[gene_id][0][0]
        
    """Here on out, refering to smaller and bigger coordinate, as for plus 
    strand the smaller coordinate is start and the bigger stop, but for minus 
    strand its the other way around."""
    #sort after smaller coordinate, specify "first" exon and remove it.
    gene_exons[gene_id].sort(key=lambda x: x[1])
    small_first_exon=gene_exons[gene_id][0][1] #smaller coord.
    gene_exons[gene_id]=gene_exons[gene_id][1::]
    #sort after bigger positions, specify last exon and remove it
    gene_exons[gene_id].sort(key=lambda x: x[2])
    big_last_exon=gene_exons[gene_id][-1][2]
    gene_exons[gene_id] = gene_exons[gene_id][0:-1]
      
    #iterate through exons of each gene_ID 
    for exons in gene_exons[gene_id]:
        #initiate count values for PSI scores
        junction3=0
        junction5=0
        splice_junction=0
            
        #reset read dictionary (used to remove mate reads at same junction)
        read_dict=dict()
        
        #Find all the bam files in the data folder.
        argument_glob=args.samples+"/**/*.bam"
        bam_file_list=glob.glob(argument_glob, recursive=True)
        #Iterate through all the alignment files.
        previous_sample_name=""
        for file in bam_file_list:
            print(file)
            samfile=pysam.AlignmentFile(file, 'rb', index_filename=file[0:-1]+"i")
            sample_name=file.split("/")[1]
            #Open bam file at the coordinates of this gene. returns iterable.
            samfile_fetch=samfile.fetch(chrom, small_first_exon, big_last_exon)
            #iterate through reads at gene's location:
            for read in samfile_fetch:
                #only use junction reads!
                if re.search(r'\d+M\d+N\d+M', read.cigarstring):
                    name=read.query_name
                    start=int(read.reference_start)
                    chrom=read.reference_name
                    read_length=read.infer_query_length()
                    
                    "Filtering pt 1"
                    # Exclude non-primary alignments
                    if read.is_secondary:
                        continue
                    #exclude second read of pair, if maps to overlapping region.
                    if name in read_dict:
                        if read_dict[name][0]<=start<=read_dict[name][1] or \
                        read_dict[name][0]<=start+read_length<=read_dict[name][1]:
                            continue
                    else:
                        read_dict[name]=[start, start+read_length]
                    
                    "get strand information"
                    if read.mate_is_reverse and read.is_read1:
                        strand="-"
                    elif read.mate_is_reverse and read.is_read2:
                        strand="+"
                    elif read.mate_is_forward and read.is_read1:
                        strand="+"
                    elif read.mate_is_forward and read.is_read2:
                        strand="-"
                        
                    "Filtering pt 2"
                    #Exclude reads on the wrong strand
                    if strand!= exons[3]:
                        continue
                    
                    #To allow for several junctions in one cigar string, 
                    #we require a loop that keeps looking for a pattern.
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
                            
                        
                        "Filtering pt 3"
                        # If junction has less than 3 bp on either side of the 
                        # intron, skip it:
                        if exon1 < 3 or exon2 < 3:
                            # update cigar string
                            current_cigar = re.sub(r'^.*?N', 'N',
                                                   current_cigar).lstrip("N")
                            current_start= exon2_start
                            continue
                        
                        "Counts:"
                        if smaller_ex[1]<exons[1] and bigger_ex[0]>exons[2]:
                            splice_junction+=1
                        elif smaller_ex[1]<exons[1] and \
                            exons[1]<=bigger_ex[1] <= exons[2]:
                            junction5+=1
                        elif bigger_ex[0]>exons[2] and \
                            exons[1]<=smaller_ex[0]<=exons[2]:
                            junction3+=1
                        
                        # update cigar string
                        current_cigar = re.sub(r'^.*?N', 'N', 
                                               current_cigar).lstrip("N")
                        current_start= exon2_start
            
            "Calculate PSI for exon, write it into file"
            #those that did not appear in reads, will be 0. They get a NAN value.
            #Also filter out low counts, as they are not significant %.
            if int((junction5+junction3+splice_junction))<=10:
                PSI="NAN"
            else:
                PSI=(junction5+junction3)/(junction5+junction3+splice_junction) 
            
            if bam_file_list.index(file)==0:
                table[gene_id]={exons[0]+"_"+str(exons[1])+"_"+str(exons[2]):[PSI]}
            else:
                table[gene_id][exons[0]+"_"+str(exons[1])+"_"+str(exons[2])].append(PSI)


#%% 5. Write into output file.
            
"open output file"
out=open(args.out, "w")
title="Location\t"
for i in range(1,len(bam_file_list)):
    title+="s"+str(i)+"\t"
out.write(title+"\n")

for gene_id in table:
    #Title for new gene, so they are clearly seperated in outputfile
    out.write("# gene id: {}, strand: {}\n".format(gene_id, gene_strand))
    for position in table[gene_id]:
        out.write("{}\t{}\n".format(position, "\t".join(table[gene_id][position])))        
out.close       
        

print("\nDone!")

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        