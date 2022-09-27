# -*- coding: utf-8 -*-
"""
date: 27-09-22
title: Variants_AS_altered.py
author: Mirjam Karlsson-Müller

Description:
    For each sample picks vcf containing the most variants and selects corresponding
    bam and gene.tsv files.
    
    The variants in the vcf files are filtered, and for those passing 
    the genotypes are extracted. If a variant does not pass filtering, it is
    saved as NEDA. For samples where the genotype at a given location is 
    unclear, gene expression data is used to predict genotype HMZR or NOEX.
    
    From GENCODE and RefSeq annotated exons are extracted and potential 
    AS events identified. For those events, PSI scores are calculated based
    on the alignment files for each sample.
    
    This one is altered to take a file with sample/folder lists.
    
Abbreviations:
    HMZA=   Genotype Homozygous Alternative Allele
    HETZ=   Heterozygous
    HMZR=   Homozygous Reference Allele
    NEDA=   Not enough data to be sure. Means there is evidence of variants, or
            alternative genotypes, but they got filtered out.
    NOEX=   no expression at this point.
    AS=     Alternative Splicing
    CE=     Casette Exon
    AA=     Alternative Acceptor
    AD=     Alternative Donor
    IR=     Intron Retention


Inputs:
    - Sample Folder
    - Database files for GENCODE and RefSeq respectively
    - output folder, where output files are stored
    - evtl. coordinates or name of a gene/region of interest
    - What type of alternative splicing event?

Instructions:
    Run in command line.
    
    #with coordinates f.e. Estrogen Receptor
    python Variants_AS.py -s ../Sample_Data/ -o ./Variants_AS_out/ -c "chr6:151690496-152103274" -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as CE
    
    #for server
    python Variants_AS.py -s /raidset/mi3258mu-se/Mirjam -o./Variants_AS_out_ESR1/ -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as CE -c "chr6:151690496-152103274"
    
    #or name f.e. Estrogen Receptor
    

    #no coordinates/name


Possible Bugs:
    - Different sample folder structure would be a problem.
    

"""


#%% Imports

import argparse
import glob
import re
import gzip
import time
import pysam
import os

#%% Time

start_time=time.time()

#%% 0. argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='VCF Parser',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     [-c] "chrX:XXXXXX-XXXXXX"',
                                 description="""Creates a genotype table out of
                                 several samples. Containing location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='file containing sample folder names containing among \
                        others, the vcf, bam and gene.tsv files.')
parser.add_argument('--out', '-o', required=True,
                    help="""Output folder, containing genotype and Psi tables.""")
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")
parser.add_argument('--gencode', '-g', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from GENCODE39 as well as gene names.""")
parser.add_argument('--refseq', '-r', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from RefSeq as well as gene names.""")
parser.add_argument('--name', '-n', type=str,
                    help="""Symbol of gene of interest. f.e. ESR1. 
                    (Symbol by HUGO Gene Nomenclature Committee).""")
parser.add_argument('--AS', '-as', type=str, required=True,
                    help="""Which type of alternative splicing event we are
                    interested in. "CE" for Casette Exons, "AA" for alternative
                    acceptors, "AD" for alternative donors, "IR" for intron
                    retention and "ALL" for all of the types. Several seperated
                    by ,.""")


args = parser.parse_args()

# Extract input coordinates, check their format.
if args.coordinates:
    if re.search(r'[a-z]{3}[MXY]?\d*:\d+-\d+', args.coordinates):
        coord = re.search(r'([a-z]{3}[MXY]?\d*):(\d+)-(\d+)', args.coordinates)
        coord_chrom = coord.group(1)
        coord_start = int(coord.group(2))
        coord_stop = int(coord.group(3))

    else:
        raise argparse.ArgumentTypeError("""The given coordinates are in the 
                                         wrong format. Input as 
                                         chrX:XXXX-XXXX.""")
        quit()


if args.AS:
    allowed_inputs=["CE", "AA", "AD", "IR", "ALL"]
    inputs=args.AS.split(",")
    for i in inputs:
        if i.upper() not in allowed_inputs:
            raise argparse.ArgumentTypeError("""The allowed abbreviations for
                                             splicing events are "CE" for 
                                             Casette Exons, "AA" for 
                                             alternative
                                             acceptors, "AD" for alternative 
                                             donors, "IR" for intron
                                             retention and "ALL" for all of 
                                             the types. Several seperated
                                             by ,.""")
            quit()

if args.name and args.coordinates:
    raise argparse.ArgumentError("""Cannot process input name and coordinates. 
                                 Please use only one or the other.""")
    quit()


#Create output directory if not there already
if not os.path.exists(args.out):
    os.makedirs(args.out)

#%% 1. Processing Database input

#outer dictionary with gene names as key
gene_dict=dict()
print("Creating Database Dictionary...", end="\r")


for file in [args.gencode, args.refseq]:
    with open(file, "r") as infile:
        for line in infile:
            # To exclude potential title lines/empty lines, formatting mistakes
            # Only takes chr[] and chr[]_random lines, in accordance with bam.
            if re.search(r"(.+)\t(?:([a-z]{3}[X,M,Y]?\d*)|([a-z]{3}[X,M,Y]?\d*)"
                         r".+_random)\t(\-?\+?)\t(\d+)\t(\d+)\t(\d+)\t([\d,]+)"
                         r"\t([\d,]+)\t(.+)", line):
                # specify groups.
                entry = re.search(r"(.+)\t([a-z]{3}[X,M,Y]?\d*).*\t"
                             r"(\-?\+?)\t(\d+)\t(\d+)\t(\d+)\t([\d,]+)"
                             r"\t([\d,]+)\t(.+)", line)
                #Assign variables to groups
                trans_ID=entry.group(1)
                chrom=entry.group(2)
                strand=entry.group(3)
                
                """To not get caught in start/stop, -/+ strand complications,
                coordinates are referred to as bigger and smaller instead."""
                
                number_exons=int(entry.group(6))
                exon_smaller=entry.group(7).split(",")[0:-1]
                exon_bigger=entry.group(8).split(",")[0:-1]
                gene_name=entry.group(9)
                if file==args.gencode:
                    db="G"
                elif file==args.refseq:
                    db="R"
                else:
                    db="You got a bug."
                
                #If choordinates are given:
                if args.coordinates:
                    #exclude entries outside of coordinates
                    if int(exon_smaller[0])<coord_start or \
                        int(exon_bigger[-1])>coord_stop or chrom!=coord_chrom:
                        continue
                
                #if name is given:
                if args.name:
                    if gene_name!=args.name:
                        continue
                
                #if its a new gene symbol, initialize inner dictionary
                if gene_name not in list(gene_dict.keys()):
                    gene_dict[gene_name]={trans_ID:[]}
                else:
                    #add transcript ID to gene_dict dictionary
                    gene_dict[gene_name][trans_ID]=[]
                
                #make entries for each exon.
                for i in range(0, number_exons-1):
                    gene_dict[gene_name][trans_ID].append([chrom, 
                                                           exon_smaller[i], 
                                                           exon_bigger[i], 
                                                           strand, db])

print("Creating Database Dictionary: Done! \n", end="\r")


#%% 2. Read off what samples are named in a file.

sample_folders=[]
with open(args.samples, "r") as samples:
    for line in samples:
        line=line.strip("\n")
        sample_folders.append(line)


vcf_file_dict=dict()
tsv_list=[]
bam_list=[]
sample_names=[]
#Find all the vcf files in the folders.
for folder in sample_folders:
    argument_glob=folder+"/**/*.vcf.gz"
    vcf_file=glob.glob(argument_glob, recursive=True)
    vcf_file_dict[vcf_file.split("/")[-5]]=vcf_file.split("/")[0]
    tsv_file="/".join(vcf_file.split("/")[0:-2])+"/t/gene.tsv"
    tsv_list.append(tsv_file)
    bam_list.append("/".join(vcf_file.split("/")[0:-2])+"/alignment.bam")
    sample_names.append(vcf_file.split("/")[-5])

"성공"

#%% 3. Read in gene_expression files.


#progress updates
total_files=len(tsv_list)
current_file=0
percentage=100*current_file/total_files
print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")
if args.name:
    gene_found=False
    starts=[]
    stops=[]

tsv_info=dict()
for file in tsv_list:
    tsv_sample_name=file.split("/")[-5]
    tsv_info[tsv_sample_name]=[]
    with open(file, "r") as tsv:
        for line in tsv:
            #if line starts with E then its a gene id
            if line.startswith("E"):
                chrom=line.split("\t")[2]
                gene_id=line.split("\t")[1]
                start=line.split("\t")[4]
                stop=line.split("\t")[5]
                fpkm=line.split("\t")[7]
                
                tsv_info[tsv_sample_name].append([chrom, start, stop, gene_id, fpkm])
                if args.name:
                    if gene_id==args.name.upper():
                        coord_chrom = chrom
                        starts.append(start)
                        stops.append(stop)
                        gene_found=True
    current_file+=1
    percentage=100*current_file/total_files
    print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

print("Reading in gene.tsv: Done!            \n", end="\r")

if args.name:
    #Catch error in name, if no match is found in tsv files.
    if gene_found==False:
        print("""The given gene symbol can not be found, please doublecheck and make
              sure there is no spaces between letters.""")
        quit()
    
    if gene_found==True:
        coord_start=min(starts)
        coord_stops=max(stops)


#%% 4. Read in variants: filter and determine genotype.

genotype_HMZR=dict()

variants=dict()
#Progress update
total_files=len(vcf_file_dict)
current_number_files=0
percentage=100*(current_number_files/total_files)
print("Filtering vcf {:.2f}%".format(percentage), end="\r")
for sample_name in sample_names:
    vcf = gzip.open(vcf_file_dict[sample_name], "rt")
    genotype_HMZR[sample_name]=[]
    #read through entries
    for line in vcf:
        #Skip headers
        if line.startswith("#"):
            continue
        chrom, position, ID, ref, alt, qual, filt, info, form, sample =line.split("\t")
        genotype=sample.split(":")[-7].strip(" ")
        
        #If coordinates or name are given, restrict area of variants.
        if args.coordinates or args.name:
            if chrom != coord_chrom:
                continue
            else:
                if int(position)< int(coord_start) or int(position)> int(coord_stop):
                    continue
        
        "Before filtering, add all variants to genotype table."
        variant_ID=chrom+"_"+position
        if variant_ID not in variants:
            variants[variant_ID]=dict()
            
        variants[variant_ID][sample_name]=[chrom, position, "ND"]
        
        "Now filter the entries."
        #Keep values with MSI<7
        if re.search(r"MSI=(\d+);", info):
            if int(re.search(r"MSI=(\d+);", info).group(1))>=7:
                continue
        else:
            continue
        
        #print("line 104")
        #Keep values with HMPOL<6
        if re.search(r"HMPOL=(\d+);", info):
            if int(re.search(r"HMPOL=(\d+);", info).group(1))>=6:
                continue
        else:
            continue
        
        #print("line 111")
        #Keep entries with GC_cont < 78%
        if re.search(r"GC_CONT=(0\.\d+);", info):
            if float(re.search(r"GC_CONT=(0\.\d+);", info).group(1))>=0.78:
                continue
        else:
            continue
        
        #print("line 118")
        #Keep variant depth >=5
        if re.search(r"VD=(\d+);", info):
            if int(re.search(r"VD=(\d+);", info).group(1))<5:
                continue
        else:
            continue
        
        #print("line 125")
        #Is a flag, so it will only be there if it applies.
        if re.search(r"low_complexity_region", info):
            continue
        
        #print("line 129")
        #6th column, filter bad quality reads.
        if qual=="." or float(qual)<55:
            continue
        #print("line 133")
        
        if re.search(r"ucsc_rep=([a-z]+);", info):
            #if re.search(r"ucsc_rep=([a-z]+);", info).group(1)=="segdup":
            continue
        
        if genotype=="0/1" or genotype=="1/0":    
            variants[variant_ID][sample_name][2]= "0/1"
        elif genotype=="1/1":
            variants[variant_ID][sample_name][2]= "1/1"
        #although it should not happen....  But it clearly is.
        elif genotype=="0/0":
            variants[variant_ID][sample_name][2]= "0/0"
            genotype_HMZR[sample_name].append(variant_ID)
        else:
            print("Invalid genotype ", genotype)
    
    "Progress updates on number of vcf files."
    current_number_files+=1
    percentage=100*(current_number_files/total_files)
    print("Filtering vcf: {:.2f}%".format(percentage),end="\r")
        
        
print("Filtering vcf: Done!            \n",end="\r")

#%% 5. Complete Genotype table

def HMZR_NOEX(variant_ID, sample):
    """
    Assigns a genotype at a specific location for a specific sample, if no
    variant is found at that location.
    
    Genotype HMZR (1/1) is given if the variant lies in a gene, which has fpkm >= 1.
    
    Otherwise genotype NE is given.

    Parameters
    ----------
    variant_ID= string
                "chrom_position"
    sample= string
            sample name f.e. S00001

    Returns
    -------
    genotype

    """
    genotype="NE"
    for gene in tsv_info[sample]:
        chrom, position= variant_ID.split("_")
        #print(gene)
        #Find gene that the variant is in
        if chrom==gene[0]:
            if int(gene[1])<int(position)<int(gene[2]):
                #Check FPKM value:
                if float(gene[4])>=1:
                    #print(gene[4])
                    genotype="0/0"
                    return genotype
        
    return genotype
            


counter=0
percentage=100*counter/len(list(variants.keys()))
print("Going through variants: {:.2f}%".format(percentage),end="\r")
#Define output file name
if args.name:
    out=args.out+"genotype_table_"+args.name+".txt"
elif args.coordinates:
    out=args.out+"genotype_table_"+args.coordinates+".txt"
else:
    out=args.out+"genotype_table.txt"
    
with open(out, "w") as out:
    out.write("Location\t"+"\t".join(sample_names)+"\n")
    for variant in variants:
        new_line=[variant]
        #iterating through sorted sample names, to always have same order.
        for sample in sample_names:
            if sample in variants[variant]:
                new_line.append(variants[variant][sample][2])
            else:
                #figure out what genotype is. HMZR or NOEX.
                genotype=HMZR_NOEX(variant, sample)
                new_line.append(genotype)
        
        out.write("\t".join(new_line)+"\n")
        counter+=1
        percentage=100*counter/len(list(variants.keys()))
        print("Going through variants: {:.2f}%".format(percentage),end="\r")

print("Going through variants: Done          \n",end="\r")


quit() #Because the rest is not ready yet, no need to run on whole genome.
#%% 6. Identify alternative splicing events in database dictionaries.

for event in inputs:
    
    "2a.) Casette Exons, Preparation"
    
    if event.upper()=="CE":
        #Remove first and last exon for each transcript, as they cannot be CE
        for gene in gene_dict:
            filtered_gene_dict={k : v for k , v in gene_dict[gene].items() 
                                if len(v) >2}
            gene_dict[gene]=filtered_gene_dict
            #print(gene_dict.keys())
            for trans_ID in gene_dict[gene]:
                #remove "first" exon, smallest start coordinate.
                gene_dict[gene][trans_ID].sort(key=lambda x: x[1])
                gene_dict[gene][trans_ID][1::]
                #remove "last" exon, biggest end coordinate.
                gene_dict[gene][trans_ID].sort(key=lambda x: x[2])
                gene_dict[gene][trans_ID][0:-1]
                
        
        """Transcript IDs are no longer required, and additionally there
        is definitely duplicate exons between R and G as well as within each
        of them. Those need to be removed."""
        
        gene_exons=dict()
        exons_db=dict()
        
        for gene in gene_dict:
            if gene_dict[gene]:
                gene_exons[gene]=[]
                for trans_ID in gene_dict[gene]:
                    for exon in gene_dict[gene][trans_ID]:
                        if exon[0:-1] not in gene_exons[gene]:
                            gene_exons[gene].append(exon[0:-1])
                        if "_".join(exon[0:-1]) in exons_db:
                            if exon[-1] in exons_db["_".join(exon[0:-1])]:
                                continue
                            else:
                                exons_db["_".join(exon[0:-1])]+=exon[-1]
                        else:
                            exons_db["_".join(exon[0:-1])]=exon[-1]


"""This is currently set up in a way that multiple AS inputs would just 
overwrite the dictionaries given by previous AS input. However, since atm
the code can only handle CE, this is a problem left to future Mirjam.

Note that rn the dictionary gene_exons is needed to output as well as 
process args.name, regardless of what AS event we are dealing with."""

#%% 7. Extract reads from bam files (SCAN-B)
"""
###############################################
#This part has not been changed since spring! #
###############################################
"""

"""Iterating through all BAM files, extracting the splicejunction reads from 
and saving them by sample of origin. """

"Initializing variables"
#Progress tracker
current_file=0
total_files=len(bam_list)
percentage=round(100*current_file/total_files,2)
print("Reading alignment files: ", "{:.2f}".format(percentage), "%", end="\r")
#Initiate read dictionary
read_dict=dict()
#Initiate sample dictionary
exons_bam=dict()
#start loop
for file in bam_list:
    sample_name=file.split("/")[-4]
    #Index file has the same name, but bai ending instead of bam.
    index_file=file[0:-1]+"i"
    #open file
    samfile=pysam.AlignmentFile(file, 'rb', index_filename=index_file)
    total_reads=samfile.count()
    #make it iterable.
    
    #if coordinates are given, only fetch that part of file
    if args.coordinates or args.name:
        samfile=samfile.fetch(coord_chrom, coord_start, coord_stop)
    else:
        samfile=samfile.fetch()
    
    exons_bam[sample_name]=[]
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
                exons_bam[sample_name].append([chrom, smaller_ex, bigger_ex, 
                                                 strand])
            
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= bigger_ex[0]
        
        #Progress Update
        current_read+=1
        percentage_reads=percentage+round(100*current_read/(total_reads*3),2)
        print("Reading alignment files: ", "{:.2f}".format(percentage_reads), "%", end="\r")
        
    current_file+=1
    percentage=round(100*current_file/total_files,2)
    print("Reading alignment files: ", "{:.2f}".format(percentage), "%", end="\r")
    #Sort entries for each sample
    exons_bam[sample_name]=sorted(exons_bam[sample_name])

print("Reading alignment files: Done!         \n", end="\r")

#%% 8. Calculating PSI scores

for event in inputs:
    
    "2a.) Casette Exons"
    
    if event.upper()=="CE":
        "open output file"
        if args.name:
            out=open(args.out+"PSI_table_"+args.name+".txt", "w")
        elif args.coordinates:
            out=open(args.out+"PSI_table_"+args.coordinates+".txt", "w")
        else:
            out=open(args.out+"PSI_table.txt", "w")
        title="Location\tDatabase\t"+"\t".join(sample_names)

            
        out.write(title+"\n")

        #Progress tracker
        total_iterations= len(gene_exons)*len(exons_bam)*len(list(gene_exons.values()))
        current_iteration=0
        percentage=round(current_iteration/total_iterations,2)
        #print("Counting Reads: ",percentage, "%", end="\r")

        for gene_id in gene_exons:
            #extract strand and chrom from first exon
            strand=gene_exons[gene_id][0][3]
            
            #Line for gene id, so that entries in table are separated.
            out.write("#"+gene_id+", "+strand+"\n")
            for exon in gene_exons[gene_id]:
                PSI_scores=[]
                chrom=exon[0]
                start=int(exon[1]) #Note that these are start and stop coor
                stop=int(exon[2]) #not start and stop of exon. Strand dependent.
                strand=exon[3]
                        
                for sample in sample_names:
                    current_iteration+=1
                    percentage=round(current_iteration/total_iterations,2)
                    #print("Counting Reads: ",percentage, "%", end="\r")
                    #initiate count values for PSI scores
                    junction3=0
                    junction5=0
                    splice_junction=0
                    
                    for read in exons_bam[sample]:
                        smaller_exon=read[1]
                        bigger_exon=read[2]
                        
                        #Exclude read if it is on the wrong strand.
                        if read[3]!=strand:
                            continue
                        
                        "Counts:"
                        if int(smaller_exon[1])<start and int(bigger_exon[0])>stop:
                            splice_junction+=1
                        elif int(smaller_exon[1])<start and \
                            start<=int(bigger_exon[1]) <= stop:
                            junction5+=1
                        elif int(bigger_exon[0])>stop and \
                            start<=int(smaller_exon[0])<=stop:
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
                loc=exon[0]+"_"+str(exon[1])+"_"+str(exon[2])
                out.write("{}\t{}\t{}\n".format(loc,exons_db[loc+"_"+strand], 
                                     "\t".join(PSI_scores)))


        out.close

print("Counting Reads: Done!         \n", end="\r")



print(genotype_HMZR)


#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))      

































