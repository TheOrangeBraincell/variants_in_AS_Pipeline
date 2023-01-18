# -*- coding: utf-8 -*-
"""
Date: Wed Jan 18 14:18:40 2023
File Name: PSI.py
Author: Mirjam Karlsson-Müller

Description:


Abbreviations:
    AS=     Alternative Splicing
    CE=     Casette Exon
    AA=     Alternative Acceptor
    AD=     Alternative Donor
    IR=     Intron Retention
    
List of Functions:
    Filter_Reads(): 
        Basic filtering done for all .bam file entries before processing. Reads with non-primary alignments
        as well as reads on the wrong strand, are filtered out.
    
    PSI_CE(): 
        Calculates the PSI score for one sample and potential CE event. Also requires gene name as input.
        Returns PSI score as string for easier writing and allowing for NAN
    
    PSI_AA():
        Calculates the PSI score for one sample and start in AA event. Also requires gene name as input.
        Returns PSI score as string for easier writing and allowing for NAN
    
    PSI_AD():
        Calculates the PSI score for one sample and stop in AD event. Also requires gene name as input.
        Returns PSI score as string for easier writing and allowing for 
    
    PSI_IR():
        Calculates the PSI score for one sample and potential IR event. Also requires gene name as input.
        Returns PSI score as string for easier writing and allowing for NAN
    
Procedure: 
    1.
    2.
    3.
    
Useage:
    
    
Possible Bugs:
"""

#%% Imports

import argparse
import glob
import re
import time
import pysam
import os
import math
from inspect import currentframe, getframeinfo

#%% Time

start_time=time.time()

#%% 0. argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Detect and Score Alternative Splicing',
                                 usage='%(prog)s -s SAMPLE-FOLDER -o OUTPUT-FOLDER \
                                     -g GENCODE-FILE -r REFSEQ-FILE \
                                         [-c] "chrX:XXXXXX-XXXXXX" -as AS-TYPE \
                                             -is INSERT-SIZE',
                                 description="""Per AS event of interest, creates
                                 a table with the PSI scores supporting said event
                                 per sample in sample folder.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing among \
                        others, the vcf, bam and gene.tsv files.')
parser.add_argument('--input', '-i', required=True,
                    help="Table with all identified AS events, created by Identify_AS.py")
parser.add_argument('--out', '-o', required=True,
                    help="""Output folder, containing PSI tables.""")
parser.add_argument('--InsertSize', '-is', type=str,
                    help="""Average Insert size plus standard deviation. Format
                    'Mean X Standard Deviation Y' """)


args = parser.parse_args()


#If we need to calculate scores for intron retention, we require an average insert size between the reads.
if args.InsertSize:
    try:
        insert_mean=float(args.InsertSize.split(" ")[1])
        insert_sd=float(args.InsertSize.split(" ")[4].strip("\n"))
    except:
        raise argparse.ArgumentTypeError("""Insert Sizes are either missing or
                                         of wrong format. Required format is:
                                             'Mean [float] Standard Deviation [float]'""")

#Create output directory if not there already
if not os.path.exists(args.out):
    os.makedirs(args.out)
    
#%% User Defined Functions

def Filter_Reads(read, gene_strand):
    skip=False
    #Exclude non-primary alignments
    if read.is_secondary:
        skip=True
    
    #Exclude reads on wrong strand
    if read.mate_is_reverse and read.is_read1:
        read_strand="-"
    elif read.mate_is_reverse and read.is_read2:
        read_strand="+"
    elif read.mate_is_forward and read.is_read1:
        read_strand="+"
    elif read.mate_is_forward and read.is_read2:
        read_strand="-"
    
    if read_strand!=gene_strand:
        skip=True
        
    return skip

def PSI_CE(sample, CE, gene):
    """
    

    Parameters
    ----------
    sample : string
        Sample name
    exon : list
        contains chromosome, start, stop, strand information for CE

    Returns
    -------
    str
        PSI score

    """
    #Needs coordinates for CE
    chrom=CE[0]
    start=int(CE[1])
    stop=int(CE[2])
    strand=CE[3]                                                                        
    
    #Find .bam file corresponding to sample name.
    for file in bam_file_list:
        if sample in file:
            bam_file=file
            index_file=file[0:-1]+"i"
            break
    
    #We only want to count each read of a read pair once. So we have a list of reads that have already been counted.
    counted=dict()
    #Note that we prioritize counting spliced reads, so those are checked first.
    
    #Initialize opening of file
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    
    #For the spliced reads, we need to fetch a bigger region and then recognise the spliced alignments.
    #We dont want to miss any reads, so we make this region the gene range the exon is in.
    spliced_reads=samfile.fetch(chrom, gene_ranges[gene][2], gene_ranges[gene][3])
    
    #initialize counters
    counter_left=0
    counter_right=0
    counter_accross=0
    for read in spliced_reads:
        #Filter
        if Filter_Reads(read, strand)==True:
            continue
        #only spliced reads
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        #Allow for several splice junctions in one read.
        current_cigar = read.cigarstring
        read_start=int(read.reference_start)
        #sum all numbers in cigar string for read length
        read_range=sum([int(i) for i in re.findall(r'\d+', read.cigarstring)])
        current_start = int(read_start)
        #Do not process spliced reads in the wrong place.
        if read_start < start and read_start+read_range < start:
            continue
        if read_start > stop and read_start+read_range > stop:
            continue
        read_name=read.query_name
    
        while re.search(r'\d+M\d+N\d+M', current_cigar):
            #assign splice junction variables
            junction = re.search(r'(\d+)M(\d+)N(\d+)M', current_cigar)
            exon1 = int(junction.group(1))
            intron = int(junction.group(2))
            exon2 = int(junction.group(3))
            exon1_start = current_start
            exon1_end = exon1_start+exon1+1  #exclusive
            exon2_start = exon1_end+intron -1 #inclusive
            exon2_end=exon2_start+exon2+1
                
            #skip alignments with less than 3 matching bases in an exon.
            if exon1<3 or exon2<3:
                # update cigar string
                current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                current_start= exon2_start
                continue
            
            #Number of spliced reads from previous to CE.
            if exon1_end< start and stop >= exon2_end>=start:
                #count!
                if read_name not in counted:
                    counter_left+=1
                    counted[read_name]="spliced"
                else:
                    if counted[read_name]=="spliced":
                        #has already been counted. proceed.
                        pass
                    else:
                        #this shouldnt happen....
                        print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                        pass
            
            #Number of spliced reads from CE to next.   
            elif start <= exon1_start<=stop and exon2_start >stop:
                #count!
                if read_name not in counted:
                    counter_right+=1
                    counted[read_name]="spliced"
                else:
                    if counted[read_name]=="spliced":
                        #has already been counted. proceed.
                        pass
                    else:
                        #this shouldnt happen....
                        print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                        pass
                    
            #Number of spliced reads accross CE.
            elif exon1_end<start and exon2_start>stop:
                #count!
                if read_name not in counted:
                    counter_accross+=1
                    counted[read_name]="spliced"
                else:
                    if counted[read_name]=="spliced":
                        #has already been counted. proceed.
                        pass
                    else:
                        #this shouldnt happen....
                        print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                        pass
                    
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    
    #C= number of reads within CE.
    #For this we fetch the part of the alignment file in the exon.
    exon_reads=samfile.fetch(chrom, start, stop)
    #The coordinates do not need to be adjusted as pysam uses also start incl and end excl and 0 based.
    counter=0
    for read in exon_reads:
        #Filter
        if Filter_Reads(read, strand)==True:
            continue
        #no spliced reads
        if re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        read_name=read.query_name
        if read_name not in counted:
            counted[read_name]="unspliced"
            #Theres no real problem if both reads are in the CE, but we dont want to count them twice. So the counter is indented.
            counter+=1
        else:
            if counted[read_name]=="spliced":
                #spliced reads > unspliced reads for counts. So we dont count this one.
                continue
            elif counted[read_name]=="unspliced":
                #We have already counted it. Ignore and proceed.
                continue
            else:
                #this shouldnt happen....
                print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                continue
    
    #counter needs to be normalized with length of the exon.
    C=counter/(stop-start)
    

    #Initialize spliced flag for PSI_CE_dict
    if counter_left!=0 or counter_right!=0:
        spliced_flag=True
    else:
        spliced_flag=False
    
    IR= counter_left+counter_right+C
    ER= counter_accross
    if ER+ IR >10:
        PSI=str(round(IR/(IR+ER),3))
    else:
        PSI="NAN"
    
    return [PSI,spliced_flag]

def PSI_AA(gene, sample, event, start):
    """
    

    Parameters
    ----------
    sample : string
        Sample name
    event : list
        contains all starts for this particular AA event.
    start: 
        start that we currently calculate PSI score for.
    Returns
    -------
    str
        PSI score

    """
    #Extract gene information
    chrom, strand, gene_start, gene_stop=gene_ranges[gene]
    
    #Find .bam file corresponding to sample name.
    for file in bam_file_list:
        if sample in file:
            bam_file=file
            index_file=file[0:-1]+"i"
            break
        
    #We only want to count each read of a read pair once per event. So we have a list of reads that have already been counted.
    counted=dict()
    for i in range(0, len(event)):
        counted[i]=dict()
    #Note that we prioritize counting spliced reads, so those are checked first.
    
    #Initialize opening of file
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    
    #Count reads
    #Find total number of spliced reads with coordinates
    spliced_counters=[0]*len(event)
    #Open bam file in gene range.
    spliced_reads=samfile.fetch(chrom, gene_start, gene_stop)
    
    #Go through reads
    for read in spliced_reads:
        #Filter
        if Filter_Reads(read, strand)==True:
            continue
        #only spliced reads
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        
        #Allow for several splice junctions in one read.
        current_cigar = read.cigarstring
        current_start = int(read.reference_start)
        read_name=read.query_name
        
        while re.search(r'\d+M\d+N\d+M', current_cigar):
            #assign splice junction variables
            junction = re.search(r'(\d+)M(\d+)N(\d+)M', current_cigar)
            exon1 = int(junction.group(1))
            intron = int(junction.group(2))
            exon2 = int(junction.group(3))
            exon1_start = current_start
            exon1_end = exon1_start+exon1+1  #exclusive
            exon2_start = exon1_end+intron -1 #inclusive
                
            #skip alignments with less than 3 matching bases in an exon.
            if exon1<3 or exon2<3:
                # update cigar string
                current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                current_start= exon2_start
                continue
            
            #Update counters based on matching starts
            if strand=="+":
                #Then the AA starts have to match the exon2_starts
                for i in range(0, len(event)):
                    if event[i]==str(exon2_start):
                        if read_name not in counted[i]:
                            spliced_counters[i]+=1
                            counted[i][read_name]="spliced"
                        else:
                            if counted[i][read_name]=="spliced":
                                #Great its already counted. skip.
                                pass
                            else:
                                #this shouldnt happen....
                                print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                                pass
            else:
                #Then the AA starts have to match the exon1_ends
                for i in range(0, len(event)):
                    #Same as for AD on the plus strand, the exon end has to be adjusted, because we calculate it to be exclusive.
                    if event[i]==str(exon1_end-1):
                        if read_name not in counted[i]:
                            spliced_counters[i]+=1
                            counted[i][read_name]="spliced"
                        else:
                            if counted[i][read_name]=="spliced":
                                #Great its already counted. skip.
                                pass
                            else:
                                #this shouldnt happen....
                                print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                                pass

            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    
    #Find difference reads
    difference_counters=[]
    difference_lengths=[]
    for i in range(0, len(event)-1):
        start1=int(event[i])
        start2=int(event[i+1])
        length=max(start1, start2)-min(start1, start2)
        difference_reads=samfile.fetch(chrom, start1, start2)
        counter=0
        
        #Filter reads
        for read in difference_reads:
            #Filter
            if Filter_Reads(read, strand)==True:
                continue
            #no spliced reads
            if re.search(r'\d+M\d+N\d+M',read.cigarstring):
                continue
            
            #Filter conditions passed:
            read_name=read.query_name
            if read_name not in counted[i]:
                counter+=1
                counted[i][read_name]="unspliced"
            else:
                if counted[i][read_name]== "unspliced":
                    #no problem. already counted. skip
                    continue
                elif counted[i][read_name]=="spliced":
                    #Takes priority. already counted. skip.
                    continue
                else:
                    #this shouldnt happen....
                    print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                    continue
            
        #Add difference reads into list
        difference_counters.append(counter)
        difference_lengths.append(length)
    
    #Note that the first start, on the minus strand, has the highest coordinate.
    start_number=event.index(start)
    if start_number!=0 and strand=="+":
        total_reads=difference_counters[start_number-1]+spliced_counters[start_number]+sum([x for i,x in enumerate(difference_counters) if i!=(start_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
        if total_reads>10:
            #Then calculate PSI
            for i, x in enumerate(difference_counters):
                difference_counters[i]=x/difference_lengths[i]
            IR=difference_counters[start_number-1]+spliced_counters[start_number]
            ER=sum([x for i,x in enumerate(difference_counters) if i!=(start_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
    elif start_number==0 and strand=="+":
        total_reads=spliced_counters[start_number]+sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
        if total_reads>10:
            #Then calculate PSI
            for i, x in enumerate(difference_counters):
                difference_counters[i]=x/difference_lengths[i]
            IR=spliced_counters[start_number]
            ER= sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
    elif start_number!=(len(event)-1) and strand=="-":
        total_reads=difference_counters[start_number-1]+spliced_counters[start_number]+sum([x for i,x in enumerate(difference_counters) if i!=(start_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
        if total_reads>10:
            #Then calculate PSI
            for i, x in enumerate(difference_counters):
                difference_counters[i]=x/difference_lengths[i]
            IR=difference_counters[start_number-1]+spliced_counters[start_number]
            ER=sum([x for i,x in enumerate(difference_counters) if i!=(start_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
    else:
        total_reads=spliced_counters[start_number]+sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)])
        if total_reads>10:
            #Then calculate PSI
            for i, x in enumerate(difference_counters):
                difference_counters[i]=x/difference_lengths[i]
            IR=spliced_counters[start_number]
            ER= sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
    
    if total_reads >10:
        PSI=str(round(IR/(IR+ER),3))
    else:
        PSI="NAN"
        
    return PSI

def PSI_AD(gene, sample, event, stop):
    """
    
    Parameters
    ----------
    sample : string
        Sample name
    event : list
        contains all starts for this particular AA event.
    start: 
        start that we currently calculate PSI score for.
    Returns
    -------
    str
        PSI score

    """
    #Extract gene information
    chrom, strand, gene_start, gene_stop=gene_ranges[gene]
    
    #Find .bam file corresponding to sample name.
    for file in bam_file_list:
        if sample in file:
            bam_file=file
            index_file=file[0:-1]+"i"
            break
    #Initialize opening of file
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    
    #We only want to count each read of a read pair once per event. So we have a list of reads that have already been counted.
    counted=dict()
    for i in range(0, len(event)):
        counted[i]=dict()
    #Note that we prioritize counting spliced reads, so those are checked first.
    
    #Count reads
    #Find total number of spliced reads with coordinates
    spliced_counters=[0]*len(event)
    #Open bam file in gene range.
    spliced_reads=samfile.fetch(chrom, gene_start, gene_stop)
    
    #Go through reads
    for read in spliced_reads:
        #Filter
        if Filter_Reads(read, strand)==True:
            continue
        #only spliced reads
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        
        #Allow for several splice junctions in one read.
        current_cigar = read.cigarstring
        current_start = int(read.reference_start)
        read_name=read.query_name
        
        while re.search(r'\d+M\d+N\d+M', current_cigar):
            #assign splice junction variables
            junction = re.search(r'(\d+)M(\d+)N(\d+)M', current_cigar)
            exon1 = int(junction.group(1))
            intron = int(junction.group(2))
            exon2 = int(junction.group(3))
            exon1_start = current_start
            exon1_end = exon1_start+exon1+1  #exclusive
            exon2_start = exon1_end+intron -1 #inclusive
                
            #skip alignments with less than 3 matching bases in an exon.
            if exon1<3 or exon2<3:
                # update cigar string
                current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                current_start= exon2_start
                continue
            
            #Update counters based on matching starts
            if strand=="+":
                #Then the AA starts have to match the exon2_starts
                for i in range(0, len(event)):
                    #The coordinates on plus strand seem to be shifted for AD, by 1. Because we calculate it to be exclusive.
                    if str(event[i])==str(exon1_end-1):
                        if read_name not in counted[i]:
                            spliced_counters[i]+=1
                            counted[i][read_name]="spliced"
                        else:
                            if counted[i][read_name]=="spliced":
                                #Great its already counted. skip.
                                pass
                            else:
                                #this shouldnt happen....
                                print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                                pass
            else:
                #Then the AA starts have to match the exon1_ends
                for i in range(0, len(event)):
                    if str(event[i])==str(exon2_start):
                        if read_name not in counted[i]:
                            spliced_counters[i]+=1
                            counted[i][read_name]="spliced"
                        else:
                            if counted[i][read_name]=="spliced":
                                #Great its already counted. skip.
                                pass
                            else:
                                #this shouldnt happen....
                                print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                                pass
            
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    
    #Find difference reads
    difference_counters=[]
    difference_lengths=[]
    for i in range(0, len(event)-1):
        stop1=int(event[i])
        stop2=int(event[i+1])
        length=max(stop1,stop2)-min(stop1,stop2)
        difference_reads=samfile.fetch(chrom, stop1, stop2)
        counter=0
        
        #Filter reads
        for read in difference_reads:
            #Filter
            if Filter_Reads(read, strand)==True:
                continue
            #no spliced reads
            if re.search(r'\d+M\d+N\d+M',read.cigarstring):
                continue
            
            #Filter conditions passed:
            read_name=read.query_name
            if read_name not in counted[i]:
                counter+=1
                counted[i][read_name]="unspliced"
            else:
                if counted[i][read_name]== "unspliced":
                    #no problem. already counted. skip
                    continue
                elif counted[i][read_name]=="spliced":
                    #Takes priority. already counted. skip.
                    continue
                else:
                    #this shouldnt happen....
                    print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                    continue
        
        #Add difference reads into list
        difference_counters.append(counter)
        difference_lengths.append(length)
        
    #Note that the first start, on the minus strand, has the highest coordinate.
    stop_number=event.index(stop)
    if stop_number!=0 and strand=="+":
        total_reads=difference_counters[stop_number-1]+spliced_counters[stop_number]+sum([x for i,x in enumerate(difference_counters) if i!=(stop_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
        if total_reads>10:
            #Then calculate PSI
            for i, x in enumerate(difference_counters):
                difference_counters[i]=x/difference_lengths[i]
            IR=difference_counters[stop_number-1]+spliced_counters[stop_number]
            ER=sum([x for i,x in enumerate(difference_counters) if i!=(stop_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
            
    elif stop_number==0 and strand=="+":
        total_reads=spliced_counters[stop_number]+sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
        if total_reads>10:
            #Then calculate PSI
            for i, x in enumerate(difference_counters):
                difference_counters[i]=x/difference_lengths[i]
            IR=spliced_counters[stop_number]
            ER= sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
    elif stop_number!=(len(event)-1) and strand=="-":
        total_reads=difference_counters[stop_number]+spliced_counters[stop_number]+sum([x for i,x in enumerate(difference_counters) if i!=(stop_number)])+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
        if total_reads>10:
            #Then calculate PSI
            for i, x in enumerate(difference_counters):
                difference_counters[i]=x/difference_lengths[i]
            IR=difference_counters[stop_number]+spliced_counters[stop_number]
            ER=sum([x for i,x in enumerate(difference_counters) if i!=(stop_number)])+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
    else:
        total_reads=spliced_counters[stop_number]+ sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
        if total_reads>10:
            #Then calculate PSI
            for i, x in enumerate(difference_counters):
                difference_counters[i]=x/difference_lengths[i]
            IR=spliced_counters[stop_number]
            ER= sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(stop_number)]) 
    
    if total_reads>10:
        PSI=str(round(IR/(IR+ER),3))
    else:
        PSI="NAN"

    return PSI
    
def PSI_IR(sample, entry, gene):
    """

    Parameters
    ----------
    sample : string
        Sample name
    entry: string separated by _
        contains strand, chromosome, exon1stop, exon2start

    Returns
    -------
    str
        PSI score

    """
    strand, chrom, exon1stop, exon2start=entry.split("_")
    exon1stop=int(exon1stop)
    exon2start=int(exon2start)
    #Find .bam file corresponding to sample name.
    for file in bam_file_list:
        if sample in file:
            bam_file=file
            index_file=file[0:-1]+"i"
            break
    
    #Check: If there is a CE in the IR we need to check its PSI score
    for gene in PSI_CE_dict:
        #Extract start and stop
        for entry in PSI_CE_dict[gene]:
            CE_start=int(entry.split("_")[1])
            CE_stop=int(entry.split("_")[2])
            #Check if within intron
            if CE_start>exon1stop and CE_stop<exon2start:
                #Then we check PSI
                CE_PSI=PSI_CE_dict[gene][entry][sample][0]
                spliced_flag=PSI_CE_dict[gene][entry][sample][1]
                
                #If this PSI=1 and spliced_flag=True, then theres no evidence for intron retention
                if CE_PSI=="1.0" and spliced_flag==True:
                    #print("It happend for "+ sample+ " CE coordinates: "+ str(CE_start)+", "+ str(CE_stop)+" IR coordinates: "+ str(exon1stop)+", "+ str(exon2start))
                    PSI="NAN"
                    return PSI
                    
    #Insert read confidence interval
    insert_confidence=[insert_mean-1.96*(insert_sd/math.sqrt(len(sample_names))),insert_mean+1.96*(insert_sd/math.sqrt(len(sample_names)))]
    
    #We only want to count each read of a read pair once per event. So we have a list of reads that have already been counted.
    counted=dict()
    #Note that we prioritize counting spliced reads, so those are checked first.
    
    #Initialize opening of file
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    
    #Get Spliced Reads
    #Open bam file in gene range.
    spliced_reads=samfile.fetch(chrom, gene_ranges[gene][2], gene_ranges[gene][3])
    
    spliced_counter=0
    for read in spliced_reads:
        #Filter
        if Filter_Reads(read, strand)==True:
            continue
        #only spliced reads
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue

        read_start=int(read.reference_start)
        read_length=int(read.infer_query_length())
        read_name=read.query_name
        #Allow for several splice junctions in one read.
        current_cigar = read.cigarstring
        current_start = int(read_start)
        
        while re.search(r'\d+M\d+N\d+M', current_cigar):
            #assign splice junction variables
            junction = re.search(r'(\d+)M(\d+)N(\d+)M', current_cigar)
            exon1 = int(junction.group(1))
            intron = int(junction.group(2))
            exon2 = int(junction.group(3))
            exon1_start = current_start
            exon1_end = exon1_start+exon1+1  #exclusive
            exon2_start = exon1_end+intron -1 #inclusive
                
            #skip alignments with less than 3 matching bases in an exon.
            if exon1<3 or exon2<3:
                # update cigar string
                current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
                current_start= exon2_start
                continue
            
            if strand=="+":
                if exon1_end<=exon1stop and exon2_start>=exon2start:
                    if read_name not in counted:
                        spliced_counter+=1
                        counted[read_name]="spliced"
                    else:
                        if counted[read_name]=="spliced":
                            #already counted. no need to double up.
                            pass
                        elif counted[read_name]=="left" or counted[read_name]=="right":
                            #this shouldnt happen....
                            print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                            pass
            
            else:
                if exon1_end<=exon1stop and exon2_start>=exon2start:
                    if read_name not in counted:
                        spliced_counter+=1
                        counted[read_name]="spliced"
                    else:
                        if counted[read_name]=="spliced":
                            #already counted. no need to double up.
                            pass
                        elif counted[read_name]=="left" or counted[read_name]=="right":
                            #this shouldnt happen....
                            print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                            pass
            
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    
    #Get Overlapping Reads
    #Fetch reads in gene, since we want them to overlap with beginning or end or IR, the reads cant be too far away from our coordinates.
    #so we use a generous range of pm 300.
    reads_genome=samfile.fetch(chrom,exon1stop-300, exon2start+300)
    overlap_counter=0

    #Only keep IR which have overlapping reads on both sides.
    left=False
    right=False
    for read in reads_genome:
        #Filter
        if Filter_Reads(read, strand)==True:
            continue
        #no spliced reads
        if re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        
        #if distance between reads too far, then this is a sign that the second read is in an exon and in between is a break. So no IR.
        if read.tlen>insert_confidence[1]:
            continue
        
        read_name=read.query_name
        read_start=int(read.reference_start)
        read_length=int(read.infer_query_length())
        #check coordinates, if they overlap the exon1end or the exon2start, it could be an overlap read.
        parts=re.search(r'((?:\d+[I,D,S,H])?)((?:\d+M)*)((?:\d+[I,D,S,H])?)((?:\d+M)?)((?:\d+[I,D,S,H])?)',read.cigarstring)
        #5 groups: first one before match, second match, third one in between matches, fourth one match, 5th following match.
        match_length=0
        after_match=0
        before_match=0
        read_stop=read_start+read_length
        #remove unaligned parts before match
        if parts.group(1):
            before_match=int(re.sub('[IDSH]', '', parts.group(1)))
        read_start=read_start+before_match
        
        #match group
        if parts.group(2):
            #add to match length
            match_length+=int(parts.group(2).strip("M"))
        
        #Potential part between matches:
        if parts.group(3):
            #Is last group if there is only one M, and middle if theres 2.
            if parts.group(4):
                #Different behavior for S, D, and I
                if re.search(r'[SI]', parts.group(3)):
                    #add to match length
                    match_length+=int(re.sub('[IDSH]', '', parts.group(3)))
                    #D does not get added, as it shortens the query.
            else:
                after_match=int(re.sub('[IDSH]', '', parts.group(3)))
        
        #Potential second match group
        if parts.group(4):
            #Add to match length
            match_length+=int(parts.group(4).strip("M"))
        
        #remove unaligned parts after match
        if parts.group(5):
            after_match=int(re.sub('[IDSH]', '', parts.group(5)))
        read_stop=read_stop-after_match
        
        #matching parts coordinates
        if read_stop-read_start!=match_length:
            print("Error calculating: ",read_stop-read_start, match_length, read.cigarstring)
            print(parts.group(1),parts.group(2),parts.group(3),parts.group(4),parts.group(5))
            quit()
        

        #Find the overlapping ones left side
        if read_start<=exon1stop and read_stop>exon1stop:
            if read_name not in counted:
                counted[read_name]="left"
                overlap_counter+=1
                left=True
            else:
                if counted[read_name]=="left":
                    #already counted. but we need a flag from left or right.
                    #flag already set if left, so we can leave it.
                    continue
                elif counted[read_name]=="right":
                    #Already counted, but we need to set the flag.
                    left=True
                elif counted[read_name]=="spliced":
                    #This should not happen, as the reads in a pair cannot pick up 2 versions of one event.
                    print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                    continue
        #right side
        elif read_start<exon2start and read_stop>=exon2start:
            if read_name not in counted:
                counted[read_name]="right"
                overlap_counter+=1
                right=True
            else:
                if counted[read_name]=="right":
                    #already counted. but we need a flag from left or right.
                    #flag already set if left, so we can leave it.
                    continue
                elif counted[read_name]=="left":
                    #Already counted, but we need to set the flag.
                    right=True
                elif counted[read_name]=="spliced":
                    #This should not happen, as the reads in a pair cannot pick up 2 versions of one event.
                    print("Theres a bug in your code. check line ", getframeinfo(currentframe()).lineno)
                    continue
    
    #If either side has no overlapping reads, then the PSI score will be NAN and the spliced reads do not need to be counted.
    if left==False or right==False:
        PSI="NAN"
        return PSI
    
    #Calculate PSI
    IR=overlap_counter
    ER=spliced_counter
    if IR+ER<=10:
        PSI="NAN"
    else:
        PSI=str(round(IR/(IR+ER),2))

    return PSI

#%% Find all samples and their alignment files

#All bam files.
argument_glob=args.samples+"/**/alignment.bam"
bam_file_list=glob.glob(argument_glob, recursive=True)

#Assuming that each bam.file is in its own sample folder:
sample_names=[i.split("/")[-4] for i in bam_file_list]

#%% Score Splicing events

#Open input file
with open(args.input,"r") as infile:
    for line in infile:
        #Find gene headers
        if line.startswith("#"):
            current_gene=line.strip("\n").split(",")
        #If there is no header, it is an AS event entry
        #These are sorted by alphabet. So it is always AA, AD, CE, IR
        else:
            
















































