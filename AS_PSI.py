# -*- coding: utf-8 -*-
"""
Date: Mon Nov  7 11:12:17 2022
File Name: AS_PSI.py
Author: Mirjam Karlsson-Müller

Description:
    
    
Abbreviations:
    AS=     Alternative Splicing
    CE=     Casette Exon
    AA=     Alternative Acceptor
    AD=     Alternative Donor
    IR=     Intron Retention
    
Procedure: 
    1.
    2.
    3.

Modules:
    

List of Functions:
    
    
Useage:

    Inputs:
        - Sample Folder
        - Database files for GENCODE and RefSeq respectively
        - output file name
        - evtl. coordinates or name of a gene/region of interest
        - What type of alternative splicing event?

    Instructions:
        Run in command line.
        
        #with coordinates f.e. Estrogen Receptor
        python AS_PSI.py -s ../Sample_Data/ -o PSI_ESR1/ -c "chr6:151656691-152129619" -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as ALL
        
        #With coordinates f.e. BRCA1 (neg strand)
        python AS_PSI.py -s ../Sample_Data/ -o PSI_BRCA1/ -c "chr17:43044295-43170245" -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as ALL
        
        
        #for server
        python AS_PSI.py -s /raidset/mi3258mu-se/Mirjam -o PSI_ESR1.tsv -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as CE -c "chr6:151690496-152103274"
        

        #no coordinates/name
        
Possible Bugs:
"""

#%% Imports

import argparse
import glob
import re
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
                    help='folder containing sample folders containing among \
                        others, the vcf, bam and gene.tsv files.')
parser.add_argument('--out', '-o', required=True,
                    help="""Output folder, containing PSI table.""")
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
        #To adjust for different inclusivity of stop/start for plus and minus strand, expand coordinate range by 1.
        coord_start = int(coord.group(2))-1
        coord_stop = int(coord.group(3))+1

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
            
    #If input is all events, make sure the code runs through all:
    if inputs[0].upper()=="ALL":
        inputs=["CE", "AA", "AD", "IR"]
if args.name and args.coordinates:
    raise argparse.ArgumentError("""Cannot process input name and coordinates. 
                                 Please use only one or the other.""")
    quit()

#Create output directory if not there already
if not os.path.exists(args.out):
    os.makedirs(args.out)
    
    
#%% Testing new stuff

"""
#How to progress bar.
n=70
print("Testing")
for i in range(0,n+1):
    progress=int(i*100/(n))
    print("#"*progress+" "*(100-progress)+" "+str(round(i/(n) *100, 2)) + "%", end="\r")

print("#"*progress+" "*(50-progress)+" "+str(round(i/(n) *100, 2)) + "%\n", end="\r")
"""


#%% User Defined Functions


def add_to(events, AA_counter, AD_counter):
    for e in events:
        
        #We dont want to add a last exon for AD or a first one for AA.
        if e=="AA":
            if entry[4]=="first" or coord_exons[key_string][4]=="first":
                continue
            #We want to save the actual exon starts (not just smaller coordinate.)
            strand=entry[3]
            if strand=="+":
                start1=int(entry[1])
                start2=int(coord_exons[key_string][1])
            
            else:
                start1=int(entry[2])
                start2=int(coord_exons[key_string][2])
            #Go through previous entries in potential_AA to make sure theres no duplicates.
            event_exists=False
            for events in potential_AA[gene]:
                if start1 in events and start2 not in events:
                    events.append(start2)
                    event_exists=True
                elif start2 in events and start1 not in events:
                    events.append(start1)
                    event_exists=True
            
            #If the event is not found, make a new one.
            if event_exists==False:
                potential_AA[gene].append([int(start1), int(start2)])
                AA_counter+=1
        
            
        if e=="AD":
            if entry[4]=="last" or coord_exons[key_string][4]=="last":
                continue
            
            #save exon stops (actual, not just bigger coordinate)
            strand=entry[3]
            if strand=="+":
                stop1=int(entry[2])
                stop2=int(coord_exons[key_string][2])
            
            else:
                stop1=int(entry[1])
                stop2=int(coord_exons[key_string][1])
            #Go through previous entries in potential_AA to make sure theres no duplicates.
            event_exists=False
            for events in potential_AD[gene]:
                if stop1 in events and stop2 not in events:
                    events.append(stop2)
                    event_exists=True
                elif stop2 in events and stop1 not in events:
                    events.append(stop1)
                    event_exists=True
                elif stop1 in events and stop2 in events:
                    event_exists=True
                    
            #If the event is not found, make a new one.
            if event_exists==False:
                potential_AD[gene].append([stop1, stop2])
                AD_counter+=1
        
    return(AA_counter, AD_counter)

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
    #print(start, stop)
    
    #Find .bam file corresponding to sample name.
    for file in bam_file_list:
        if sample in file:
            bam_file=file
            index_file=file[0:-1]+"i"
            break
    
    #Initialize opening of file
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    
    #C= number of reads within CE.
    #For this we fetch the part of the alignment file in the exon.
    exon_reads=samfile.fetch(chrom, start, stop)
    #The coordinates do not need to be adjusted as pysam uses also start incl and end excl and 0 based.
    
    read_dict=dict()
    counter=0
    for read in exon_reads:
        #Filter
        # Exclude non-primary alignments
        if read.is_secondary:
            continue
        #exclude second read of pair, if maps to overlapping region.
        read_name=read.query_name
        read_start=int(read.reference_start)
        read_length=int(read.infer_query_length())
        if read_name in read_dict:
            if read_dict[read_name][0]<=read_start<=read_dict[read_name][1] or \
            read_dict[read_name][0]<=read_start+read_length<=read_dict[read_name][1]:
                continue
        else:
            read_dict[read_name]=[read_start, read_start+read_length] 
        counter+=1
    
    #counter needs to be normalized with length of the exon.
    C=counter/(stop-start)
    #print(counter)
    #quit()
    
    #For the spliced reads, we need to fetch a bigger region and then recognise the spliced alignments.
    #We dont want to miss any reads, so we make this region the gene range the exon is in.
    spliced_reads=samfile.fetch(chrom, gene_ranges[gene][2], gene_ranges[gene][3])
    
    #initialize counters
    counter_left=0
    counter_right=0
    counter_accross=0
    #reset read dictionary
    read_dict=dict()
    
    for read in spliced_reads:
        #only spliced reads
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        #exclude second read of pair, if maps to overlapping region.
        read_name=read.query_name
        read_start=int(read.reference_start)
        read_length=int(read.infer_query_length())
        if read_name in read_dict:
            if read_dict[read_name][0]<=read_start<=read_dict[read_name][1] or \
            read_dict[read_name][0]<=read_start+read_length<=read_dict[read_name][1]:
                continue
        else:
            read_dict[read_name]=[read_start, read_start+read_length] 
        
        "Get strand information, exclude reads on wrong strand"
        if read.mate_is_reverse and read.is_read1:
            read_strand="-"
        elif read.mate_is_reverse and read.is_read2:
            read_strand="+"
        elif read.mate_is_forward and read.is_read1:
            read_strand="+"
        elif read.mate_is_forward and read.is_read2:
            read_strand="-"
        
        if read_strand!=strand:
            continue
        
        #Allow for several splice junctions in one read.
        current_cigar = read.cigarstring
        current_start = int(read_start)
        #Do not process spliced reads in the wrong place.
        """
        if read_start < start and read_start+read_length < start:
            continue
        if read_start > stop and read_start+read_length > stop:
            continue
        """
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
            
            #if exon1_end < start:
            #    print(exon2_start, start)
            
            #D= Number of spliced reads from previous to CE.
            #print(exon2_start, start)
            if exon1_end< start and stop >= exon2_end>=start:
                #print(exon1_end, exon2_start)
                counter_left+=1
            #E= Number of spliced reads from CE to next.   
            elif start <= exon1_start<=stop and exon2_start >stop:
                counter_right+=1
                #print(exon1_end, exon2_start)
            #F= number of spliced reads accross CE.
            elif exon1_end<start and exon2_start>stop:
                counter_accross+=1
                #print(exon1_end, exon2_start)
                
            
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    
    IR= counter_left+counter_right+C
    ER= counter_accross
    if ER+ IR >10:
        PSI=str(round(IR/(IR+ER),3))
    else:
        PSI="NAN"
    
    return PSI

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
    #Initialize opening of file
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    
    #Count reads
    #Find total number of spliced reads with coordinates
    spliced_counters=[0]*len(event)
    #Open bam file in gene range.
    spliced_reads=samfile.fetch(chrom, gene_start, gene_stop)
    
    #Go through reads
    read_dict=dict()
    for read in spliced_reads:
        #only spliced reads
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        #exclude second read of pair, if maps to overlapping region.
        read_name=read.query_name
        read_start=int(read.reference_start)
        read_length=int(read.infer_query_length())
        if read_name in read_dict:
            if read_dict[read_name][0]<=read_start<=read_dict[read_name][1] or \
            read_dict[read_name][0]<=read_start+read_length<=read_dict[read_name][1]:
                continue
        else:
            read_dict[read_name]=[read_start, read_start+read_length] 
        
        "Get strand information, exclude reads on wrong strand"
        if read.mate_is_reverse and read.is_read1:
            read_strand="-"
        elif read.mate_is_reverse and read.is_read2:
            read_strand="+"
        elif read.mate_is_forward and read.is_read1:
            read_strand="+"
        elif read.mate_is_forward and read.is_read2:
            read_strand="-"
        
        if read_strand!=strand:
            continue
        
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
            
            #if sample=="S000001":
            #    print(start, exon1_end)
            #Update counters based on matching starts
            if strand=="+":
                #Then the AA starts have to match the exon2_starts
                for i in range(0, len(event)):
                    if event[i]==str(exon2_start):
                        spliced_counters[i]+=1
            else:
                #Then the AA starts have to match the exon1_ends
                for i in range(0, len(event)):
                    #Same as for AD on the plus strand, the exon end has to be adjusted, because we calculate it to be exclusive.
                    if event[i]==str(exon1_end-1):
                        spliced_counters[i]+=1
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    
    #Find difference reads
    difference_counters=[]
    for i in range(0, len(event)-1):
        start1=int(event[i])
        start2=int(event[i+1])
        length=max(start1, start2)-min(start1, start2)
        difference_reads=samfile.fetch(chrom, start1, start2)
        counter=0
        
        read_dict=dict()
        #Filter reads
        for read in difference_reads:
            #no spliced reads
            if re.search(r'\d+M\d+N\d+M',read.cigarstring):
                continue
            #exclude second read of pair, if maps to overlapping region.
            read_name=read.query_name
            read_start=int(read.reference_start)
            read_length=int(read.infer_query_length())
            if read_name in read_dict:
                if read_dict[read_name][0]<=read_start<=read_dict[read_name][1] or \
                read_dict[read_name][0]<=read_start+read_length<=read_dict[read_name][1]:
                    continue
            else:
                read_dict[read_name]=[read_start, read_start+read_length] 
            
            #Get strand information, exclude reads on wrong strand
            if read.mate_is_reverse and read.is_read1:
                read_strand="-"
            elif read.mate_is_reverse and read.is_read2:
                read_strand="+"
            elif read.mate_is_forward and read.is_read1:
                read_strand="+"
            elif read.mate_is_forward and read.is_read2:
                read_strand="-"
            
            if read_strand!=strand:
                continue
            
            #Filter conditions passed:
            counter+=1
        #Add difference reads into list
        difference_counters.append(round(counter/length,2))
        
    
    #Note that the first start, on the minus strand, has the highest coordinate.
    start_number=event.index(start)
    if start_number!=0 and strand=="+":
        IR=difference_counters[start_number-1]+spliced_counters[start_number]
        ER=sum([x for i,x in enumerate(difference_counters) if i!=(start_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
    elif start_number==0 and strand=="+":
        IR=spliced_counters[start_number]
        ER= sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
    elif start_number!=(len(event)-1) and strand=="-":
        IR=difference_counters[start_number-1]+spliced_counters[start_number]
        ER=sum([x for i,x in enumerate(difference_counters) if i!=(start_number-1)])+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
    else:
        IR=spliced_counters[start_number]
        ER= sum(difference_counters)+sum([x for i,x in enumerate(spliced_counters) if i!=(start_number)]) 
    
    if ER+ IR >10:
        PSI=str(round(IR/(IR+ER),3))
    else:
        PSI="NAN"
        
    #if sample=="S000001":
    #    print(start, spliced_counters, difference_counters)
    return PSI

def PSI_AD(gene, sample, exon, stop):
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
    
    #Count reads
    #Find total number of spliced reads with coordinates
    spliced_counters=[0]*len(event)
    #Open bam file in gene range.
    spliced_reads=samfile.fetch(chrom, gene_start, gene_stop)
    
    #Go through reads
    read_dict=dict()
    for read in spliced_reads:
        #only spliced reads
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        #exclude second read of pair, if maps to overlapping region.
        read_name=read.query_name
        read_start=int(read.reference_start)
        read_length=int(read.infer_query_length())
        if read_name in read_dict:
            if read_dict[read_name][0]<=read_start<=read_dict[read_name][1] or \
            read_dict[read_name][0]<=read_start+read_length<=read_dict[read_name][1]:
                continue
        else:
            read_dict[read_name]=[read_start, read_start+read_length] 
        
        "Get strand information, exclude reads on wrong strand"
        if read.mate_is_reverse and read.is_read1:
            read_strand="-"
        elif read.mate_is_reverse and read.is_read2:
            read_strand="+"
        elif read.mate_is_forward and read.is_read1:
            read_strand="+"
        elif read.mate_is_forward and read.is_read2:
            read_strand="-"
        
        if read_strand!=strand:
            continue
        
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
            #if sample=="S000001":
            #    print(stop, exon2_start)
            #Update counters based on matching starts
            if strand=="+":
                #Then the AA starts have to match the exon2_starts
                for i in range(0, len(event)):
                    #The coordinates on plus strand seem to be shifted for AD, by 1. Because we calculate it to be exclusive.
                    if str(event[i])==str(exon1_end-1):
                        spliced_counters[i]+=1
            else:
                #Then the AA starts have to match the exon1_ends
                for i in range(0, len(event)):
                    if str(event[i])==str(exon2_start):
                        spliced_counters[i]+=1
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
        
        read_dict=dict()
        #Filter reads
        for read in difference_reads:
            #no spliced reads
            if re.search(r'\d+M\d+N\d+M',read.cigarstring):
                continue
            #exclude second read of pair, if maps to overlapping region.
            read_name=read.query_name
            read_start=int(read.reference_start)
            read_length=int(read.infer_query_length())
            if read_name in read_dict:
                if read_dict[read_name][0]<=read_start<=read_dict[read_name][1] or \
                read_dict[read_name][0]<=read_start+read_length<=read_dict[read_name][1]:
                    continue
            else:
                read_dict[read_name]=[read_start, read_start+read_length] 
            
            #Get strand information, exclude reads on wrong strand
            if read.mate_is_reverse and read.is_read1:
                read_strand="-"
            elif read.mate_is_reverse and read.is_read2:
                read_strand="+"
            elif read.mate_is_forward and read.is_read1:
                read_strand="+"
            elif read.mate_is_forward and read.is_read2:
                read_strand="-"
            
            if read_strand!=strand:
                continue
            
            #Filter conditions passed:
            counter+=1
        #Add difference reads into list
        difference_counters.append(counter)
        difference_lengths.append(length)
    if sample=="S000003":
        print(event, difference_counters, spliced_counters)
        
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
    
    if sample=="S000003" and total_reads>10:
        print(IR, ER)
    if total_reads>10:
        PSI=str(round(IR/(IR+ER),3))
    else:
        PSI="NAN"
    if sample=="S000003":
        print(PSI)
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
    
    #Initialize opening of file
    samfile=pysam.AlignmentFile(bam_file, 'rb', index_filename=index_file)
    
    #Get Overlapping Reads
    #Fetch reads in gene
    reads_genome=samfile.fetch(chrom,gene_ranges[gene][2], gene_ranges[gene][3])
    read_dict=dict()
    overlap_counter=0
    
    for read in reads_genome:
        #Filter. First exclude non-primary alignments
        if read.is_secondary:
            continue
        #no spliced reads
        if re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        
        #Exclude second read of pair if they map to the same region.
        read_name=read.query_name
        read_start=int(read.reference_start)
        read_length=int(read.infer_query_length())
        if read_name in read_dict:
            if read_dict[read_name][0]<=read_start<=read_dict[read_name][1] or \
            read_dict[read_name][0]<=read_start+read_length<=read_dict[read_name][1]:
                continue
        else:
            read_dict[read_name]=[read_start, read_start+read_length] 
            
        #Get strand information, exclude reads on wrong strand
        if read.mate_is_reverse and read.is_read1:
            read_strand="-"
        elif read.mate_is_reverse and read.is_read2:
            read_strand="+"
        elif read.mate_is_forward and read.is_read1:
            read_strand="+"
        elif read.mate_is_forward and read.is_read2:
            read_strand="-"
        
        if read_strand!=strand:
            continue
        
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
        
        pass_on_read=False
        #Filter out reads that are within AD/AA events
        for event in potential_AA[gene]:
            if strand=="+":
                if read_start>=min(event):
                    pass_on_read=True
            if strand=="-":
                if read_stop<=max(event):
                    pass_on_read=True
        for event in potential_AD[gene]:
            if strand=="+":
                if read_stop<=max(event):
                    pass_on_read=True
            if strand=="-":
                if read_start>=min(event):
                    pass_on_read=True
        if pass_on_read==True:
            continue
        
        #Filter out reads that do not have their pair within an intron
        
        #Find the overlapping ones left side
        if read_start<=exon1stop and read_stop>exon1stop:
            overlap_counter+=1
        #right side
        elif read_start<exon2start and read_stop>=exon1stop:
            overlap_counter+=1
    
    #Get Spliced Reads
    #Open bam file in gene range.
    spliced_reads=samfile.fetch(chrom, gene_ranges[gene][2], gene_ranges[gene][3])
    
    #Go through reads
    read_dict=dict()
    spliced_counter=0
    for read in spliced_reads:
        #only spliced reads
        if not re.search(r'\d+M\d+N\d+M',read.cigarstring):
            continue
        #exclude second read of pair, if maps to overlapping region.
        read_name=read.query_name
        read_start=int(read.reference_start)
        read_length=int(read.infer_query_length())
        if read_name in read_dict:
            if read_dict[read_name][0]<=read_start<=read_dict[read_name][1] or \
            read_dict[read_name][0]<=read_start+read_length<=read_dict[read_name][1]:
                continue
        else:
            read_dict[read_name]=[read_start, read_start+read_length] 
        
        "Get strand information, exclude reads on wrong strand"
        if read.mate_is_reverse and read.is_read1:
            read_strand="-"
        elif read.mate_is_reverse and read.is_read2:
            read_strand="+"
        elif read.mate_is_forward and read.is_read1:
            read_strand="+"
        elif read.mate_is_forward and read.is_read2:
            read_strand="-"
        
        if read_strand!=strand:
            continue
        
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
                if exon1_end-1<=exon1stop and exon2_start>=exon2start:
                    spliced_counter+=1
            
            else:
                if exon1_end<=exon1stop and exon2_start>=exon2start:
                    spliced_counter+=1
            
            # update cigar string
            current_cigar = re.sub(r'^.*?N', 'N', current_cigar).lstrip("N")
            current_start= exon2_start
    #if sample=="S000001":
    #    print(entry, overlap_counter, spliced_counter)
    #Calculate PSI
    IR=overlap_counter
    ER=overlap_counter+spliced_counter
    if IR+ER<=10:
        PSI="NAN"
    else:
        PSI=str(round(IR/(IR+ER),2))
    
    return PSI

#%% 1. Process Database input

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
                
                #if its a new gene symbol, initialize inner dictionary
                if gene_name not in list(gene_dict.keys()):
                    gene_dict[gene_name]={trans_ID:[]}
                    
                else:
                    #add transcript ID to gene_dict dictionary
                    gene_dict[gene_name][trans_ID]=[]
                
                
                #make entries for each exon.
                for i in range(0, number_exons):
                    if i==0:
                        if strand=="+":
                            position= "first"
                        else:
                            position="last"
                    elif i==number_exons-1:
                        if strand=="+":
                            position="last"
                        else:
                            position="first"
                    else:
                        position="middle"
                    gene_dict[gene_name][trans_ID].append([chrom, 
                                                           exon_smaller[i], 
                                                           exon_bigger[i], 
                                                           strand, position,db])

print("Creating Database Dictionary: Done! \n", end="\r")

#Make gene ranges dictionary to extract right read region later.
gene_ranges=dict()
for gene in gene_dict:
    starts=[]
    stops=[]
    for trans_ID in gene_dict[gene]:
        starts.append(int(gene_dict[gene][trans_ID][0][1]))
        stops.append(int(gene_dict[gene][trans_ID][-1][2]))
        chrom=gene_dict[gene][trans_ID][0][0]
        strand=gene_dict[gene][trans_ID][0][3]
    #Take smallest start and biggest stop as range of gene.
    gene_ranges[gene]=[chrom, strand, min(starts), max(stops)]

#print(gene_ranges)
#%% SAMPLE NAMES?

#All bam files.
argument_glob=args.samples+"/**/alignment.bam"
bam_file_list=glob.glob(argument_glob, recursive=True)

#Assuming that each bam.file is in its own sample folder:
sample_names=[i.split("/")[-4] for i in bam_file_list]


#%% Identify alternative splicing events.

#To not find AA and AD twice:
AA_AD_dict=False

for event in inputs:
    
    "Casette Exons"
    
    if event.upper()=="CE":
        print("\nFinding Casette Exons in annotated Transcripts: \n", end="\r")
        #if a transcript has 2 or less exons, there can be no CE.
        filtered_gene_dict=dict()
        for gene in gene_dict:
            filtered_gene_dict[gene]={k : v for k , v in gene_dict[gene].items() 
                                if len(v) >2}
        print("#"*25+ " "*75+ "  25 %", end="\r")
        #print(filtered_gene_dict)
        #Remove first and last exon for each transcript, as they cannot be CE
        new_gene_dict=dict()
        for gene in filtered_gene_dict:
            if len(filtered_gene_dict[gene].items())>0:
                new_gene_dict[gene]=dict()
                for trans_ID in filtered_gene_dict[gene]:
                    new_gene_dict[gene][trans_ID]=[i for i in filtered_gene_dict[gene][trans_ID] if i[4]=="middle"]  
                
        CE_gene_dict=new_gene_dict
        print("#"*50+ " "*50+ " 50 %", end="\r")
        #print(CE_gene_dict)
        #Transcript IDs are no longer required, and additionally there
        #is definitely duplicate exons between R and G as well as within each
        #of them. Those need to be removed.
        
        #Creates dictionary with potential casette exons per gene.
        gene_exons=dict()
        #Creates a dictionary where each exon is listed with its database (in case needed later)
        exons_db=dict()
        
        CE_counter=0
        for gene in CE_gene_dict:
            if CE_gene_dict[gene]:
                gene_exons[gene]=[]
                for trans_ID in CE_gene_dict[gene]:
                    for exon in CE_gene_dict[gene][trans_ID]:
                        if exon[0:-1] not in gene_exons[gene]:
                            gene_exons[gene].append(exon[0:-1])
                            CE_counter+=1
                        if "_".join(exon[0:-1]) in exons_db:
                            if exon[-1] in exons_db["_".join(exon[0:-1])]:
                                continue
                            else:
                                exons_db["_".join(exon[0:-1])]+=exon[-1]
                        else:
                            exons_db["_".join(exon[0:-1])]=exon[-1]
            
        print("#"*100+ " 100 %\n", end="\r")
        print("Calculating PSI scores for Casette Exons: \n", end="\r")
        out=open(args.out+"PSI_CE.tsv", "w")
        #Make header for output file
        out.write("#PSI Table CE \n")
        out.write("Location\t"+"\t".join(sample_names)+"\n")
        #progress counter
        exons_done=0
        #Go through all potential casette exons, and calculate their PSI scores.
        for gene in gene_exons:
            #Extract strand information for header (for each gene)
            strand=gene_ranges[gene][1]
            
            #Sort exons based on start coordinate.
            gene_exons[gene].sort(key=lambda x: int(x[1]))
            #Create gene header
            out.write("#"+gene+", "+strand+"\n")
            for entry in gene_exons[gene]:
                PSI_scores=[]
                for sample in sample_names:
                    PSI_scores.append(PSI_CE(sample, entry, gene))
                    
                exons_done+=1
                out.write("_".join(entry[0:3])+"\t" +"\t".join(PSI_scores)+"\n")
                #Progress update
                percentage=round(exons_done/CE_counter *100, 2)
                print("#"*int(percentage)+" "*(100-int(percentage))+" " + str(percentage)+" %", end="\r")
        #Finish the bar.
        print("#"*int(percentage)+" "*(100-int(percentage))+" " + str(percentage)+" %\n", end="\r")
        out.close()
        print("Writing Casette Exon coordinates into .bed file ... ", end="\r")
        "Print .bed file, can be removed later."
        out=open(args.out+"CE.bed", "w")
        for gene in gene_exons:
            strand=gene_exons[gene][0][3]
            
            for entry in gene_exons[gene]:
                chrom, start, stop, strand = entry[0:4]
                name="CE_"+chrom+"_"+start+"_"+stop
                out.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, start, stop, name, strand))
        out.close()
        print("Writing Casette Exon coordinates into .bed file: Done! \n", end="\r")
    
    #"Alternative Donors/Acceptors"
    elif event.upper()=="AD" or event.upper()=="AA" or event.upper()=="IR":
        if AA_AD_dict ==False:
            print("\nFinding Alternative Donors/Acceptors: \n", end="\r")
            #Creates dictionary with exons per gene.
            gene_exons=dict()
            
            #exon count
            total_exon_count=0
            #Remove transcript ids and duplicate exons.
            for gene in gene_dict:
                if gene_dict[gene]:
                    gene_exons[gene]=[]
                    for trans_ID in gene_dict[gene]:
                        for exon in gene_dict[gene][trans_ID]:
                            if exon[0:-1] not in gene_exons[gene]:
                                gene_exons[gene].append(exon[0:-1])
                                total_exon_count+=1
            
            #Go through all exons per gene to find alternative donors/acceptors.
            potential_AA=dict()
            potential_AD=dict()
            #key=start, value=stop
            coordinates=dict()
            #key=start_stop, value=entry.
            coord_exons=dict()
            #initiate progress counter
            c=0
            #For progress updates during PSI.
            AA_counter=0
            AD_counter=0
            for gene in gene_exons:
                #Initiate potential AA/AD for each gene.
                potential_AA[gene]=[]
                potential_AD[gene]=[]
                #Go through exons to find potentials
                for entry in gene_exons[gene]:
                    start=int(entry[1])
                    stop=int(entry[2])
                    strand=entry[3]
                    
                    #print(start, stop)
                    
                    """If the start or stop is shared with another exon, then we found potential AS.
                    But one of the exons can also be completely within another exon, or completely
                    surrounding the other one. Or having one coordinate within the other. That also counts"""
                    
                    for coordinate in coordinates:
                        key_string=str(coordinate)+"_"+str(coordinates[coordinate])
                        #if the start coordinate lies within a previous exon.
                        #print(coordinates[coordinate], stop)
                        if coordinate<=start<coordinates[coordinate] and coordinates[coordinate]!=stop:
                            #Then either they share start, and we only have one event
                            if start==coordinate:
                                # if + AD, if - AA.
                                if strand=="+":
                                    AA_counter, AD_counter=add_to(["AD"], AA_counter, AD_counter)
                                else:
                                    #print(entry, coord_exons[key_string])
                                    AA_counter, AD_counter=add_to(["AA"], AA_counter, AD_counter)
                            #or we have both events
                            else:
                                #both
                                AA_counter, AD_counter=add_to(["AA","AD"], AA_counter, AD_counter)
                                
                        #if the stop coordinate lies within a previous exon.
                        elif coordinate+1< stop <= coordinates[coordinate] and start!= coordinate:
                            #either they share stop, and we have one event
                            if stop==coordinates[coordinate]:
                                #if - AD, if + AA.
                                if strand=="-":
                                    AA_counter, AD_counter=add_to(["AD"], AA_counter, AD_counter)
                                else:
                                    AA_counter, AD_counter=add_to(["AA"], AA_counter, AD_counter)
                            #Or we have both events.
                            else:
                                #both
                                AA_counter, AD_counter=add_to(["AA","AD"], AA_counter, AD_counter)
                        # if this exon contains a previous exon completely then we have both
                        elif coordinate > start and coordinates[coordinate] < stop:
                            #both
                            AA_counter, AD_counter=add_to(["AA","AD"], AA_counter, AD_counter)
                    
                    #Add new coordinates to compare new entries to.
                    coord_exons[str(start)+"_"+str(stop)]=entry
                    coordinates[start]=stop
                    # Progress bar
                    c+=1
                    percentage= round(c/total_exon_count * 100, 2)
                    print("#"*int(percentage)+ " "*(100-int(percentage))+ " "+ str(percentage)+ " %", end="\r")
            #Final bar
            print("#"*int(percentage)+ " "*(100-int(percentage))+ " "+ str(percentage)+ " %\n", end="\r")
            AA_AD_dict=True

        # Calculate psi scores for potential AA, if asked for.
        if event=="AA":
            print("Calculating PSI scores for alternative acceptors: \n", end="\r")
            out=open(args.out+"PSI_AA.tsv", "w")
            out.write("#PSI Table AA \n")
            out.write("Event_Location\t"+"\t".join(sample_names)+"\n")
            #Progress
            counter=0
            for gene in potential_AA:
                #If there is exons for this gene:
                if len(potential_AA[gene])!=0:
                    #Extrat strand information for header (for each gene)
                    strand=gene_ranges[gene][1]
                
                    #Create gene header
                    out.write("#"+gene+", "+strand+"\n")
                    #sort events
                    for event in potential_AA[gene]:
                        event.sort()
                    potential_AA[gene].sort(key=lambda x: int(x[0]))
                    
                    #numerate events
                    event_ID=0
                    for event in potential_AA[gene]:
                        counter+=1
                        event_ID+=1
                        for starts in event:
                            PSI_scores=[]
                            for sample in sample_names:
                                PSI_scores.append(PSI_AA(gene, sample, event, starts))
                    
                            out.write(str(event_ID)+"_"+str(starts)+"\t"+ "\t".join(PSI_scores)+"\n")
                        #Progress update
                        percentage= round(counter/AA_counter *100, 2)
                        print("#"*int(percentage)+ " "*(100-int(percentage))+ " "+ str(percentage)+ " %", end="\r")
            
            if AA_counter!=0:
                print("#"*int(percentage)+ " "*(100-int(percentage))+ " "+ str(percentage)+ " %\n", end="\r")
            else:
                print("No Alternative Acceptor events found! \n")
            out.close()
            
            print("Printing alternative acceptor coordinates into .bed file...", end="\r")
            "Print .bed file, can be removed later."
            out=open(args.out+"AA.bed", "w")
            for gene in potential_AA:
                if len(potential_AA[gene])!=0:
                    #Extrat strand information for header (for each gene)
                    strand=gene_ranges[gene][1]
                    chrom=gene_ranges[gene][0]
                
                for entry in potential_AA[gene]:
                    smallest_start=min(entry)
                    biggest_start=max(entry)
                    name="AA_"+chrom+"_"+str(smallest_start)+"_"+str(biggest_start)
                    out.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, str(smallest_start), str(biggest_start), name, strand))
            out.close()
            print("Printing alternative acceptor coordinates into .bed file: Done! \n", end="\r")

        elif event=="AD":
            print("Calculating PSI scores for Alternative Donors: \n", end="\r")
            out=open(args.out+"PSI_AD.tsv", "w")
            out.write("#PSI Table AD \n")
            out.write("Event_Location"+"\t".join(sample_names)+"\n")
            #Progress
            counter=0
            for gene in potential_AD:
                #If there is exons for this gene:
                if len(potential_AD[gene])!=0:
                    #Extrat strand information for header (for each gene)
                    strand=gene_ranges[gene][1]
                
                    #Create gene header
                    out.write("#"+gene+", "+strand+"\n")
                    #sort events
                    for event in potential_AD[gene]:
                        event.sort()
                    potential_AD[gene].sort(key=lambda x: int(x[0]))
                    
                    #numerate events
                    event_ID=0
                    for event in potential_AD[gene]:
                        counter+=1
                        event_ID+=1
                        for stops in event:
                            PSI_scores=[]
                            for sample in sample_names:
                                PSI_scores.append(PSI_AD(gene, sample, event, stops))
                    
                            out.write(str(event_ID)+"_"+str(stops)+"\t"+ "\t".join(PSI_scores)+"\n")
                        #Progress update
                        percentage= round(counter/AD_counter *100, 2)
                        print("#"*int(percentage)+ " "*(100-int(percentage))+ " "+ str(percentage)+ " %", end="\r")
            
            if AD_counter!=0:
                print("#"*int(percentage)+ " "*(100-int(percentage))+ " "+ str(percentage)+ " %\n", end="\r")  
            else:
                print("No Alternative Donor events found! \n")
            out.close()

            print("Printing alternative donor coordinates into .bed file ...", end="\r")
            "Print .bed file, can be removed later."
            out=open(args.out+"AD.bed", "w")
            for gene in potential_AD:
                if len(potential_AD[gene])!=0:
                    #Extrat strand information for header (for each gene)
                    strand=gene_ranges[gene][1]
                    chrom=gene_ranges[gene][0]
                
                for entry in potential_AD[gene]:
                    smallest_stop=min(entry)
                    biggest_stop=max(entry)
                    name="AD_"+chrom+"_"+str(smallest_stop)+"_"+str(biggest_stop)
                    out.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, str(smallest_stop), str(biggest_stop), name, strand))
            out.close()
            print("Printing alternative donor coordinates into .bed file: Done! \n", end="\r")

        else:
            #Intron Retention, because invalid inputs are checked at the beginning of the script to save time
            #Retention of introns is usually not annotated. So we do the same as for CE.
            #And check for every gap between exons, if there is read showing intron retention.
            #That means at this point we save all potential gaps.
            #Format: chr_start1_stop1_start2_stop2
            print("\nFinding potential intron retention events... ", end="\r")
            IR_coord=dict()
            #IR counter for PSI progress
            IR_counter=0
            #Remove transcript ids and duplicate exons.
            for gene in gene_dict:
                if gene_dict[gene]:
                    IR_coord[gene]=[]
                    for trans_ID in gene_dict[gene]:
                        for i in range(0, len(gene_dict[gene][trans_ID])-1):
                            chrom=gene_dict[gene][trans_ID][i][0]
                            strand=gene_dict[gene][trans_ID][i][3]
                            IR=strand+"_"+chrom+"_"+gene_dict[gene][trans_ID][i][2]+"_"+gene_dict[gene][trans_ID][i+1][1]
                            if IR not in IR_coord[gene]:                        
                                IR_coord[gene].append(IR)
                                IR_counter+=1
            print("Finding potential intron retention events: Done! \n", end="\r")
            
            out=open(args.out+"PSI_IR.tsv", "w")
            print("Calculating PSI scores for intron retention events \n", end="\r")
            out.write("#PSI Table IR \n")
            out.write("Location\tPosition_Transcript\t"+"\t".join(sample_names)+"\n")
            #progress counter
            counter=0
            for gene in IR_coord:
                #If there is exons for this gene:
                if len(IR_coord[gene])!=0:
                    #Extrat strand information for header (for each gene)
                    strand=IR_coord[gene][0].split("_")[0]
                
                    #Create gene header
                    out.write("#"+gene+", "+strand+"\n")
                    #sort by first exons start.
                    IR_coord[gene].sort(key=lambda x: int(x.split("_")[2]))
                for exon in IR_coord[gene]:
                    PSI_scores=[]
                    for sample in sample_names:
                        PSI_scores.append(PSI_IR(sample, exon, gene))
                
                    out.write("_".join(exon.split("_")[1::])+"\t"+ "\t".join(PSI_scores)+"\n")
                    counter+=1
                    percentage=round(counter/IR_counter *100, 2)
                    print("#"*int(percentage)+ " "*(100-int(percentage))+ " "+ str(percentage)+ " %", end="\r")
            
            if IR_counter!=0:
                print("#"*int(percentage)+ " "*(100-int(percentage))+ " "+ str(percentage)+ " %\n", end="\r")
            else:
                print("No potential IR events found! \n", end="\r")
            out.close()
            "Print .bed file, can be removed later."
            print("Printing intron retention coordinates into .bed file ... ", end="\r")
            out=open(args.out+"IR.bed", "w")
            for gene in IR_coord:
                strand=IR_coord[gene][0].split("_")[0]
                
                for entry in IR_coord[gene]:
                    chrom = entry.split("_")[1]
                    start=entry.split("_")[2]
                    stop=entry.split("_")[3]
                    name="IR_"+chrom+"_"+start+"_"+stop
                    out.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, start, stop, name, strand))
            out.close()
            print("Printing intron retention coordinates into .bed file: Done! \n", end="\r")

#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  

















