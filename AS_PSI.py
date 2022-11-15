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
        python AS_PSI.py -s ../Sample_Data/ -o PSI_scores/ -c "chr6:151690496-152103274" -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as CE
        
        #With coordinates f.e. SYNE1 (neg strand)
        python AS_PSI.py -s ../Sample_Data/ -o PSI_scores/ -c "chr6:152121515-152652710" -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as ALL
        
        
        #for server
        python AS_PSI.py -s /raidset/mi3258mu-se/Mirjam -o PSI_ESR1.tsv -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as CE -c "chr6:151690496-152103274"
        

        #no coordinates/name
        
Possible Bugs:
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

#%% User Defined Functions


def PSI_CE(sample, exon):
    """
    

    Parameters
    ----------
    sample : string
        Sample name
    exon : list
        contains chromosome, start, stop, strand information.

    Returns
    -------
    str
        PSI score

    """
    
    return "NAN"

def PSI_AA(sample, exon):
    """
    

    Parameters
    ----------
    sample : string
        Sample name
    exon : list
        contains chromosome, start, stop, strand information.

    Returns
    -------
    str
        PSI score

    """
    
    return "NAN"

def PSI_AD(sample, exon):
    """
    

    Parameters
    ----------
    sample : string
        Sample name
    exon : list
        contains chromosome, start, stop, strand information.

    Returns
    -------
    str
        PSI score

    """
    
    return "NAN"

def PSI_IR(sample, exon):
    """
    

    Parameters
    ----------
    sample : string
        Sample name
    exon : list
        contains chromosome, start, stop, strand information.

    Returns
    -------
    str
        PSI score

    """
    
    return "NAN"

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
                    if i==0:
                        if strand=="+":
                            position= "first"
                        else:
                            position="last"
                    elif i==number_exons-2:
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

#%% SAMPLE NAMES?

#All bam files.
argument_glob=args.samples+"/**/alignment.bam"
bam_file_list=glob.glob(argument_glob, recursive=True)

#Assuming that each bam.file is in its own sample folder:
sample_names=[i.split("/")[-4] for i in bam_file_list]


#%% Identify alternative splicing events.


for event in inputs:
    
    "Casette Exons"
    
    if event.upper()=="CE":
        out=open(args.out+"PSI_CE.tsv", "w")
        #if a transcript has 2 or less exons, there can be no CE.
        filtered_gene_dict=dict()
        for gene in gene_dict:
            filtered_gene_dict[gene]={k : v for k , v in gene_dict[gene].items() 
                                if len(v) >2}
        
        #Remove first and last exon for each transcript, as they cannot be CE
        new_gene_dict=dict()
        for gene in filtered_gene_dict:
            if len(filtered_gene_dict[gene].items())>0:
                new_gene_dict[gene]=dict()
                for trans_ID in filtered_gene_dict[gene]:
                    new_gene_dict[gene][trans_ID]=[i for i in filtered_gene_dict[gene][trans_ID] if i[4]=="middle"]  
                
        CE_gene_dict=new_gene_dict
        #print(CE_gene_dict)
        """Transcript IDs are no longer required, and additionally there
        is definitely duplicate exons between R and G as well as within each
        of them. Those need to be removed."""
        
        #Creates dictionary with potential casette exons per gene.
        gene_exons=dict()
        #Creates a dictionary where each exon is listed with its database (in case needed later)
        exons_db=dict()
        
        for gene in CE_gene_dict:
            if CE_gene_dict[gene]:
                gene_exons[gene]=[]
                for trans_ID in CE_gene_dict[gene]:
                    for exon in CE_gene_dict[gene][trans_ID]:
                        if exon[0:-1] not in gene_exons[gene]:
                            gene_exons[gene].append(exon[0:-1])
                        if "_".join(exon[0:-1]) in exons_db:
                            if exon[-1] in exons_db["_".join(exon[0:-1])]:
                                continue
                            else:
                                exons_db["_".join(exon[0:-1])]+=exon[-1]
                        else:
                            exons_db["_".join(exon[0:-1])]=exon[-1]

        #Make header for output file
        out.write("#PSI Table CE \n")
        out.write("Location\t"+"\t".join(sample_names)+"\n")
        #Go through all potential casette exons, and calculate their PSI scores.
        for gene in gene_exons:
            #Extract strand information for header (for each gene)
            strand=gene_exons[gene][0][3]
            
            #Create gene header
            out.write("#"+gene+", "+strand+"\n")
            for entry in gene_exons[gene]:
                PSI_scores=[]
                chrom=exon[0]
                start=int(exon[1]) #Note that these are start and stop coor
                stop=int(exon[2]) #not start and stop of exon. Strand dependent.
                strand=exon[3]
                
                for sample in sample_names:
                    PSI_scores.append(PSI_CE(sample, entry))
                
                out.write("_".join(entry[0:3])+"\t" +"\t".join(PSI_scores)+"\n")
        
        out.close()
    
    #"Alternative Donors/Acceptors"
    elif event.upper()=="AD" or event.upper()=="AA":
        #Creates dictionary with exons per gene.
        gene_exons=dict()
        
        #Remove transcript ids and duplicate exons.
        for gene in gene_dict:
            if gene_dict[gene]:
                gene_exons[gene]=[]
                for trans_ID in gene_dict[gene]:
                    for exon in gene_dict[gene][trans_ID]:
                        if exon[0:-1] not in gene_exons[gene]:
                            gene_exons[gene].append(exon[0:-1])
        
        #Go through all exons per gene to find alternative donors/acceptors.
        potential_AA=dict()
        potential_AD=dict()
        #key=start, value=stop
        coordinates=dict()
        #key=start_stop, value=entry.
        coord_exons=dict()
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
                    #If the start= start of another exon or within range of another exon its an event.
                    if start >= coordinate and start <=coordinates[coordinate]:
                        #if start == start of another exon, its only AD (on plus). otherwise its both.
                        if start==coordinate:
                            if strand=="+" and entry[4]!="last" and coord_exons[key_string][4]!="last":
                                if entry not in potential_AD[gene]:
                                    potential_AD[gene].append(entry)
                                if coord_exons[key_string] not in potential_AD[gene]:
                                    potential_AD[gene].append(coord_exons[key_string])
                            elif strand=="-" and entry[4]!="first":
                                if entry not in potential_AA[gene]:
                                    potential_AA[gene].append(entry)
                                if coord_exons[key_string] not in potential_AA[gene]:
                                    potential_AA[gene].append(coord_exons[key_string])
                        else:
                            if strand=="+":
                                if entry[4]!="last" and coord_exons[key_string]!="last":
                                    if entry not in potential_AD[gene]:
                                        potential_AD[gene].append(entry)
                                    if coord_exons[key_string] not in potential_AD[gene]:
                                        potential_AD[gene].append(coord_exons[key_string])
                            if strand=="-":
                                if entry[4]!="first" and coord_exons[key_string]!="first":
                                    if entry not in potential_AA[gene]:
                                        potential_AA[gene].append(entry)
                                    if coord_exons[key_string] not in potential_AA[gene]:
                                        potential_AA[gene].append(coord_exons[key_string])
                                    
                            
                    #If the stop lays anywhere within the previous exon, its an event too.
                    elif stop >= coordinate and stop <=coordinates[coordinate]:
                        #If stop==stop of another exon it is only AA (+). otherwise its both.
                        if stop==coordinates[coordinate]:
                            if strand=="+" and entry[4]!="first" and coord_exons[key_string][4!="first"]:
                                if entry not in potential_AA[gene]:
                                    potential_AA[gene].append(entry)
                                if coord_exons[key_string] not in potential_AA[gene]:
                                    potential_AA[gene].append(coord_exons[key_string])
                            elif strand=="-" and entry[4]!="last" and coord_exons[key_string][4]!="last":
                                if entry not in potential_AD[gene]:
                                    potential_AD[gene].append(entry)
                                if coord_exons[key_string] not in potential_AD[gene]:
                                    potential_AD[gene].append(coord_exons[key_string])
                        else:
                            if strand=="-":
                                if entry[4]!="last" and coord_exons[key_string]!="last":
                                    if entry not in potential_AD[gene]:
                                        potential_AD[gene].append(entry)
                                    if coord_exons[key_string] not in potential_AD[gene]:
                                        potential_AD[gene].append(coord_exons[key_string])
                            if strand=="+":
                                if entry[4]!="first" and coord_exons[key_string]!="first":
                                    if entry not in potential_AA[gene]:
                                        potential_AA[gene].append(entry)
                                    if coord_exons[key_string] not in potential_AA[gene]:
                                        potential_AA[gene].append(coord_exons[key_string])
                                        
                    #If the new exon contains a previous exon, its also an event.
                    elif stop>coordinates[coordinate] and start<coordinate:
                        if strand=="+":
                            if entry[4]!="first" and coord_exons[key_string][4]!="first":
                                if entry not in potential_AA[gene]:
                                    potential_AA[gene].append(entry)
                                if coord_exons[key_string] not in potential_AA[gene]:
                                    potential_AA[gene].append(coord_exons[key_string])
                        
                        if strand=="-":
                            if entry[4]!="last" and coord_exons[key_string][4]!="last":
                                if entry not in potential_AD[gene]:
                                    potential_AD[gene].append(entry)
                                if coord_exons[key_string] not in potential_AD[gene]:
                                    potential_AD[gene].append(coord_exons[key_string])
                    #Note that the new exon being contained in a previous one does not need its own statement, as it will be caught by statement 1.
                
                #Add new coordinates to compare new entries to.
                coord_exons[str(start)+"_"+str(stop)]=entry
                coordinates[start]=stop
        """
        print("AA-events")
        for gene in potential_AA:
            for i in potential_AA[gene]:
                print(i[0]+":"+i[1]+"-"+i[2])
        
        print("AD-events")
        for gene in potential_AD:
            for i in potential_AD[gene]:
                print(i[0]+":"+i[1]+"-"+i[2])
        
        """
        #print(potential_AA)
        #print(potential_AD)
        # Calculate psi scores for potential AA, if asked for.
        if event=="AA":
            out=open(args.out+"PSI_AA.tsv", "w")
            out.write("#PSI Table AA \n")
            out.write("Location\tPosition_Transcript\t"+"\t".join(sample_names)+"\n")
            for gene in potential_AA:
                #If there is exons for this gene:
                if len(potential_AA[gene])!=0:
                    #Extrat strand information for header (for each gene)
                    strand=potential_AA[gene][0][3]
                
                    #Create gene header
                    out.write("#"+gene+", "+strand+"\n")
                for exon in potential_AA[gene]:
                    PSI_scores=[]
                    for sample in sample_names:
                        PSI_scores.append(PSI_AA(sample, exon))
                
                    out.write("_".join(exon[0:3])+"\t"+exon[4]+"\t"+ "\t".join(PSI_scores)+"\n")
            out.close()
        else:
            out=open(args.out+"PSI_AD.tsv", "w")
            out.write("#PSI Table AD \n")
            out.write("Location\tPosition_Transcript\t"+"\t".join(sample_names)+"\n")
            for gene in potential_AD:
                #If there is exons for this gene:
                if len(potential_AD[gene])!=0:
                    #Extrat strand information for header (for each gene)
                    strand=potential_AD[gene][0][3]
                
                    #Create gene header
                    out.write("#"+gene+", "+strand+"\n")
                for exon in potential_AD[gene]:
                    PSI_scores=[]
                    for sample in sample_names:
                        PSI_scores.append(PSI_AD(sample, exon))
                
                    out.write("_".join(exon[0:3])+"\t"+exon[4]+"\t"+ "\t".join(PSI_scores)+"\n")
            out.close()

    else:
        #Intron Retention, because invalid inputs are checked at the beginning of the script to save time
        #Retention of introns is usually not annotated. So we do the same as for CE.
        #And check for every gap between exons, if there is read showing intron retention.
        #That means at this point we save all potential gaps.
        #Format: chr_start1_stop1_start2_stop2
        
        IR_coord=dict()
        #Remove transcript ids and duplicate exons.
        for gene in gene_dict:
            if gene_dict[gene]:
                IR_coord[gene]=[]
                for trans_ID in gene_dict[gene]:
                    for i in range(0, len(gene_dict[gene][trans_ID])-1):
                        chrom=gene_dict[gene][trans_ID][i][0]
                        strand=gene_dict[gene][trans_ID][i][3]
                        IR=strand+"_"+chrom+"_"+"_".join(gene_dict[gene][trans_ID][i][1:3])+"_"+"_".join(gene_dict[gene][trans_ID][i+1][1:3])
                        if IR not in IR_coord[gene]:                        
                            IR_coord[gene].append(IR)
        """        
        #Find IR events.
        IR_coord=dict()
        #key=start, value=stop
        coordinates=dict()
        for gene in gene_exons:
            IR_coord[gene]=[]
            #Go through exons to find potentials
            for entry in gene_exons[gene]:
                start=int(entry[1])
                stop=int(entry[2])
                strand=entry[3]
                
                #Flags
                Start_Found=False
                Stop_Found=False
                exon1=""
                exon2=""
                #go through previously saved coordinates. Check if start and stop lay in two different exons.
                for starts, stops in coordinates:
                    key_string=str(starts)+"_"+str(stops)
                    if Start_Found==True and Stop_Found==True:
                        if exon1!=exon2 and exon1!="" and exon2!="":
                            IR_coord[gene].append(entry)
                            Start_Found=False
                            Stop_Found=False
                    if start >=starts and start <=stops:
                        Start_Found=True
                        exon1= key_string
                        continue
                    elif stop>=starts and stop <=stops:
                        Stop_Found=True
                        exon2=key_string
                        continue
                    
        """
        
        out=open(args.out+"PSI_IR.tsv", "w")
        out.write("#PSI Table IR \n")
        out.write("Location\tPosition_Transcript\t"+"\t".join(sample_names)+"\n")
        for gene in IR_coord:
            #If there is exons for this gene:
            if len(IR_coord[gene])!=0:
                #Extrat strand information for header (for each gene)
                strand=IR_coord[gene][0].split("_")[0]
            
                #Create gene header
                out.write("#"+gene+", "+strand+"\n")
            for exon in IR_coord[gene]:
                PSI_scores=[]
                for sample in sample_names:
                    PSI_scores.append(PSI_IR(sample, exon))
            
                out.write("_".join(exon.split("_")[1::])+"\t"+ "\t".join(PSI_scores)+"\n")
        out.close()



















