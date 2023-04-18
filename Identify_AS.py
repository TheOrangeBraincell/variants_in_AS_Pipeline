# -*- coding: utf-8 -*-
"""
Date: Wed Jan 18 09:22:39 2023
File Name: Identify_AS.py
Author: Mirjam Karlsson-Müller

Description:
    Finds AS events in annotated genes from GENCODE and RefSeq.
    Returns Table sorted by gene, containing all AS event locations and types.
    Also returns .bed files if wished for, one per type event, with event location entries, possible to view in f.e. igv.

Abbreviations:
    AS=     Alternative Splicing
    CE=     Casette Exon
    AA=     Alternative Acceptor
    AD=     Alternative Donor
    IR=     Intron Retention
    
List of Functions:
    add_to(): 
        Adds exon to AD and/or AA event. If the event this exon belongs to does not exist yet, creates new one.
        Specific cases for plus and negative strand.
    
Procedure: 
    1. Creates a dictionary with genes, transcripts and their annotated exons from databases GENCODE and RefSeq.
    2. Going through this dictionary, finds potential AS events:
        AD/AA: Find exons with overlap, put them into AA, AD dictionary according to type of overlap and strand.
        CE: Every exon thats not first or last is a potential CE.
        IR: Every space between two exons in the same transcript can be a potential IR event.

    
Useage:
    Inputs:
        - Database files for GENCODE and RefSeq respectively
        - output file name for table of eventlocations/types
        - opt. coordinates of a gene/region of interest
        - What type of alternative splicing event are we identifying? Multiple possible in listformat. Or ALL for all AS types.
        - If IR is one of the AS events, requires mean and sd of intron size.
        - Whether bed file is wished for
    
    Outputs:
        - table of event locations and types, as well as gene information and transcript information.
        - (optional) a .bed file containing the location of AS events.


    Instructions:
        Run in command line. For example.
        
        #with coordinates f.e. Estrogen Receptor
        python variants_in_AS_Pipeline/Identify_AS.py -o AS_events_ESR1.tsv -c "chr6:151656691-152129619" -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as ALL 
        
        #With coordinates f.e. BRCA1 (neg strand)
        python variants_in_AS_Pipeline/Identify_AS.py -o AS_events_BRCA1.tsv -c "chr17:43044295-43170245" -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as ALL 
        
       
    
Possible Bugs:
    -not compatible with other database tables, as specifically tailored to the two used.
    Instructions on how to retrieve the right ones from the ucsc table browser can be found
    on the github.
"""


#%% Imports

import argparse
import re
import time



#%% Time

start_time=time.time()

#%% 0. argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Identify Alternative Splicing',
                                 usage='%(prog)s -o OUTPUT-FILE \
                                     -g GENCODE-FILE -r REFSEQ-FILE \
                                         [-c] "chrX:XXXXXX-XXXXXX" -as AS-TYPE',
                                 description="""Per AS event of interest, creates
                                 a table with the PSI scores supporting said event
                                 per sample in sample folder.""")

parser.add_argument('--out', '-o', required=True,
                    help="""Output file, containing events and type per gene.""")
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")
parser.add_argument('--gencode', '-g', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from GENCODE39 as well as gene names.""")
parser.add_argument('--refseq', '-r', required=True,
                    help="""tsv file containing bed file information on 
                    annotated exons from RefSeq as well as gene names.""")
parser.add_argument('--AS', '-as', type=str, required=True,
                    help="""Which type of alternative splicing event we are
                    interested in. "CE" for Casette Exons, "AA" for alternative
                    acceptors, "AD" for alternative donors, "IR" for intron
                    retention and "ALL" for all of the types. Several seperated
                    by ,.""")
parser.add_argument('--bed', '-b', type=bool,
                    help="If an output file bed file for each type of event is wished for,\
                        set to True")


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
    inputs=[i.upper() for i in inputs]
    for i in inputs:
        if i not in allowed_inputs:
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
    #if IR is in the events, but CE is not, then CE still has to run! Because needed for PSI of IR
    if "IR" in inputs and "CE" not in inputs:
        print("Casette exons (CE) will also be identified as needed for PSI score calculation of intron retention (IR).")
        inputs.append("CE")
        
#Sort inputs so that CE will be scored before IR.
inputs=sorted(inputs)



#%% User Defined Functions

def add_to(events):
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
            for events in potential_AA:
                if start1 in events[0] and start2 not in events[0]:
                    events[0].append(start2)
                    event_exists=True
                elif start2 in events and start1 not in events:
                    events[0].append(start1)
                    event_exists=True
            
            #If the event is not found, make a new one.
            if event_exists==False:
                potential_AA.append([[int(start1), int(start2)], entry[-2], entry[-1]])
        
            
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
            for events in potential_AD:
                if stop1 in events[0] and stop2 not in events[0]:
                    events[0].append(stop2)
                    event_exists=True
                elif stop2 in events[0] and stop1 not in events[0]:
                    events[0].append(stop1)
                    event_exists=True
                elif stop1 in events[0] and stop2 in events[0]:
                    event_exists=True
                    
            #If the event is not found, make a new one.
            if event_exists==False:
                potential_AD.append([[stop1, stop2], entry[-2], entry[-1]])

        


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

    
#%% Identify alternative splicing events

out=open(args.out,"w")
out.write("Location\tAS_Type\tGene_Name\tDatabase\tGencode_Transcript_IDs\tRefseq_Transcript_IDs\n")

#initiate progress counter
total_count=len(list(gene_ranges.keys()))*len(inputs)
current_count=0
print("Finding AS events: {:.2f}%".format(0), end="\r")

#So that the output is sorted by gene not by event
for gene in gene_dict:
    #gene header for output file
    out.write("#"+gene+" ,"+str(gene_ranges[gene][0])+" ,"+str(gene_ranges[gene][1])+
              " ,"+str(gene_ranges[gene][2])+" ,"+str(gene_ranges[gene][3])+"\n")
    #Go through events
    for event in inputs:
        "Alternative Donors/Acceptors"
        if event.upper()=="AD" or event.upper()=="AA":
            
            #Remove transcript ids and duplicate exons.
            if gene_dict[gene]:
                #list of exons for current gene
                gene_exons=[]
                for trans_ID in gene_dict[gene]:
                    for exon in gene_dict[gene][trans_ID]:

                        #exon is now: [chrom, small, big, strand, position, db]
                        #Because db and transid info wont match necessarily, we only check for up to -2
                        no_match=True
                        for index,entry in enumerate(gene_exons):
                            if exon[0:-1] == entry[0:-2]:
                                no_match=False
                                #if its already there, update db and Transcript ID

                                #update db according to whether this one is included already.
                                #The if statement will catch it only if its the opposite database by itself at this point in the entry.
                                if exon[-1] not in entry[-2]:
                                    entry[-2]=entry[-2]+exon[-1]
                                    #Then the dictionary at the end of entry will also have - for the new database.
                                    entry[-1][exon[-1]]=[trans_ID]
                                else:
                                    #there should already be a transcript id from that database in the dictionary at the last spot. append.
                                    entry[-1][exon[-1]].append(trans_ID)
                                #Update gene_exons entry
                                gene_exons[index]=entry
                        #If there was no match found, then we have a completely new exon and append it to gene_exons
                        if no_match==True:
                            if exon[-1]=="G":
                                temp_list=[i for i in exon]
                                temp_list.append({exon[-1]:[trans_ID], "R":"-"})
                                gene_exons.append(temp_list)
                            elif exon[-1]=="R":
                                temp_list=[i for i in exon]
                                temp_list.append({exon[-1]:[trans_ID], "G":"-"})
                                gene_exons.append(temp_list)
                            else:
                                print("WTF? ",  exon)

            #Go through all exons per gene to find alternative donors/acceptors.
            #Initiate potential AA/AD for each gene.
            potential_AA=[]
            potential_AD=[]
            #Dont need to save coordinates from different gene.
            #key=start, value=stop
            coordinates=dict()
            #key=start_stop, value=entry.
            coord_exons=dict()
            #Go through exons to find potentials
            for entry in gene_exons:
                start=int(entry[1])
                stop=int(entry[2])
                strand=entry[3]
                
                """If the start or stop is shared with another exon, then we found potential AS.
                But one of the exons can also be completely within another exon, or completely
                surrounding the other one. Or having one coordinate within the other. That also counts. So we find any type of overlap"""
                
                for coordinate in coordinates:
                    key_string=str(coordinate)+"_"+str(coordinates[coordinate])
                    #if the start coordinate lies within a previous exon.
                    #print(coordinates[coordinate], stop)
                    if coordinate<=start<coordinates[coordinate] and coordinates[coordinate]!=stop:
                        #Then either they share start, and we only have one event
                        if start==coordinate:
                            # if + AD, if - AA.
                            if strand=="+":
                                add_to(["AD"])
                            else:
                                #print(entry, coord_exons[key_string])
                                add_to(["AA"])
                        #or we have both events
                        else:
                            #both
                            add_to(["AA","AD"])
                            
                    #if the stop coordinate lies within a previous exon.
                    elif coordinate+1< stop <= coordinates[coordinate] and start!= coordinate:
                        #either they share stop, and we have one event
                        if stop==coordinates[coordinate]:
                            #if - AD, if + AA.
                            if strand=="-":
                                add_to(["AD"])
                            else:
                                add_to(["AA"])
                        #Or we have both events.
                        else:
                            #both
                            add_to(["AA","AD"])
                    # if this exon contains a previous exon completely then we have both
                    elif coordinate > start and coordinates[coordinate] < stop:
                        #both
                        add_to(["AA","AD"])
                
                #Add new coordinates to compare new entries to.
                coord_exons[str(start)+"_"+str(stop)]=entry
                coordinates[start]=stop
                        
                    
            #Write all potential AA/AD into output file for each gene, then reset dictionaries
            if event=="AA":
                #Sort entries after location and coordinates
                for aa_event in potential_AA:
                    aa_event[0].sort()
                potential_AA.sort(key=lambda x: int(x[0][0]))
                
                #numerate events
                event_ID=0
                for aa_event in potential_AA:
                    event_ID+=1
                    for starts in aa_event[0]:
                        out.write("{}_{}_{}_{}\t{}\t{}\t{}\t{}\t{}\n".format(event_ID, gene_ranges[gene][0],
                                                             gene_ranges[gene][1],
                                                             starts, "AA", gene, aa_event[1], ",".join(aa_event[2]["G"]), ",".join(aa_event[2]["R"])))
                if args.bed==True:
                    bed=open(args.out+"AA.bed", "w")
                    if len(potential_AA)!=0:
                        #Extract strand information for header (for each gene)
                        strand=gene_ranges[gene][1]
                        chrom=gene_ranges[gene][0]
                    
                    for entry in potential_AA:
                        smallest_start=min(entry)
                        biggest_start=max(entry)
                        name="AA_"+chrom+"_"+str(smallest_start)+"_"+str(biggest_start)
                        bed.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, str(smallest_start), str(biggest_start), name, strand))
                    bed.close()
            #Same for AD
            elif event=="AD":
                #Sort entries after location and coordinates
                for ad_event in potential_AD:
                    ad_event[0].sort()
                potential_AD.sort(key=lambda x: int(x[0][0]))
                
                #numerate events
                event_ID=0
                for ad_event in potential_AD:
                    event_ID+=1
                    #print(ad_event)
                    for stops in ad_event[0]:
                        out.write("{}_{}_{}_{}\t{}\t{}\t{}\t{}\t{}\n".format(event_ID, gene_ranges[gene][0],
                                                             gene_ranges[gene][1],
                                                             stops, "AD", gene, ad_event[1], ",".join(ad_event[2]["G"]), ",".join(ad_event[2]["R"])))
                if args.bed==True:
                    bed=open(args.out+"AD.bed", "w")
                    if len(potential_AD)!=0:
                        #Extract strand information for header (for each gene)
                        strand=gene_ranges[gene][1]
                        chrom=gene_ranges[gene][0]
                    
                    for entry in potential_AD:
                        smallest_stop=min(entry)
                        biggest_stop=max(entry)
                        name="AA_"+chrom+"_"+str(smallest_stop)+"_"+str(biggest_stop)
                        bed.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, str(smallest_stop), str(biggest_stop), name, strand))
                    bed.close()  
            
        
        #"Casette Exons"
        elif event.upper()=="CE":
            #Remove first and last exon for each transcript, as they cannot be CE
            CE_exons=[]
            for trans_ID in gene_dict[gene]:
                #If a transcript has 2 or less exons, there is no CE
                if len(gene_dict[gene][trans_ID])<=2:
                    continue
                for exon in gene_dict[gene][trans_ID]:
                    if exon[4]=="middle":
                        #Check if they are already in CE exons
                        no_match=True
                        for index,entry in enumerate(CE_exons):
                            if exon[0:-1] == entry[0:-2]:
                                no_match=False
                                #if its already there, update db and Transcript ID

                                #update db according to whether this one is included already.
                                #The if statement will catch it only if its the opposite database by itself at this point in the entry.
                                if exon[-1] not in entry[-2]:
                                    entry[-2]=entry[-2]+exon[-1]
                                    #Then the dictionary at the end of entry will also have - for the new database.
                                    entry[-1][exon[-1]]=[trans_ID]
                                else:
                                    #there should already be a transcript id from that database in the dictionary at the last spot. append.
                                    entry[-1][exon[-1]].append(trans_ID)
                                #Update gene_exons entry
                                CE_exons[index]=entry
                        #If there was no match found, then we have a completely new exon and append it to gene_exons
                        if no_match==True:
                            if exon[-1]=="G":
                                temp_list=[i for i in exon]
                                temp_list.append({exon[-1]:[trans_ID], "R":"-"})
                                CE_exons.append(temp_list)
                            elif exon[-1]=="R":
                                temp_list=[i for i in exon]
                                temp_list.append({exon[-1]:[trans_ID], "G":"-"})
                                CE_exons.append(temp_list)
                            else:
                                print("WTF? ",  exon)
                        
            #Sort exons after start coordinate
            CE_exons.sort(key=lambda x: int(x[1]))
            #Write potential CE into output for this gene.
            for exon in CE_exons:
                out.write("{}_{}_{}_{}\t{}\t{}\t{}\t{}\t{}\n".format(exon[0], exon[3],
                                              exon[1], exon[2], "CE", gene, exon[-2], ",".join(exon[-1]["G"]), ",".join(exon[-1]["R"])))
            #If bed file wanted
            if args.bed==True:
                bed=open(args.out+"CE.bed", "w")
                for exon in CE_exons:
                    chrom, start, stop, strand = entry[0:4]
                    name="CE_"+chrom+"_"+start+"_"+stop
                    bed.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, start, stop, name, strand))
                bed.close()
        else:
            #Intron Retention, because invalid inputs are checked at the beginning of the script to save time
            #Retention of introns is usually not annotated. So we do the same as for CE.
            #And check for every gap between exons, if there is read showing intron retention.
            #That means at this point we save all potential gaps.
            
            #Remove transcript ids and duplicate exons.
            IR_coord=[]
            for trans_ID in gene_dict[gene]:
                for i in range(0, len(gene_dict[gene][trans_ID])-1):
                    chrom=gene_dict[gene][trans_ID][i][0]
                    strand=gene_dict[gene][trans_ID][i][3]
                    db=gene_dict[gene][trans_ID][i][-1]
                    IR=[chrom,strand,gene_dict[gene][trans_ID][i][2],gene_dict[gene][trans_ID][i+1][1], db, trans_ID]
                    #Check if already in list
                    no_match=True
                    for entry in IR_coord:
                        if IR[0:-2]== entry[0:-2]:
                            #theres a match
                            no_match=False
                            #if its already there, update db and Transcript ID
                            #if the db is not in the entries db, then it needs to be added to db and dict.
                            if IR[-2] not in entry[-2]:
                                #add the db strings
                                entry[-2]=entry[-2]+IR[-2]
                                #replace "-" in dictionary, initialize transcript list for new db
                                entry[-1][IR[-2]]=[IR[-1]]
                            else:
                                #just add the trans ID to the dict, as db is already ok.
                                entry[-1][IR[-2]].append(IR[-1])

                    #If there was no match found, then we have a completely new exon and append it to gene_exons
                    if no_match==True:
                        temp_list=IR[0:-1]
                        if IR[-2]=="G":
                            temp_list.append({"G":[IR[-1]], "R":"-"})
                            IR_coord.append(temp_list)
                        elif IR[-2]=="R":
                            temp_list.append({"R":[IR[-1]], "G":"-"})
                            IR_coord.append(temp_list)
                        else:
                            print("WTF ", IR)
                        
            #Sort IR entries by start of intron
            IR_coord.sort(key=lambda x: x[2])
            #Write IR entries into output file
            for IR in IR_coord:
                out.write("{}_{}_{}_{}\t{}\t{}\t{}\t{}\t{}\n".format(IR[0], IR[1], IR[2], IR[3], "IR", gene, IR[-2], ",".join(IR[-1]["G"]), ",".join(IR[-1]["R"])))
            
            #If bed file wished for
            if args.bed==True:
                bed=open(args.out+"IR.bed", "w")
                for gene in IR_coord:
                    strand=IR_coord[0].split("_")[0]
                    
                    for entry in IR_coord:
                        chrom = entry.split("_")[1]
                        start=entry.split("_")[2]
                        stop=entry.split("_")[3]
                        name="IR_"+chrom+"_"+start+"_"+stop
                        bed.write("{}\t{}\t{}\t{}\t.\t{}\n".format(chrom, start, stop, name, strand))
                bed.close()
                    
    #Printing progress
    current_count+=1
    percentage=100*(current_count/total_count)
    print("Finding AS events: {:.2f}%".format(percentage), end="\r")
                
print("Finding AS events: Done!     \n", end="\r")

out.close()

#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))  



































































    
