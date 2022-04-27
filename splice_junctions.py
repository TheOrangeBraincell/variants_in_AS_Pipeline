# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 15:18:11 2022

@author: mirja
"""


import argparse
import re
import pysam

#%%
"Functions"

def which_strand(flag):
    """Takes a SCANB sam file, determines based on sam flag, on which strand 
    the read is."""
    binary_flag=bin(int(flag))
    if str(binary_flag)[-8:-4] in ["0110", "1001", "b110"]:
        return "-"
    elif str(binary_flag)[-8:-4] in["1010", "0101", "b101"]:
        return "+"

#%%

"1. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Find Splicesites',
                                 usage='%(prog)s -s SAM-INPUT -o OUTPUT ',
                                 description="""Returns coordinates of unique  
                                              splicesites in SAM file.""")

parser.add_argument('--bamfile', '-b', required=True,
                    help='BAM file, containing alignment information.')

parser.add_argument('--database', '-db', required=True,
                    help="""BED file containing already known/annotated splice
                    junctions from the ucsc.""")
parser.add_argument('--out', '-o', required=True,
                    help="""Output bed file, where comparison data should be 
                    printed to.""")
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")

args=parser.parse_args()

#Extract input coordinates:
if args.coordinates:
    if  bool(re.search(r'[a-z]{3}[MXY]?\d*:\d+-\d+', args.coordinates))==False:
        raise argparse.ArgumentTypeError("""The given coordinates are in the 
                                         wrong format. Input as 
                                         chrX:XXXX-XXXX.""")
        quit()
    else:
        coord=re.search(r'([a-z]{3}[MXY]?\d*):(\d+)-(\d+)', args.coordinates)
        coord_chrom=coord.group(1)
        coord_start=int(coord.group(2))
        coord_stop=int(coord.group(3))
        

#%%

"2a. Reading in ucsc bed file into a dictionary"

db_junctions=dict()
# exons_plus=[]
# exons_minus=[]
# #introns_plus=[]
# #introns_minus=[]
with open(args.database, 'r') as db:
    for line in db:
        #To exclude potential title lines/empty lines, formatting mistakes
        #Only takes chr[] and chr[]_random lines, in accordance with bam.
        if bool(re.search(r"([a-z]{3}[X,M,Y]?\d*)\t(\d+)\t(\d+)\t([A-Z]+\d+"
                          r".*\d*)\t0\t(\-?\+?)\t\d+\t\d+\t0\t(\d+)\t([\d,]+)"
                          r"\t([\d,]+)", line))==True or bool(re.search(
                            r"([a-z]{3}[X,M,Y]?\d*).+_random\t(\d+)\t(\d+)\t"
                            r"([A-Z]+\d+\.*\d*)\t0\t(\-?\+?)\t\d+\t\d+\t0\t"
                            r"(\d+)\t([\d,]+)\t([\d,]+)", line))==True:
            #specify what the groups in the line correspond to.
            entry=re.search(r"([a-z]{3}[X,M,Y]?\d*).*\t(\d+)"
                        r"\t(\d+)\t([A-Z]+\d+\.*\d*)\t0\t(\-?\+?)\t\d+\t"
                        r"\d+\t0\t(\d+)\t([\d,]+)\t([\d,]+)", line)
            chrom=entry.group(1)
            start=int(entry.group(2))
            stop=entry.group(3)
            gene_name=entry.group(4)
            strand=entry.group(5)
            number_exons=entry.group(6)
            exon_size=entry.group(7)
            rel_start=entry.group(8)
            
            """To be able to compare it to bam file, we require single splice 
            junctions between 2 exons."""
            exon_size=exon_size.split(",")
            rel_start=rel_start.split(",")
            current_start=start
            for i in range(0, int(number_exons)-1):
                #exon1_start=current_start+1
                exon1_end=current_start+int(exon_size[i])-1
                exon2_start=start+int(rel_start[i+1])+1
                #exon2_end=exon2_start+exon_size[i+1]-1
                "For 5."
                # #Make exon/intron lists, seperately per strand.
                # if strand=="-":
                #     exons_minus.append([exon1_end-1, exon1_start-1],
                #                        [exon2_end-1, exon2_start-1])
                #     #introns_minus.append([exon1_start-1, exon2_end-1])
                # else:
                #     exons_plus.append([exon1_start, exon1_end], [exon2_start, exon2_end])
                #     #introns_plus.append([exon1_end+1, exon2_start-1, strand])
                """if we have input coordinates, check if the exons are in the 
                right region. If they arent, skip them."""
                if args.coordinates:
                    if exon1_end<coord_start or exon2_start>coord_stop or \
                        chrom!=coord_chrom:
                        continue
                    
                """"For each pair of relevant exons, we write the information 
                into the dict""" 
                name=str(exon1_end)+"_"+str(exon2_start)
                curr_rel_start="0,"+ str(exon2_start-current_start)
                #update current_start with exon2_start
                current_start=exon2_start-1
                #write information into dictionary.
                db_junctions[name]=[chrom, exon1_end, exon2_start, 
                                    gene_name, strand, "2", 
                                    exon_size[i:i+2], curr_rel_start]
                        
                
"""2.b. If we are only interested in one specific gene, then we do not need the 
splice junctions found on the other strand. So in those cases, we check on what
strand the majority of the exons can be found and then exclude exons on the 
other."""

if args.coordinates:
    plus_count=0
    minus_count=0
    #count pluses and minuses.
    for name in db_junctions:
        if db_junctions[name][4]=="+":
            plus_count+=1
        else:
            minus_count+=1
    #check for majority
    if plus_count>minus_count:
        strand="+"
    else:
        strand="-"               
    #remove entries on the wrong strand from the dictionary.
    wrong_strand=[]
    for name in db_junctions:
        if db_junctions[name][4]!=strand:
            wrong_strand.append(name)
    for name in wrong_strand:
        del db_junctions[name]


#%% Identify splicesites in bam file
"""3. Identify splicesites in bam file."""

"3a. Read in information"

#Splicejunction has alignment pattern:
pattern_junction=re.compile(r'\d+M\d+N\d+M')
#initialize dictionary for information
read_junction=dict()
#initialize dictionary for counts of strands. 
read_counts_junction=dict()

#Transform bam file into samfile:
if args.bamfile:
    sam=pysam.AlignmentFile(args.bamfile, "rb")
    if args.coordinates:
        sam=sam.fetch(coord_chrom, coord_start, coord_stop)
    for line in sam:
        line=line.to_string()
        if line.startswith("@"):
            continue
        #if the line follows the pattern containing one or more splice junctions:
        if bool(re.search(r"(\d+)\t([a-z]{3}\d+)\t(\d+)\t(\d+)\t"
                          r"(\d+M\d+N\d+M)*\t.+YT:Z:([A-Z]{2})", 
                          line))==True:
            entries=re.search(r"(\d+)\t([a-z]{3}\d+)\t(\d+)\t(\d+)\t"
                              r"(\d+M\d+N\d+M)\t.+YT:Z:([A-Z]{2})",line)
            #Assign variables to entries in line.
            flag=entries.group(1)
            chrom=entries.group(2)
            start=entries.group(3)
            cigar=entries.group(5)
            pair=entries.group(6)
            
            #We only include concordant pairs:
            if pair!="CP":
                continue
            
            """To allow for several junctions in one cigar string, we require a loop that keeps 
            looking for a pattern."""
            current_cigar=cigar
            current_length=0
            while bool(re.search(pattern_junction, current_cigar))==True:
                junction=re.search(r'(\d+)M(\d+)N(\d+)M',current_cigar)
                exon1=int(junction.group(1))
                intron=int(junction.group(2))
                exon2=int(junction.group(3))
                exon1_end=int(start)+exon1+current_length-1
                exon2_start=exon1_end+intron +1
                
                """If junction has less than 3 bp on either side of the intron, 
                remove it:"""
                if exon1<3 or exon2<3:
                    #update cigar string
                    current_cigar=re.sub(r'^.*?N','N', 
                                         current_cigar).lstrip("N")
                    current_length+=exon2_start
                    continue
                
                #if there's input coordinates, exons need to be in the region.
                if exon1_end<coord_start or exon2_start>coord_stop:
                    #update cigar string
                    current_cigar=re.sub(r'^.*?N','N', 
                                         current_cigar).lstrip("N")
                    current_length+=exon2_start
                    continue
                
                #Exclude non-primary alignments (flag 256)
                if len(bin(int(flag)))>=11 and str(bin(int(flag)))[-9]=="1":
                    #update cigar string
                    current_cigar=re.sub(r'^.*?N','N', 
                                         current_cigar).lstrip("N")
                    current_length+=exon2_start
                    continue
                
                #Generate information for dictionary:
                name=str(exon1_end -1)+"_"+str(exon2_start)
                strand=which_strand(flag)
                #Add strand information to strand count dictionary:
                if name in read_counts_junction:
                    if strand=="-":
                        read_counts_junction[name][1]+=1
                    else:
                        read_counts_junction[name][0]+=1
                else:
                    if strand=="-":
                        read_counts_junction[name]=[0,1]
                    else:
                        read_counts_junction[name]=[1,0]
                
                #Add the information into the dictionary:
                read_junction[name]=[chrom, exon1_end-1, exon2_start, strand, 
                                     "0,"+ str(exon2_start-exon1_end)]
                #update cigar string
                current_cigar=re.sub(r'^.*?N','N', 
                                     current_cigar).lstrip("N")
                current_length+=exon2_start
              
"""
#To look at read_junction dictionary and strand counts in a table
print("name\t+\t-")
for name in read_counts_junction:
    print(name+"\t"+str(read_counts_junction[name][0])+"\t"+str(
        read_counts_junction[name][1])+"\t"+read_junction[name][3]+
        "\t"+read_junction[name][6])
"""

"3b Check for each junction, on which strand the majority of its reads lays"

for name in read_junction:
    if read_counts_junction[name][0]>read_counts_junction[name][1]:
        read_junction[name][3]="+"
    else:
        read_junction[name][3]="-"
        
#%% Write bed output with known and novel splicesites.
            
"""4. Summarize both dictionaries into a table of splice junctions, into bed 
format. Note that to be able to display the bed file in the genome browser,
we set exon sizes to 1, and relative starts to 0, intron size. This way there
is also uniformity with the novel entries, as they do not have information
such as exon sizes. The exon sizes etc, will be used further down to classify
the splice junctions."""


with open(args.out, 'w') as out:
    "Write headline for displaying the file in genome browser."
    if args.coordinates:
        out.write("browser position "+ coord_chrom+ ":"+str(coord_start)+"-"+
                  str(coord_stop)+"\ntrack name=\"Splice Junctions\""
                  +"description=\"Splice Junctions for specific region with" 
                  +" RGB\" visibility=2 itemRgb=\"On\"\n")
    else:
        out.write("track name=\"Splice Junctions\""
                  +"description=\"Splice Junctions for genome with" 
                  +"RGB\" visibility=2 itemRgb=\"On\"\n")
    
    #make combined name list, to iterate through.
    names=sorted(list(set(list(db_junctions.keys())+list(
        read_junction.keys()))))
    exon_count="2"
    block_size="1,1" #Because we do not know for the novel splicesites.
    for name in names:
        if name in db_junctions:
            if name in read_junction:
                RGB="0,0,255"
            else:
                RGB="0,255,0"
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"\
                      .format( db_junctions[name][0], db_junctions[name][1], 
                       db_junctions[name][2], db_junctions[name][3], 0,
                       db_junctions[name][4], db_junctions[name][1], 
                        db_junctions[name][2], RGB, exon_count, block_size, 
                        db_junctions[name][7]))
        else:
            RGB="255,0,0"
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"\
                      .format(read_junction[name][0], read_junction[name][1], 
                              read_junction[name][2], name, 0, 
                              read_junction[name][3], read_junction[name][1], 
                              read_junction[name][2], RGB, exon_count, 
                              block_size, read_junction[name][4]))
            



#%% Classifying annotated splice junctions (Commented because atm not used.)

# """5. To classify novel splice junctions, we want to create a dataset
# with information about the annotated splice junctions"""

# #From the exons and introns lists in 2.

# "5a. classify AA and AD"

# def AA_AD(exons):
#     """ finds alternative donors and acceptors using a list of exon 
#     coordinates. The output will be written into two dictionaries."""
#     starts=set()
#     ends=set()
#     aa=dict()
#     ad=dict()
#     for coordinates in exons:
#         current_start=coordinates[0]
#         current_end=coordinates[1]
#         if current_end in ends:
#             #ad found!
#             if current_end in ad and current_start in ad[current_end]:
#                 continue
#             elif current_end in ad:
#                 ad[current_end].append(current_start)
#             else:
#                 ad[current_end]=[current_start]
#                 #also add pair with which the end matched.
#                 ad[current_end].append(coor[0] for coor in exons if 
#                                        coor[1]==current_end)
#         elif current_start in starts:
#             #aa found!
#             if current_start in aa and current_end in aa[current_start]:
#                 continue
#             elif current_start in aa:
#                 aa[current_start].append(current_end)
#             else:
#                 aa[current_start]=[current_end]
#                 #also add pair with which the start matched.
#                 aa[current_start].append(coor[1] for coor in exons if
#                                          coor[0]==current_start)
#         starts.add(current_start)
#         ends.add(current_end)
#     return aa, ad, starts, ends

# "For positive strand"
# aa_plus, ad_plus, starts_plus, ends_plus=AA_AD(exons_plus)

# "For negative strand"
# aa_minus, ad_minus, starts_minus, ends_minus=AA_AD(exons_minus)


    
# "5b. classify CE"
# #based on introns lists.

# "positive strand"
# ce_plus=dict()
# for starts, ends in zip(starts_plus, ends_plus):
#     if starts>introns_plus[0] and ends<introns_plus[1]:
#         #CE found.
#         ce_plus[starts]=ends

# "negative strand"
# ce_minus=dict()
# for starts, ends in zip(starts_minus, ends_minus):
#     if starts < introns_minus[0] and ends>introns_minus[1]:
#         #CE found
#         ce_plus[starts]=ends



# "5c write results into output file"

# """
# #output sorted by starting position of junction (or lowest starting position.)

# #idk how to do this elegantly. so we gonna do it on tuesday ._.  
# #merging dictionaries doesnt work, starts might be included in different SE.
# #Start list, sorted? Then check for each of the 6 dictionaries if it is in,
# #if yes, write it to annotate file. That would be a lot of code though. 

# #also i did not save chromosome or strand info with it... only exons. 
# #wait. my dictionaries of aa, and ad only contain exons. not junctions. 
# #just exons..... so i classified exons, not splice junctions.
# #is that what i want?

# #should i link the exon classification to the splice junction the exon is
# #involved in?
# what if two exons involved in one splice junction are classified differently? 
# # can i classify the splice junctions directly? Follow up in lab journal..
# """

# with open("annotate.txt", 'w') as annotate:
#     annotate.write("#Chrom\tstart(s)\tstop(s)\tstrand\ttype")




        

























