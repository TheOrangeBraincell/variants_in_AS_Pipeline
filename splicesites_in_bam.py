# -*- coding: utf-8 -*-
"""
date: 29.03.22
name: splicesites_in_bam.py
author: Mirjam Karlsson-Müller


"""


#%%
"1. For now: Take sam file, find unique splice sites:"

import argparse
import re
import subprocess

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
"Functions"

def which_strand(flag):
    binary_flag=bin(int(flag))
    if str(binary_flag)[-8:-4] in ["0110", "1001", "b110"]:
        return "-"
    elif str(binary_flag)[-8:-4] in["1010", "0101", "b101"]:
        return "+"
       

#%%

"1. Find splice sites in SAM file, save their information in dictionary"
#Extract information from SAM:
format_splicesite= re.compile(r'\d+M\d+N\d+M')
splicesite=dict()
#counts per read for junctions:
read_counts_junction=dict()
if args.bamfile:
    if args.coordinates:
        command=['samtools','view','-h','-F','0x100', args.bamfile,
                 args.coordinates]
    else:
        command=['samtools','view','-h','-F','0x100', args.bamfile]
        
    p=subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                       stderr=subprocess.STDOUT)
    output=p.communicate(input=args.bamfile.encode())[0]
    sam=output.decode('utf-8')
    sam=sam.split("\n")
    for line in sam:
        if line.startswith("@"):
            continue
        else:
            if bool(re.search(r"(\d+)\t([a-z]{3}\d+)\t(\d+)\t(\d+)\t"
                              r"(\d+M\d+N\d+M)*\t.+YT:Z:([A-Z]{2})", 
                              line))==True:
                entries=re.search(r"(\d+)\t([a-z]{3}\d+)\t(\d+)\t(\d+)\t"
                                  r"(\d+M\d+N\d+M)\t.+YT:Z:([A-Z]{2})",line)
                #Determines whether we have a forward or reverse read
                flag=entries.group(1)
                chrom=entries.group(2)
                start=entries.group(3)
                cigar=entries.group(5)
                pair=entries.group(6)
                #Only include concordant pairs:
                if pair=="CP":
                    #Find coordinates for junction:
                    current_cigar=cigar
                    current_length=0
                    #To allow for several junctions in one cigar string.
                    while bool(re.search(format_splicesite,
                                         current_cigar))==True:
                        """Only process part of cigar string which matches the
                        format (softclips removed)"""
                        junction=re.search(r'(\d+)M(\d+)N(\d+)M',current_cigar)
                        exon1=int(junction.group(1))
                        intron=int(junction.group(2))
                        exon2=int(junction.group(3))
                        exon1_end=int(start)+exon1+current_length-1
                        exon2_start=exon1_end+intron +1 #+1 to fit helenas.
                        #Remove junctions with less than 3 bp on either side:
                        if exon1>2 and exon2>2:
                            """if coordinates are given, only include reads in
                            that region:"""
                            if args.coordinates:
                                #If the whole junction is in range of gene:
                                if exon1_end>coord_start \
                                    and exon2_start<coord_stop:
                                        #Generate data for bed file
                                        exon_count="2"
                                        block_sizes="1,1"
                                        name=str(exon1_end -1)+"_"+str(exon2_start)
                                        #Exclude non-primary alignments:
                                        if len(bin(int(flag)))>=11:
                                            #Non-primary has flag 256. 
                                            if str(bin(int(flag)))[-9]=="1":
                                                continue
                                        strand=which_strand(flag)
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
                                        splicesite[name]=[chrom, exon1_end -1, #exon1_end -1 to fit helenas
                                                  exon2_start, strand, exon_count,
                                                  block_sizes, "0,"+
                                                  str(exon2_start-exon1_end)]
                            else:
                                #Generate data for bed file
                                exon_count="2"
                                block_sizes="1,1"
                                #Flag for being in database, default =False
                                in_db=False
                                name=str(exon1_end -1)+"_"+str(exon2_start)
                                splicesite[name]=[chrom, exon1_end -1, #exon1_end -1 to fit helenas
                                          exon2_start, strand, exon_count,
                                          block_sizes, "0,"+
                                          str(exon2_start-exon1_end)]
                                
                        #remove used part of the string
                        current_cigar=re.sub(r'^.*?N','N', 
                                             current_cigar).lstrip("N")
                        current_length+=exon2_start

"""
print("name\t+\t-")
for name in read_counts_junction:
    print(name+"\t"+str(read_counts_junction[name][0])+"\t"+str(read_counts_junction[name][1])+"\t"+
          splicesite[name][3]+"\t"+splicesite[name][6])
"""

"Check for splice junctions found, what strand the majority aligns with"
for name in splicesite:
    if read_counts_junction[name][0]>read_counts_junction[name][1]:
        splicesite[name][3]="+"
    else:
        splicesite[name][3]="-"

"""
for name in splicesite:
    print("{}\t{}\t{}\t{}\t{}\t".format(name, splicesite[name][1], splicesite[name][2],
                                        splicesite[name][3],splicesite[name][6]))
"""
#%%
"2. Read in known splice junctions."

db_junctions=dict()
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
            """The data points from the bam file are in pairs, so it is 
            necessary to split the exon groups into pairs to be able to compare
            the two"""
            exon_size=exon_size.split(",")
            rel_start=rel_start.split(",")
            current_start=start
            for i in range(0, int(number_exons)-1):
                exon1_end=current_start+int(exon_size[i])-1
                exon2_start=start+int(rel_start[i+1])+1
                #if input coordinates, restrict region of interest:
                if args.coordinates:
                    if exon1_end>coord_start and exon2_start<coord_stop and \
                        chrom==coord_chrom:
                        name=str(exon1_end)+"_"+str(exon2_start)
                        #-1 because exon2_start has a +1 to exclude end of intron.
                        curr_rel_start="0,"+ str(exon2_start-exon1_end-1) #str(exon2_start-current_start)
                        #update current_start with exon2_start
                        current_start=exon2_start-1
                        #write information into dictionary.
                        db_junctions[name]=[chrom, exon1_end, exon2_start, 
                                            gene_name, strand, "2", 
                                            exon_size[i:i+2], curr_rel_start]
                else:
                    name=str(exon1_end)+"_"+str(exon2_start)
                    curr_rel_start="0,"+ str(exon2_start-current_start)
                    #update current_start with exon2_start
                    current_start=exon2_start-1
                    #write information into dictionary.
                    db_junctions[name]=[chrom, exon1_end, exon2_start, 
                                        gene_name, strand, "2", 
                                        exon_size[i:i+2], curr_rel_start]
                    
"""if we are only interested in a gene, find out which strand it is on,
by checking on what strand the majority of the known splice junctions on this
region is."""
plus_count=0
minus_count=0
for name in db_junctions:
    if db_junctions[name][4]=="+":
        plus_count+=1
    else:
        minus_count+=1
if plus_count>minus_count:
    strand="+"
else:
    strand="-"


#%%

"3. Write the results in BED file"
with open(args.out, 'w') as out:
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
    names=sorted(list(set(list(db_junctions.keys())+list(splicesite.keys()))))
    for name in names:
        if name in db_junctions:
            if name in splicesite:
                RGB="0,0,255"
            else:
                RGB="0,255,0"
            if args.coordinates:
                #Exclude entries from less frequent strand on given region.
                if db_junctions[name][4]==strand:
                    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        db_junctions[name][0], db_junctions[name][1], 
                        db_junctions[name][2], db_junctions[name][3], 0, 
                        strand, db_junctions[name][1],
                        db_junctions[name][2], RGB, db_junctions[name][5], 
                        ','.join(db_junctions[name][6]), db_junctions[name][7]))
            else:
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    db_junctions[name][0], db_junctions[name][1], 
                    db_junctions[name][2], db_junctions[name][3], 0, 
                    db_junctions[name][4], db_junctions[name][1],
                    db_junctions[name][2], RGB, db_junctions[name][5], 
                    ','.join(db_junctions[name][6]), db_junctions[name][7]))
        else:
            RGB="255,0,0"
            if args.coordinates:
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    splicesite[name][0], splicesite[name][1], splicesite[name][2], 
                    name, 0, splicesite[name][3], splicesite[name][1], 
                    splicesite[name][2], RGB, splicesite[name][4], 
                    splicesite[name][5], splicesite[name][6]))  
            else:
                """strand ="?", because it is not clear based on the bam file
                what strand the read is from. """
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    splicesite[name][0], splicesite[name][1], splicesite[name][2], 
                    name, 0, splicesite[name][3], splicesite[name][1], 
                    splicesite[name][2], RGB, splicesite[name][4], 
                    splicesite[name][5], splicesite[name][6])) 
      

                                
