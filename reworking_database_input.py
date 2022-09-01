# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 08:15:01 2022

@author: mirjam

#With coordinates
python reworking_database_input.py -s ../Sample_Data/ -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -o CE_ESR1_cor.tsv -c "chr6:151650000-152110000" -as CE 

#with gene symbol
 python reworking_database_input.py -s ../Sample_Data/ -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -o CE_ESR1_sym.tsv -n ESR1 -as CE


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
                                     -g GENCODE-TSV -r REFSEQ-TSV \
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
parser.add_argument('--AS', '-as', type=str, required=True,
                    help="""Which type of alternative splicing event we are
                    interested in. "CE" for Casette Exons, "AA" for alternative
                    acceptors, "AD" for alternative donors, "IR" for intron
                    retention and "ALL" for all of the types. Several seperated
                    by ,.""")

args = parser.parse_args()

# Extract input coordinates:
if args.coordinates:
    if bool(re.search(r'[a-z]{3}[MXY]?\d*:\d+-\d+', args.coordinates))== False:
        raise argparse.ArgumentTypeError("""The given coordinates are in the 
                                         wrong format. Input as 
                                         chrX:XXXX-XXXX.""")
        quit()
    else:
        coord = re.search(r'([a-z]{3}[MXY]?\d*):(\d+)-(\d+)', args.coordinates)
        coord_chrom = coord.group(1)
        coord_start = int(coord.group(2))
        coord_stop = int(coord.group(3))

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


#%% 1. Database Input.

#outer dictionary with gene names as key
gene_dict=dict()

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

print("Creating Database Dictionary: Done!")

#%% 2. Alternative Splicing events pt 1

for event in inputs:
    
    "2a.) Casette Exons, Preparation"
    
    if event.upper()=="CE":
        
        #Remove first and last exon for each transcript, as they cannot be CE
        for gene in gene_dict:
            #if a transcript only has <=2 exons, there can be no CE.
            filtered_gene_dict={k : v for k , v in gene_dict[gene].items() 
                                if len(v) >2}
            gene_dict=filtered_gene_dict
            print(gene_dict.keys())
            for trans_ID in gene_dict[gene]:
                #remove "first" exon, smallest start coordinate.
                gene_dict[gene][trans_ID].sort(key=lambda x: x[1])
                first_exon=gene_dict[gene][trans_ID][0]
                gene_dict[gene][trans_ID][1::]
                #remove "last" exon, biggest end coordinate.
                gene_dict[gene][trans_ID].sort(key=lambda x: x[2])
                last_exon=gene_dict[gene][trans_ID][-1]
                gene_dict[gene][trans_ID][0:-1]
                
        
        """Transcript IDs are no longer required, and additionally there
        is definitely duplicate exons between R and G as well as within each
        of them. Those need to be removed."""
        
        gene_exons=dict()
        exons_db=dict()
        
        for gene in gene_dict:
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

#%% 3. Extract reads from bam files (SCAN-B)
"""
###############################################
#This part has not been changed since spring! #
###############################################
"""

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
    elif args.name:
        samfile=samfile.fetch(first_exon[0], int(first_exon[1]), int(last_exon[2]))
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


#%% 4. Alternative Splicing events pt 2:

"this next part is in no way adjusted for whats to come"
"IT will have to be reworked as soon as theres more AS"
"Also solve the output file problem. How to mark different events?"

for event in inputs:
    
    "2a.) Casette Exons"
    
    if event.upper()=="CE":
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
                    
                    for read in sample_dict[sample]:
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
                out.write("{}\t{}\n".format(exon[0]+"_"+str(exon[1])+"_"+str(exon[2]), 
                                     "\t".join(PSI_scores)))


        out.close
                     
                        

        



print("Counting Reads: Done!         \n", end="\r")


#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))
        














