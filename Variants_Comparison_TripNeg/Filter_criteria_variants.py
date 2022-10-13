# -*- coding: utf-8 -*-
"""
Date: Thu Oct 13 11:04:58 2022
File Name: Filter_criteria_variants.py
Author: Mirjam Karlsson-Müller

Description: 
    To find out which variants are filtered out at which criteria, a new
    output table is created with additional columns for each criteria and 
    information whether this specific variant passes the criteria or not.
    
    Note that theres only one table returned here, instead of a filtered and
    unfiltered version as for Variants_AS_altered.py (on which this script is
    based on.)
    This is again meant to be run on the 252 samples that are shared between
    our cohort and the TNBC cohort.
    
    
List of Functions:
    
Procedure: 
    1. For shorter runtime, genes are put in a nested dictionary based on what
    chromosome and what are on the chromosome they lay on. To this nested 
    dictionary structure we add the exon information for each gene. so we have
    chromosomes: areas : genes : exons.
    2. Find out which RNA files to use: for each sample choose the vcf file
    containing the most variants, as well as its corresponding gene.tsv file
    containing the gene expression data.
    3. Read in said gene expression information (gene: FPKM)
    4. Read in variants from chosen .vcf file, checking whether the variant
    is exonic, intronic or intergenic and if its in a gene add the corresponding
    FPKM value.
    5. Go through filter criteria and add for each criteria information on whether
    the variant passes it or not.
    
    
    
Useage:
    python Filter_criteria_variants.py -s ../../Sample_Data/ -db ../Database/gencode.v38.annotation.gff3.gz -c ../Database/chromosome_sizes.txt -o filter_criteria_vcf.tsv
    
Possible Bugs:
    
    
"""


#%% Imports

import argparse
import glob
import re
import gzip
import time

#%% Time

start_time=time.time()

#%% 0. argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='VCF Filter Criteria',
                                 usage='%(prog)s -s INPUT-FOLDERS -o OUTPUT \
                                     -c CHROMOSOME FILE -db ANNOTATED EXON GTF',
                                 description="""Creates a genotype table out of
                                 several samples. Containing location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='file containing sample folder names containing among \
                        others, the vcf, bam and gene.tsv files.')
parser.add_argument('--chromosome', '-c', type=str,
                    help="""GENCODE chromosome sizes""")
parser.add_argument('--database', '-db', type=str,
                    help="""GENCODE annotated exons hg38""")
parser.add_argument('--output', '-o', required=True,
                    help="Output file")


args = parser.parse_args()

#%% Custom Functions

def Location_FPKM(chrom, position, sample):
    """

    Parameters
    ----------
    variant_position : int
        where a variant is positioned on its chromosome

    Returns
    -------
    Whether the variant lies within an exon, intron or between genes.

    """
    #If the position is not found within the coordinates of any gene, its between genes.
    location="Intergenic"
    fpkm="NA"
    #To save computation time, go through areas instead of whole chromosome.
    for area in chromosomes[chrom]:
        #Check if variant is in area.
        if int(position)>= int(area.split("-")[0]) and int(position)<=int(area.split("-")[1]):
            #If it is, check if variant is within gene.
            for gene in chromosomes[chrom][area]:
                if int(position)>= int(gene_coordinates[gene][0]) and int(position)<= int(gene_coordinates[gene][1]):
                    #If it is, then the location is definitely intron, if not found in exon.
                    location="Intron"
                    if gene in tsv_info[sample]:
                        fpkm=tsv_info[sample][gene]
                    #Go through exons
                    for exon in chromosomes[chrom][area][gene]:
                        if int(position)>=int(exon[0]) and int(position)<=int(exon[1]):
                            #If position is within exon coordinates, set location to exon.
                            location="Exon"
                            break
                    break
            break
    
    return(fpkm, location)



#%% 1. Make chromosome dictionary

print("Making chromsome dictionary...", end="\r")

chrom_sizes=dict()
with open(args.chromosome, "r") as chromosome:
    for line in chromosome:
        chrom, size= line.strip("\n").split("\t")
        chrom_sizes[chrom]=int(size)
        
        
#Split the sizes into areas (easier to read.)
chromosomes=dict()
for chrom in chrom_sizes:
    i=0
    while i<chrom_sizes[chrom]:
        start=i
        stop=i+1000000
        if stop>=chrom_sizes[chrom]:
            stop=chrom_sizes[chrom]
        if chrom not in chromosomes:
            chromosomes[chrom]=dict()
        
        chromosomes[chrom][str(start)+"-"+str(stop)]=dict()
        i=stop

print("Making chromsome dictionary: Done! \n", end="\r")
# ready to add genes and their coord.

#%% 2. exon coordinates using gtf file.

print("Extracting gene and exon coordinates...", end="\r")
db = gzip.open(args.database, "rt")

annotated_exons=dict()
gene_coordinates=dict()
for line in db:
    if line.startswith("#"):
        continue
    chrom, typ, start, stop, strand, info=list(map(line.split("\t").__getitem__,[0, 2, 3, 4, 6, 8]))
    gene_name=re.search(r"gene_name=(.+?);", info).group(1)
    
    if typ=="gene":
        gene_found=False
        for area in chromosomes[chrom]:
            """To account for the case where a gene overlaps the border between
            regions, we sort the genes to regions based only on their "smaller"
            coordinate"""
            if int(start)>=int(area.split("-")[0]) and int(start)<=int(area.split("-")[1]):
                chromosomes[chrom][area][gene_name]=[]
                gene_found=True
                gene_coordinates[gene_name]=[start, stop]
                
        if gene_found==False:
            print("gene ", gene_name, " has not been found on its chromosome.")
    
    if typ=="exon":
        exon_found=False
        for area in chromosomes[chrom]:
            for gene in chromosomes[chrom][area]:
                if gene==gene_name:
                    chromosomes[chrom][area][gene].append([start, stop])
                    exon_found=True
                    break
            if exon_found==True:
                break
print("Extracting gene and exon coordinates Done! \n", end="\r")


#%% 2. Determining which version of sample files to use.

"""Find the vcf file with the most variants per sample, choose corresponding
 gene.tsv files """

#Find all the vcf files in the data folder.
argument_glob=args.samples+"/**/*.vcf.gz"
vcf_file_list=glob.glob(argument_glob, recursive=True)

#If no vcf files are found, quit the program.
if len(vcf_file_list)==0:
    print("""There were no vcf files found in the input folder. Please make 
          sure to use the right input folder. The vcf files can be in any 
          subfolder of the input folder.""")
    quit()

"Just counting, no processing. To avoid redundance."
most_variants=dict()
sample_names=[]
#Progress update
total_files=len(vcf_file_list)
current_number_files=0
percentage=100*(current_number_files/total_files)
print("Counting variants: {:.2f}%".format(percentage), end="\r")
for file in vcf_file_list:
    currently_highest=0
    vcf = gzip.open(file, "rt")
    sample_name=re.search(r"/(S\d+)/",file).group(1)
    if sample_name not in sample_names:    
        sample_names.append(sample_name)
        most_variants[sample_name]=""
    #count lines
    counter=0
    for line in vcf:
        #Skip headers
        if line.startswith("#"):
            continue
        else:
            counter+=1
            
    if counter> currently_highest:
        currently_highest=counter
        most_variants[sample_name]=file
        

    "Progress updates on number of vcf files."
    current_number_files+=1
    percentage=100*(current_number_files/total_files)
    print("Counting variants: {:.2f}%".format(percentage),end="\r")
        
      
print("Counting variants: Done!            \n",end="\r")

sample_names=sorted(sample_names)
            
"Extract filepaths of corresponding bam and gene.tsv"

tsv_list=[]

for sample in most_variants:
    tsv_file="/".join(most_variants[sample].split("/")[0:-2])+"/t/gene.tsv"
    tsv_list.append(tsv_file)

#%% 3. Read in gene_expression files, put in dictionary with one entry per sample.


#progress updates
total_files=len(tsv_list)
current_file=0
percentage=100*current_file/total_files
print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")


tsv_info=dict()
for file in tsv_list:
    tsv_sample_name=file.split("/")[-5]
    tsv_info[tsv_sample_name]=dict()
    with open(file, "r") as tsv:
        for line in tsv:
            #if line starts with E then its a gene id
            if line.startswith("E"):
                gene_name=line.split("\t")[1]
                fpkm=line.split("\t")[7]
                tsv_info[tsv_sample_name][gene_name]=fpkm
                
    current_file+=1
    percentage=100*current_file/total_files
    print("Reading in gene.tsv: {:.2f}%".format(percentage),end="\r")

print("Reading in gene.tsv: Done!            \n", end="\r")

#%% 3. Read in variants and filter.
#Find all the vcf files in the data folder.
argument_glob=args.samples+"/**/*.vcf.gz"
vcf_file_list=glob.glob(argument_glob, recursive=True)

#If no vcf files are found, quit the program.
if len(vcf_file_list)==0:
    print("""There were no vcf files found in the input folder. Please make 
          sure to use the right input folder. The vcf files can be in any 
          subfolder of the input folder.""")
    quit()

#Progress update
total_files=len(vcf_file_list)
current_number_files=0
percentage=100*(current_number_files/total_files)
print("Extracting variants {:.2f}%".format(percentage), end="\r")
with open(args.output, "w") as out:
    out.write("Sample\tchrom\tposition\tref\talt\tgenotype\tdbsnp\tLocation\tFPKM\
              MSI\tHMPOL\tGC\tVD\tlow_complexity\tquality\tucsc_rep\n")
    for file in vcf_file_list:
        vcf = gzip.open(file, "rt")
        sample_name=file.split("/")[-5]
        #read through entries
        for line in vcf:
            #Skip headers
            if line.startswith("#"):
                continue
            chrom, position, ID, ref, alt, qual, filt, info, form, sample =line.split("\t")
            genotype=sample.split(":")[-7].strip(" ")
            
            #The tnbc data removes germline variants. So we mark them for comparison.
            if re.search(r"dbsnp_ID=.+?;", info):
                dbsnp="Yes"
            else:
                dbsnp="No"
            
            #We sort out unlocated entries.
            if chrom not in chromosomes:
                continue
            
            #check whether variant is within an exon, intron or between genes.
            fpkm,location=Location_FPKM(chrom,position, sample_name)
            
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                sample_name, chrom, position, ref, alt, genotype, dbsnp, 
                location, fpkm))
            
            "Now filter the entries."
            #Keep values with MSI<7
            #If we find the information, 
            if re.search(r"MSI=(\d+);", info):
                #check if fulfills failing criteria
                if int(re.search(r"MSI=(\d+);", info).group(1))>=7:
                    out.write("\tFail")
                #Else write a pass for this criteria
                else:
                    out.write("\tPass")
            #If the information cannot be found, counts as Fail.
            else:
                out.write("\tFail")
            #Same principle for next 3 criteria.
            
            #Keep values with HMPOL<6
            if re.search(r"HMPOL=(\d+);", info):
                if int(re.search(r"HMPOL=(\d+);", info).group(1))>=6:
                    out.write("\tFail")
                else:
                    out.write("\tPass")
            else:
                out.write("\tFail")
            
            #Keep entries with GC_cont < 78%
            if re.search(r"GC_CONT=(0\.\d+);", info):
                if float(re.search(r"GC_CONT=(0\.\d+);", info).group(1))>=0.78:
                    out.write("\tFail")
                else:
                    out.write("\tPass")
            else:
                out.write("\tFail")
            
            #Keep variant depth >=5
            if re.search(r"VD=(\d+);", info):
                if int(re.search(r"VD=(\d+);", info).group(1))<5:
                    out.write("\tFail")
                else:
                    out.write("\tPass")
            else:
                out.write("\tFail")
            
            #Is a flag, so it will only be there if it applies.
            if re.search(r"low_complexity_region", info):
                out.write("\tFail")
            else:
                out.write("\tPass")
            
            #6th column, filter bad quality reads.
            if qual=="." or float(qual)<55:
                out.write("\tFail")
            else:
                out.write("\tPass")
            
            if re.search(r"ucsc_rep=([a-z]+);", info):
                #if re.search(r"ucsc_rep=([a-z]+);", info).group(1)=="segdup":
                out.write("\tFail\n")
            else:
                out.write("\tPass\n")
            
        
        
        "Progress updates on number of vcf files."
        current_number_files+=1
        percentage=100*(current_number_files/total_files)
        print("Extracting variants: {:.2f}%".format(percentage),end="\r")
            
print("Extracting variants: Done!            \n",end="\r")





#%% End time

print("Run time: {:.2f} seconds.".format(time.time()-start_time))    
