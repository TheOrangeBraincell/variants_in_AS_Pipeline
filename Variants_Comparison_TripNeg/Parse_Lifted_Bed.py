# -*- coding: utf-8 -*-
"""
name: Parse_Lifted_Bed.py
author: Mirjam Karlsson-Müller
date: 30.9.22

Description: Reformats the output bed from LiftOver of ucsc.

Usage:
python Parse_Lifted_Bed.py -v hglft_genome_34d4d_aad970.bed -db ../Database/gencode.v38.annotation.gff3.gz -c ../Database/chromosome_sizes.txt -o tnbc_variants.tsv -s ../../Sample_Data/

"""

#%% Imports
import argparse
import re
import gzip
import glob

#%% argparse

parser = argparse.ArgumentParser(prog='Parse TNBC')

parser.add_argument('--variants', '-v', required=True,
                    help='variant calling file tnbc')
parser.add_argument('--database', '-db', type=str,
                    help="""GENCODE annotated exons hg38""")
parser.add_argument('--out', '-o', required=True,
                    help="Output file")
parser.add_argument('--chromosome', '-c', type=str,
                    help="""GENCODE chromosome sizes""")
parser.add_argument('--samples', '-s', required=True,
                    help='file containing sample folder names containing among \
                        others, the vcf, bam and gene.tsv files.')


args = parser.parse_args()

#%% 1. Read gtf database file



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


# ready to add genes and their coord.
#%%exon coordinates.

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

#%% 2. Determining which version of gene.tsv sample files to use.

"""Find the vcf file with the most variants per sample, choose corresponding
bam and gene.tsv files """

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
    

#%% Reading in gene expression infos to link to gene names.
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

#%% 2. Read tnbc data and check wether they are in exon etc.


with open(args.variants,"r") as bed, open(args.out, "w") as out:
    out.write("Sample\tchrom\tposition\tref\talt\tgenotype\tlocation\tFPKM\n")
    for line in bed:
        chrom, start, stop, name= line.strip("\n").split("\t")[0:4]
        sample, genotype, ref, alt= name.split("_")
        
        #Check whether the variant is in an exon.
        location="Intergenic"
        fpkm="NA"
        #We sort out unlocated entries.
        if chrom not in chromosomes:
            continue
        for area in chromosomes[chrom]:
            if int(stop)>= int(area.split("-")[0]) and int(stop)<=int(area.split("-")[1]):
                for gene in chromosomes[chrom][area]:
                    if int(stop)>= int(gene_coordinates[gene][0]) and int(stop)<= int(gene_coordinates[gene][1]):
                        location="Intron"
                        if sample in tsv_info:
                            if gene in tsv_info[sample]:    
                                fpkm=tsv_info[sample][gene]
                        for exon in chromosomes[chrom][area][gene]:
                            if int(stop)>=int(exon[0]) and int(stop)<=int(exon[1]):
                                location="Exon"
                                break
                        break
                break
        
        
        
        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, chrom, stop, 
                                                    ref, alt, genotype, location, fpkm))
        
        