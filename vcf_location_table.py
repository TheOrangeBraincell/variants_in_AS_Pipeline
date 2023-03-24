# -*- coding: utf-8 -*-
"""
Date: Thu Nov  3 13:01:32 2022
File Name: vcf_location_table.py
Author: Mirjam Karlsson-Müller

Description:
    Parses through all variant calling files from a cohort and makes a list of all
    variants found (filtered and unfiltered).
    Outputs a table sample x location.
    
List of Functions:
    
Procedure: 
    1. Argparse
    2. Read in variant calling files, generate table (incomplete) of all found locations/variants.
    3. Print said table into a file.
    
Useage: Run in command line
    
    With coordinates for f.e. estrogen receptor:
        python vcf_location_table.py -s ../Sample_Data/ -o location_table_ESR1.tsv -c "chr6:151656691-152129619"
    
    For whole genome:
        python vcf_location_table.py -s ../Sample_Data/ -o location_table_WG.tsv
    
    Server:
        python vcf_location_table.py -s /raidset/mi3258mu-se/Mirjam -o location_table_WG.tsv
    
Possible Bugs:
    
    
"""

#%% Imports
import argparse
import glob
import gzip
import os
import re
import time

#%% Start Timer

start_time=time.time()


#%% Argparse

"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='vcf location table',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     [-c] "chrX:XXXXXX-XXXXXX"',
                                 description="""Creates a location of variants table out of
                                 several samples. Containing location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing the vcf files.')
parser.add_argument('--out', '-o', required=True,
                    help="""Output file, containing vcf location table.""")
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")


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

""" Not useful for batch runs
#Check if output file already exists.
while True:
    if os.path.isfile(args.out):
        #If it does, ask user if they want to overwrite the output.
        answer=str(input("The output file "+ args.out+ " already exists. Would you like to overwrite it? (Y/N) "))
        if answer.upper()=="Y":
            break
        elif answer.upper()=="N":
            new_output=str(input("Please enter a new output name: "))
            if all([e in [" ", "", "\t", "\n"] for e in new_output]):
                print("This is not a valid output name. Lets try this again.")
            else:
                output=new_output
        else:
            print("This is not a valid response. Please answer with Y (yes) or N (no).")
    else:
        #There is no problem with this output file name. We proceed with the code.
        break
"""
#%% Read vcf file

#Find all the vcf files in the data folder.
argument_glob=args.samples+"/*/*/*/*/*.vcf.gz"
vcf_file_list=glob.glob(argument_glob, recursive=True)

#If no vcf files are found, quit the program.
if len(vcf_file_list)==0:
    print("""There were no vcf files found in the input folder. Please make 
          sure to use the right input folder. The vcf files can be in any 
          subfolder of the input folder.""")
    quit()


chromosomes=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X", "Y", "M"]
if args.coordinates:
    chrom_index=chromosomes.index(coord_chrom.strip("chr"))
#Initialize dictionaries
variants=dict()
#Sample list
sample_names=[]
#Progress update
total_files=len(vcf_file_list)
current_number_files=0
percentage=100*(current_number_files/total_files)
print("Reading vcf {:.2f}%".format(percentage), end="\r")
for file in vcf_file_list:
    vcf = gzip.open(file, "rt")
    sample_name=re.search(r"/(S\d+)/",file).group(1)
    if sample_name not in sample_names:    
        sample_names.append(sample_name)
    #read through entries
    for line in vcf:
        #Skip headers
        if line.startswith("#"):
            continue
        chrom, position, ID, ref, alt, qual, filt, info, form, sample =line.split("\t")
        genotype=sample.split(":")[-7].strip(" ")
        #We only want single point mutations
        if len(alt)>1 or len(ref)>1:
            continue
        
        #If coordinates or name are given, restrict area of variants.
        if args.coordinates:
            #we can assume that the chromosomes at least are sorted by appearance. Which means if we are at a chromosome after the one we look for, 
            #we can stop reading the file.
            if chrom != coord_chrom:
                if chrom_index<chromosomes.index(chrom.strip("chr")):
                    break
                else:
                    continue
            
            else:
                if int(position)< int(coord_start):
                    continue
                elif int(position)> int(coord_stop):
                    break
        
        "Before filtering, add all variants to genotype table."
        variant_ID=chrom+"_"+position
        #add to variant dict
        if variant_ID not in variants:
            variants[variant_ID]=dict()

        if sample_name not in variants[variant_ID]:
            variants[variant_ID][sample_name]=[ref, {alt:"ND"}]
        
        if alt not in variants[variant_ID][sample_name][1]:
            variants[variant_ID][sample_name][1][alt]="ND"
        
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
            variants[variant_ID][sample_name][1][alt]="0/1"
        elif genotype=="1/1":
            variants[variant_ID][sample_name][1][alt]="1/1"
        #although it should not happen....  But it clearly is.
        elif genotype=="0/0":
            variants[variant_ID][sample_name][1][alt]="0/0"
        else:
            print("Invalid genotype ", genotype)
    
    "Progress updates on number of vcf files."
    current_number_files+=1
    percentage=100*(current_number_files/total_files)
    print("Reading vcf: {:.2f}%".format(percentage),end="\r")
        
        
print("Reading vcf: Done!            \n",end="\r")

sample_names=sorted(sample_names)

#sort variants after chromosome and location.
variant_keys=list(variants.keys())
sorted_chr_variants=[]

#chromosomes have a bit of a specific order, so we do it manually...
for chromosome in chromosomes:
    for key, value in variants.items():
        if key.startswith("chr"+chromosome):
            sorted_chr_variants.append(key)

#For positions i dont need to do it manually.
sorted_variants=sorted(sorted_chr_variants, key=lambda x: int(x.split("_")[1]))

sorted_dict=dict()
for v in sorted_variants:
    sorted_dict[v]=variants[v]

#to not have to change everything below, we overwrite the variable.
variants=sorted_dict
#delete big thingies.
del sorted_dict

#%% Write output file

counter=0
invalid_counter=0
percentage=100*counter/len(list(variants.keys()))
print("Writing Location Table: {:.2f}%".format(percentage),end="\r")

with open(args.out, "w") as out:
    #Create Header
    if args.coordinates:
        out.write("#Variant Location Table for variants in range "+ args.coordinates + "\n")
    else:
        out.write("#Variant Location Table for variants in whole genome.\n")
    out.write("#Identifiers in first column in format chr_position_reference nucleotide_(alternative nucleotide(s))\n")
    out.write("#\n#Genotype Identifiers:\n#  0/0 = Homozygous Reference Allele\n#  0/1 = Heterozgyous Alternative Allele\n#  1/1 = Homozygous Alternative Allele\n")
    out.write("Location\t"+"\t".join(sample_names)+"\n")
    
    for variant in variants:
        invalid_variant=False
        #Starting string for file: want variant_ID, ref and alt in info
        #Check if ref the same for all samples. otherwise exclude variant
        ref_list=[]
        for sample in variants[variant]:
            if variants[variant][sample][0] not in ref_list:
                ref_list.append(variants[variant][sample][0])
        if len(ref_list)!=1:
            continue
        
        #Otherwise we analyse the alternative nucleotide(s) and their genotype
        #Check how many alts we have
        alt_bases=[]
        for sample in variants[variant]:
            for alt in variants[variant][sample][1]:
                if alt not in alt_bases:
                    alt_bases.append(alt)
        #Now write the whole thing.
        new_line=[variant+ "_"+ ref_list[0]+"_("+"_".join(alt_bases)+")"]
        #Go through sorted samples, to always have the same order.
        for sample in sample_names:
            if sample in variants[variant]:
                #check what genotype to append
                #Iterate through alt bases
                genotype_changed=False
                #if length of alternative alleles for one sample is 1, then we can just read them off
                if len(list(variants[variant][sample][1].keys()))==1:
                    for alt in variants[variant][sample][1]:
                        #According to genotype
                        if variants[variant][sample][1][alt]=="ND":
                            genotype="ND"
                        if variants[variant][sample][1][alt]=="0/0":
                            genotype="0/0"
                        elif variants[variant][sample][1][alt]=="0/1":
                            genotype="0/"+str(alt_bases.index(alt)+1)
                        elif variants[variant][sample][1][alt]=="1/1":
                            genotype=str(alt_bases.index(alt)+1)+"/"+str(alt_bases.index(alt)+1)

                        new_line.append(genotype)
                #If length is 2, then we need to  merge them.
                elif len(list(variants[variant][sample][1].keys()))==2:
                    
                    g1, g2= variants[variant][sample][1].values()
                    a1, a2= variants[variant][sample][1].keys()
                    #the same and homozygous reference
                    if g1==g2 and g1=="0/0":
                        genotype="0/0"
                    #not the same but genotype 1 is 0/0
                    elif g1=="0/0":
                        genotype="0/"+str(alt_bases.index(a2)+1)
                    #same but for genotype 2
                    elif g2=="0/0":
                        genotype="0/"+str(alt_bases.index(a1)+1)
                    elif g1=="ND" and g2=="ND":
                        genotype="ND"
                    elif g1=="ND":
                        genotype=g2
                    elif g2=="ND":
                        genotype=g1
                    #Now if both are 0/1 or 1/1 (which is also nonsense...)
                    else:
                        str(alt_bases.index(a1)+1)+"/"+str(alt_bases.index(a2)+1)
                    new_line.append(genotype)

                #if length is 3 or more, then i dont understand biology... that shouldnt happen. so we just print an error.
                
                elif len(list(variants[variant][sample][1].keys()))>2:
                    #print("3 alternative alleles for one sample and posiiton. That makes no sense. This is variant " +variant+ " and sample " + sample +"\n The alt bases are ", list(variants[variant][sample][1].keys())  )
                    genotype="ND"
                
            else:
                #Put a spaceholder. Genotype will be determined using gene expression data by next script in pipeline.
                new_line.append("-")
        #Skip to new variant if there is 3 or more alleles.        
        if invalid_variant==True:
            invalid_counter+=1
            
            continue
        
        #If the new line only has "-" and "ND" then theres no point having it in the output file. as no alternative genotype passed filtering.
        for i in new_line[1:]:
            if i not in ["-", "ND"]:
                #Proof that line is important. If it never finds any that isnt - or ND, it wont write it.
                out.write("\t".join(new_line)+"\n")
                break
            
        counter+=1
        percentage=100*counter/len(list(variants.keys()))
        print("Writing Location Table: {:.2f}%".format(percentage),end="\r")

print("Writing Location Table: Done          \n",end="\r")
        

#%% Stop Timer

print("Run time: {:.2f} seconds.".format(time.time()-start_time))      
