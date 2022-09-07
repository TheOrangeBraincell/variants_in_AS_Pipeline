# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:40:14 2022

python AF_vs_TR.py -s ../Sample_Data/ -c "chr6:1-170805979" 

@author: mirjam
"""

"""The goal is to get an idea how high the gene expression/how many reads there
need to be for us to be able to make a reliable genotype prediction. 
(Though we will eventually probably consult the DNA data)"""

#%% Imports

import argparse
import glob
import re
import gzip
import matplotlib.pylab as plt


#%% Argparse
"0. Setting up argparse, handling input parameters"

parser = argparse.ArgumentParser(prog='Genotype Plots',
                                 usage='%(prog)s -s INPUT-FOLDER -o OUTPUT \
                                     [-c] "chrX:XXXXXX-XXXXXX"',
                                 description="""Creates a variant table out of
                                 several samples. Containing location x sample.""")

parser.add_argument('--samples', '-s', required=True,
                    help='folder containing sample folders containing among \
                        others, the vcf files.')
parser.add_argument('--coordinates', '-c', type=str,
                    help="""Start and stop coordinates of region of interest,
                    as well as chromosome. Format: chr[]:start-stop""")
parser.add_argument('--out', '-o', required=True,
                    help="""Output txt file.""")

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


#%% Idea 3: allele fractions/total read depth.

#need only variant files.

#Find all the vcf files in the data folder.
argument_glob=args.samples+"/**/*.vcf.gz"
vcf_file_list=glob.glob(argument_glob, recursive=True)

#If no vcf files are found, quit the program.
if len(vcf_file_list)==0:
    print("""There were no vcf files found in the input folder. Please make 
          sure to use the right input folder. The vcf files can be in any 
          subfolder of the input folder.""")
    quit()

#Files are zipped though.
#Progress update
total_files=len(vcf_file_list)
current_number_files=0
percentage=100*(current_number_files/total_files)
print("Filtering vcf {:.2f}%".format(percentage), end="\r")
variants=dict()
sample_names=[]
for file in vcf_file_list:
    sample_name=file.split("/")[-5]
    if sample_name not in sample_names:
        sample_names.append(sample_name)
    vcf = gzip.open(file, "rt")
    #read through entries
    for line in vcf:
        #Skip headers
        if line.startswith("#"):
            continue
        chrom, position, ID, ref, alt, qual, filt, info, form, sample =line.split("\t")
        genotype=sample.split(":")[-7].strip(" ")
        
        #If coordinates are given:
        if args.coordinates:
            if chrom != coord_chrom:
                continue
            else:
                if int(position)< int(coord_start) or int(position)> int(coord_stop):
                    continue
    
        #print(info)
        "Filter entries based on info string"
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
        

        "Add coordinates, for those who pass filtering"
        if sample_name in variants:
            variants[sample_name].append([re.search(r"DP=(\d+);", info).group(1),re.search(r"AF=(\d+.\d+);", info).group(1), chrom, position])
        else:
            variants[sample_name]=[[re.search(r"DP=(\d+);", info).group(1),re.search(r"AF=(\d+.\d+);", info).group(1), chrom, position]]
        
        
        
    current_number_files+=1
    percentage=100*(current_number_files/total_files)
    print("Filtering vcf: {:.2f}%".format(percentage),end="\r")

print("Filtering vcf: Done!            \n",end="\r")




with open(args.out, "w") as out:
    out.write("Sample\tDP\tAF\tchrom\tpos\n")
    for sample_name in variants:
        for variant in variants[sample_name]:
            out.write(sample_name+"\t"+"\t".join(variant)+"\n")

""" Plots are to be made in R
#make plot
x=[]
y=[]

for sample in variants:
    for variant in variants[sample]:
        x.append(float(variant[0]))
        y.append(float(variant[1]))
        if int(variant[0])>100 and float(variant[1])<0.3:
            print(sample, variant)



#print(x)
#print("\n\n")
#print(y)

plt.plot(x,y,"ro",mfc='none')
plt.ylabel("Allele Fraction")
plt.xlabel("Total reads")
plt.title("Allele Fraction vs total reads")
plt.savefig("AlleleFraction_TotalReads.png")

plt.plot(x,y,"ro", mfc='none')
axis=[i for i in range(1,int(max(x)))]
plt.plot(axis,[5/i for i in axis])
plt.ylabel("Allele Fraction")
plt.xlabel("Total reads")
plt.xlim(0,50)
plt.ylim(0,1)
plt.title("Allele Fraction vs total reads, xlim")
plt.savefig("AlleleFraction_TotalReads_xlim.png")

"""






