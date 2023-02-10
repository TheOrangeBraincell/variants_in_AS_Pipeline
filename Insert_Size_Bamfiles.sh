#!/bin/bash -l

#SBATCH -A sens2019019 -p core -n 1 -t 14:00 --mail-user=mi3258mu-s@student.lu.se --mail-type=END,FAIL --J Insert_Size

module load python/3.8.7
start=`date +%s`


#Getting insert size stats from all bam files in subfolders
#for i in ../Sample_Data/*/*/*/*.bam; do samtools stats $i | grep -A 1 "insert size average"; done
#
#New Plan 
#
#1 Get 3' UTR coordinates from gencode gff files: i.e. last exons with size >1000 bp. Samtools wants bedfile as input.
#So the entries are written into bedfile.
#
cat Database/gencode.v39.annotation.gff3| grep "three_prime_UTR" | awk 'BEGIN{OFS="\t"} ($5-$4>1000) {print ($1 OFS $4 OFS $5)}'>3_UTR.bed 
#
#2 Extract and filter reads from those regions with samtools (Uniquely aligned, proper pair, map quality)
#3 Take subsample of reads for each sample ID.
#4 Average and variance for all of the insert sizes of those subsamples of reads.
#
for i in ../Sample_Data/*/*/*/*.bam; do samtools view $i -L Database/3_UTR.bed -h -f PROPER_PAIR -F UNMAP,SECONDARY,QCFAIL --subsample 0.25 | head -n 500000| samtools stats | grep -A 1 "insert size average"; done >average_insert.txt #
#
#5 Now the averages have to be averaged (beautiful word) and the standard deviations as well
awk 'BEGIN{FS="\t"} (NR %2 == 1) {s+=$3} (NR %2 == 0) {sd+=$3^2} END{print("Mean", s/(NR/2) , "Standard Deviation", sqrt(sd/(NR/2)))}' average_insert.txt
#Mean X Standard Deviation Y
