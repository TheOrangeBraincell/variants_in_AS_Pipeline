# Developing a method to study the effect of germline variants and somatic mutations in breast cancer patients on mRNA splicing

 Master Thesis in Bioinformatics, 60 cp.
 Mirjam Karlsson-Müller
 
#### Abbrevations
* AS = Alternative Splicing
* PSI = Percent Spliced In
* CE = Casette Exon
* IR = Intron Retention
* AD = Alternative Donor
* AA = Alternative Acceptor
 
## The Pipeline
This is a pipeline developed to analyze the effect of genetic variants in breast cancer patients on mRNA splicing. It can be divided into steps.

1. Collecting variant locations from all samples' .vcf files.
2. Determining genotypes for samples without variant at specific variant location.
3. Identifying AS events based on GENCODE and RefSeq gene annotation files. Scoring these events based on samples' .bam files.
4. Linking Events to PSI scores.
5. Statistical Testing

### Requirements
The RNA data comes from the SCAN-B cohort. Alignment, variant calling and gene expression files were used.
To write the pipeline, both Python (v.3.9.7) and R (v.4.1.1) were used.
 
#### Python Packages
* argparse
* glob
* re
* gzip
* time
* os
* pysam
* math
* currentframe, getframeinfo from inspect

pysam (v0.19.0), has to be installed separately via conda.
```
conda install -c bioconda pysam 
```

#### Input Files
* Human Genome .tsv file GENCODE, containing annotated exons from UCSC Table Browser https://genome.ucsc.edu/cgi-bin/hgTables. (hg 38)
 * Human Genome .tsv file RefSeq, containing annotated exons from UCSC Table Browser https://genome.ucsc.edu/cgi-bin/hgTable. (hg 38)
 * sample set of alignment files (.bam/.bai) from the SCAN-B initiative.
 * sample set of variant calling files (.vcf) from the SCAN-B initiative.
 * sample set of gene expression information files (gene.tsv) from the SCAN-B initiative.
 
 
### 1. Parsing through variant calling files (.vcf): Extracting Locations
For this step all .vcf files for all samples to be investigated are needed. The script *vcf_location_table.py* can be used on the whole genome or on a specific set of coordinates.
 
```
usage: vcf location table -s INPUT-FOLDER -o OUTPUT [-c] "chrX:XXXXXX-XXXXXX" 

Creates a location of variants table out of several samples. Containing location x sample.                 

optional arguments:
-h, --help                                   show this help message and exit
--samples SAMPLES, -s SAMPLES                folder containing sample folders containing the vcf files.
--out OUT, -o OUT                            Output file, containing vcf location table.
--coordinates COORDINATES, -c COORDINATES    Start and stop coordinates of region of interest, as well as  chromosome. Format: chr[]:start-stop
``` 
Run with for example:

```
#with coordinates
python vcf_location_table.py -s ../Sample_Data/ -o location_table_ESR1.tsv -c "chr6:151690496-152103274"

#without coordinates
python vcf_location_table.py -s ../Sample_Data/ -o location_table_WG.tsv
```
The output file is a location x sample table, containing genotypes where given by the .vcf files, and spaceholders "-" where there is no matching entry for a sample in a .vcf file.

### 2. Parsing through gene expression files: Assigning Genotypes
For this step, all gene expression information files for all samples in the location table created in step 1, are needed. These files contain the FPKM per gene and sample. Based on that value, the missing genotypes in the location table are assigned. This is done with the script *assigning_genotypes.py*.

```
usage: assigning genotypes -s SAMPLE-FOLDER -o OUTPUT -i INPUT

Creates a genotype table out of a location table. Containing genotypes for location x sample.

optional arguments:
-h, --help                       show this help message and exit
--samples SAMPLES, -s SAMPLES    folder containing sample folders containing gene.tsv files.
--input INPUT, -i INPUT          Input file containing location table.
--out OUT, -o OUT                Output file containing genotypes. 
```
Run with for example:
```
python assigning_genotypes.py -s ../Sample_Data/ -i location_table_ESR1.tsv -o genotype_table.tsv 
```

The result is a complete location x sample table, containing the genotype of every sample at every location.

### 3. Identifying annotated Alternative splicing events and Scoring them based on Alignment files

Based on GENCODE (v 39) and RefSeq annotation files (hg v38) alternative splicing (AS) events are identified. For the coordinates of those events, reads are extracted from the samples' alignment files and PSI scores calculated.

```
usage: AS_PSY.py -s SAMPLE-FOLDER -o OUTPUT-FOLDER -g GENCODE-FILE -r REFSEQ-FILE [-c] "chrX:XXXXXX-XXXXXX"  -as AS-TYPE -is INSERT-SIZE

Per AS event of interest, creates a table with the PSI scores supporting said event per sample in sample folder.

optional arguments:

-h, --help                                   show this help message and exit
--samples SAMPLES, -s SAMPLES                folder containing sample folders containing among others, the vcf, bam and gene.tsv files.
--out OUT, -o OUT                            Output folder, containing PSI table.
--coordinates COORDINATES, -c COORDINATES    Start and stop coordinates of region of interest, as well as chromosome. Format: chr[]:start-stop
--gencode GENCODE, -g GENCODE                tsv file containing bed file information on annotated exons from GENCODE39 as well as gene names.
--refseq REFSEQ, -r REFSEQ                   tsv file containing bed file information on annotated exons from RefSeq as well as gene names.
--AS AS, -as AS                              Which type of alternative splicing event we are interested in. "CE" for Casette Exons, "AA" for alternative acceptors, "AD" for alternative donors, "IR" for intron retention and "ALL" for all of the types. Several seperated by ,.
--InsertSize INSERTSIZE, -is INSERTSIZE      Average Insert size plus standard deviation. Format 'Mean X Standard Deviation Y. Only needed if one of the AS of interest is IR.'
```
Run with for example:
```
#with coordinates f.e. Estrogen Receptor
python variants_in_AS_Pipeline/AS_PSI.py -s ../Sample_Data/ -o PSI_ESR1/ -c "chr6:151656691-152129619" -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as ALL -is "Mean 231.467 Standard Deviation 92.8968"
```

The result is a separate table for each AS event the script was run for containing PSI scores for event location by Sample. The script also outputs a .bed file per AS event, that can be opened in f.e. igv to see the coordinates of the AS events.

## Comparing Variants found based on RNA and DNA data
 
 
### Requirements
 The RNA variants come from the SCAN-B cohort and the DNA variants can be found on https://data.mendeley.com/datasets/2mn4ctdpxp/3 (article for corresponding study at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6859071/). The comparison was made based on 248 overlapping samples.

Both Python (v.3.9.7) and R (v.4.1.1) were used.

#### Python Packages
* argparse
* glob
* re
* gzip
* time

#### R Package:
* Tidyverse

#### Other Software:
* samtools (v.1.15.1) installed via conda.

### 1. Parsing DNA Variants

The DNA data set used different sample names. So first step was to find the corresponding samples in our data and then extract all the information from the substitution and indel file correspondingly. This can be done with script *Parse_TNBC.py* run in command line.

```
usage: python Parse_TNBC.py -v VARIANTS -i IDS -o OUTPUT

optional arguments:                                                
-h, --help                           show this help message and exit
--variants VARIANTS, -v VARIANTS     variant calling file tnbc
--ids IDS, -i IDS                    Pa_ids linking to the sample names we use
--out OUT, -o OUT                    Output .bed file containing columns: chromosome, start position, stop position, info (sample+ref+alt)
```

### 2. LiftOver

The DNA data set was also using hg37 for their alignments, whereas we use hg38. So to compare the positions, they had to be "lifted over" using "LiftOver" from ucsc (https://genome.ucsc.edu/cgi-bin/hgLiftOver). The entries that were unable to be transferred to hg38 were saved in a separate file.

### 3. Creating .tsv file using gene expression data
The now correctly positioned variants, are flagged based on whether they lie in an exon, an intron or between genes. If they are in a gene, the corresponding FPKM value is extracted from the gene.tsv files from the RNA-Seq data.
```
usage: python Parse_Lifted_Bed.py -v BEDFILE -db GENCODE -c CHROMOSOME-SIZES -o OUTPUT -s SAMPLE-FOLDER -vt VARIANT-TYPE

optional arguments:
-h, --help                                    show this help message and exit
--variants VARIANTS, -v VARIANTS              variant calling file tnbc
--database DATABASE, -db DATABASE             GENCODE annotated exons hg38
--out OUT, -o OUT     Output file    
--chromosome CHROMOSOME, -c CHROMOSOME        GENCODE chromosome sizes
--samples SAMPLES, -s SAMPLES                 file containing sample folder names containing among others, the vcf, bam and gene.tsv files.
--varianttype VARIANTTYPE, -vt VARIANTTYPE    'indels' or 'substitutions'
```
The result is a .tsv file contaning columns: Sample, chromosome, position, reference base, alternative base, genotype, Location (Exon, Intron, Intergenic) and FPKM.
### 4. Extracting variants from RNA data's matching samples

```
usage: Variants_AS_altered -s INPUT-FOLDER -c CHROMOSOME-SIZES -db GENCODE

optional arguments:

-h, --help                                show this help message and exit
--samples SAMPLES, -s SAMPLES             file containing sample folder names containing among others, the vcf, bam and gene.tsv files.
--chromosome CHROMOSOME, -c CHROMOSOME    GENCODE chromosome sizes
--database DATABASE, -db DATABASE         GENCODE annotated exons hg38 
```
The result is a .tsv file contaning columns: Sample, chromosome, position, reference base, alternative base, genotype, Location (Exon, Intron, Intergenic) and FPKM.

### 5. Plotting in R
To find out what percentage of variants found in DNA data can also be found in RNA data, the script *tnbc_comparison.R* was used. This was not run from command line but in R Studio. The resulting plots can be found in *variants_in_AS_Pipeline/Variants_Comparison_TripNeg/* (Note that some of the plots were slightly altered outside of R to make the labels readable.)


### 6. Filter Criteria 
To see which filter criteria removed matching variants between DNA and RNA data, the script *Filter_criteria_variants.py* was used, which shows exactly which of the filter criteria each variant entry fails (if any). Note that this step was only done for substitutions not for indels.

```
usage: VCF Filter Criteria -s INPUT-FOLDERS -o OUTPUT -c CHROMOSOME FILE -db ANNOTATED EXON GTF

optional arguments:
-h, --help                                show this help message and exit
--samples SAMPLES, -s SAMPLES             file containing sample folder names containing among others, the vcf, bam and gene.tsv files.
--chromosome CHROMOSOME, -c CHROMOSOME    GENCODE chromosome sizes
--database DATABASE, -db DATABASE         GENCODE annotated exons hg38
--output OUTPUT, -o OUTPUT                Output file
```

### 7. Sequencing Depth
Using the script *Sequencing_Depth.R* the sequencing depth for variants in regions with high FPKM was extracted. This was done by calling samtools depth from within the R script. The resulting dataframe was then evaluated with script *Subs_depth_plots.R*. Note that this step was only done for substitutions not for indels.