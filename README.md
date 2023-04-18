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
3. Identifying AS events based on GENCODE and RefSeq gene annotation files. 
4. Scoring these events based on samples' .bam files.
5. Preparing genotype tables and PSI tables for statistics.
6. Statistical Testing

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
* Human Genome .tsv file GENCODE, containing annotated exons from UCSC Table Browser https://genome.ucsc.edu/cgi-bin/hgTables. (hg 38) (You can find screenshots of the table settings in the "How_to_get_right_Database_Tables".)
 * Human Genome .tsv file RefSeq, containing annotated exons from UCSC Table Browser https://genome.ucsc.edu/cgi-bin/hgTable. (hg 38) (You can find screenshots of the table settings in the "How_to_get_right_Database_Tables".)
 * sample set of alignment files (.bam/.bai) from the SCAN-B initiative.
 * sample set of variant calling files (.vcf) from the SCAN-B initiative.
 * sample set of gene expression information files (gene.tsv) from the SCAN-B initiative.
 
#### Computational Power
The pipeline was run on a cluster for sensitive data in Uppsala (https://www.uppmax.uu.se/resources/systems/the-bianca-cluster/). It was sped up by running the steps only for set of 10,699 genes with sufficient expression (FPKM>10).
 
### 1. Parsing through variant calling files (.vcf): Extracting Locations
For this step all .vcf files for all samples to be investigated are needed. The script *vcf_location_table.py* can be used on the whole genome or on a specific set of coordinates. The vcf file list is created with bash.
```
find [Sample Directory] -name "variants-annotated.vcf.gz" > vcf_file_list.txt
```
 The script is then according to:
```
usage: vcf location table -s INPUT-FOLDER -o OUTPUT [-c] "chrX:XXXXXX-XXXXXX" 

Creates a location of variants table out of several samples. Containing location x sample.                 

optional arguments:
-h, --help                                   show this help message and exit
--samples SAMPLES, -s SAMPLES                file containing the path to the vcf files.
--out OUT, -o OUT                            Output file, containing vcf location table.
--coordinates COORDINATES, -c COORDINATES    Start and stop coordinates of region of interest, as well as  chromosome. Format: chr[]:start-stop
``` 
Run with for example:

```
#with coordinates
python vcf_location_table.py -s vcf_file_list.txt -o location_table_ESR1.tsv -c "chr6:151690496-152103274"

#without coordinates
python vcf_location_table.py -s vcf_file_list.txt -o location_table_WG.tsv
```
The output file is a location x sample table, containing genotypes where given by the .vcf files, and spaceholders "-" where there is no matching entry for a sample in a .vcf file. Note that the entries are sorted and filtered before writing them into files. That way we only keep information on variants that have at least one alternative genotype that passes filtering.

### 2. Parsing through gene expression files: Assigning Genotypes
For this step, all gene expression information files for all samples in the location table created in step 1, are needed. These files contain the FPKM per gene and sample. To speed up the genotyping script, they were collected in a gene x sample table, containing the FPKM values. Based on that table, the missing genotypes in the location table are assigned. This is done with the script *genotype.py*. For that we require an input file giving us the coordinates of each annotated gene found in the RefSeq and GENCODE database. This we create with *gene_ranges.py*:

```
usage: Create gene range tsv -o OUTPUT-FILE -g GENCODE-FILE -r REFSEQ-FILE [-c] "chrX:XXXXXX-XXXXXX"

Create a tsv table that contains information on the chromosome, strand, min and max coordinate for every gene found in RefSeq or GENCODE.

optional arguments:
-h, --help                                   show this help message and exit
--out OUT, -o OUT                            Output file, containing gene ranges.
--coordinates COORDINATES, -c COORDINATES    Start and stop coordinates of region of interest, as well as chromosome. Format: chr[]:start-stop
--gencode GENCODE, -g GENCODE                tsv file containing bed file information on annotated exons from GENCODE39 as well as gene names.
--refseq REFSEQ, -r REFSEQ                   tsv file containing bed file information on annotated exons from RefSeq as well as gene names.
```
Now we are ready to run *genotypes.py*.
```
python genotype.py --help
usage: assigning genotypes -s INPUT-FOLDER -o OUTPUT -i LOCATION-TABLE [-c] "chrX:XXXXXX-XXXXXX" -r RANGE-TSV -f FPKM-TABLE

Creates a genotype table out of a location table. Containing genotypes for location x sample.

optional arguments:
  -h, --help                                   show this help message and exit
  --fpkm FPKM, -f FPKM                         table containing fpkm values for all genes and samples
  --input INPUT, -i INPUT                      Input file containing location table.
  --out OUT, -o OUT                            Output file containing genotypes.
  --coordinates COORDINATES, -c COORDINATES    Start and stop coordinates of region of interest, as well as chromosome. Format: chr[]:start-stop
  --ranges RANGES, -r RANGES                   File containing gene ranges created with refseq and gencode.
```
Run with for example:
```
python genotype.py -s /Sample_Data/ -i location_table_ESR1.tsv -o genotype_table_ESR1.tsv -r gene_ranges.tsv -c "chr6:151690496-152103274" -f fpkm_table.tsv
```

The result is a complete location x sample table, containing the genotype of every sample at every location. 

### 3. Identifying alternative splicing events based on gene annotation
Based on GENCODE (v 39, hg38) and RefSeq annotation files (hg v38) alternative splicing (AS) events are identified. 

```
usage: Identify Alternative Splicing -o OUTPUT-FILE                                      -g GENCODE-FILE -r REFSEQ-FILE                                          [-c] "chrX:XXXXXX-XXXXXX" -as AS-TYPE

Per AS event of interest, creates a table with the PSI scores supporting said event per sample in sample folder.

optional arguments:
  -h, --help            show this help message and exit
  --out OUT, -o OUT     Output file, containing events and type per gene.
  --coordinates COORDINATES, -c COORDINATES
                        Start and stop coordinates of region of interest, as well as chromosome. Format: chr[]:start-
                        stop
  --gencode GENCODE, -g GENCODE
                        tsv file containing bed file information on annotated exons from GENCODE39 as well as gene
                        names.
  --refseq REFSEQ, -r REFSEQ
                        tsv file containing bed file information on annotated exons from RefSeq as well as gene names.
  --AS AS, -as AS       Which type of alternative splicing event we are interested in. "CE" for Casette Exons, "AA"
                        for alternative acceptors, "AD" for alternative donors, "IR" for intron retention and "ALL"
                        for all of the types. Several seperated by ,.
  --bed BED, -b BED     If an output file bed file for each type of event is wished for, set to True
```
Run with for example:
```
#with coordinates f.e. Estrogen Receptor
 python Identify_AS.py -o AS_events_ESR1.tsv -c "chr6:151656691-152129619" -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -as ALL
```

The result is a tab separated file containing alternative splicing events and their coordinates as well as information on what transcript ID and gene name they are found in.

### 4. Scoring the identified alternative splicing events based on RNA reads
We now score the alternative splice events identified in the previous step with PSI scores using the script *PSI.py*. 

```
usage: Score Alternative Splicing -s SAMPLE-FOLDER -o OUTPUT-FOLDER -g GENCODE-FILE -r REFSEQ-FILE [-c] "chrX:XXXXXX-XXXXXX" -as AS-TYPE -is INSERT-SIZE

Creates a table with the PSI scores supporting said event per sample in sample folder.

optional arguments:
  -h, --help            show this help message and exit
  --samples SAMPLES, -s SAMPLES
                        folder containing sample folders containing among others, the vcf, bam and gene.tsv files.
  --input INPUT, -i INPUT
                        Table with all identified AS events, created by Identify_AS.py
  --out OUT, -o OUT     Output file, containing PSI tables.
  --InsertSize INSERTSIZE, -is INSERTSIZE
                        Average Insert size plus standard deviation. Format 'Mean X Standard Deviation Y'
  --coordinates COORDINATES, -c COORDINATES
                        Start and stop coordinates of region of interest, as well as chromosome. Format: chr[]:start-
                        stop
  --AS AS, -as AS       Which type of alternative splicing event we are interested in. "CE" for Casette Exons, "AA"
                        for alternative acceptors, "AD" for alternative donors, "IR" for intron retention and "ALL"
                        for all of the types. Several seperated by ,.
```
Run with for example:
```
python PSI.py -s /Sample_Data -i AS_events_ESR1.tsv -c "chr6:151656691-152129619" -is "Mean 231.467 Standard Deviation 92.8968" -as ALL -o ESR1_PSI.tsv
```
The result is a tab separated file containing the location of the event and its type in the first column, and the Samples PSI scores in the remaining columns. 

### 5. Statistical testing of results
Before statistical testing, the PSI tables and genotype tables were filtered to only include entries that are relevant for testing. This was done with the *Prep_Genotype_Stats.py* and *Prep_PSI_Stats.py* scripts respectively.

```
python Prep_Genotype_Stats.py INPUT-FILE OUTPUT-FILE
python Prep_PSI_Stats.py INPUT-FILE OUTPUT-FILE
```
The resulting tables have the same layout, but only rows passing filtering conditions.
Note that the statistics script was run per gene, in an attempt to not overwhelm R.
```
Rscript Statistics.R PSI-TABLE GENOTYPE-TABLE GENE-NAME
```
This produces a tsv file which contains the gene, variant location, AS_event location and p values. The p values are then corrected with the script *FDR.R* using the total number of tests as input.
```
Rscript FDR.R NUMBER-TESTS INPUT-FILE
```
The resulting file is identical to the previous, but with an additional column with q values, the for multiple testing corrected p values.


### Novel Splice Sites (Not a part of the pipeline)
The script *splice_junctions.py* has been written to identify novel splice sites based on patient .bam files. This is how to run it on a bam file:
```
 python splice_junctions.py -b alignment.bam -db hg38_GENCODE39.bed -o test_out.txt -c "chr6:151690496-152103274"
```
Due to the time constraint, this approach has not been pursued further. 

## Comparing Variants found based on RNA and DNA data
 
 This is not a part of the pipeline, but decisions for the pipeline have been made based on this comparison. So it is included here with stepwise instructions.
 
### Requirements
 The RNA variants come from the SCAN-B cohort and the DNA variants can be found on https://data.mendeley.com/datasets/2mn4ctdpxp/3 (article for corresponding study at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6859071/). The comparison was made based on 248 overlapping samples.

Both 
(v.3.9.7) and R (v.4.1.1) were used.

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