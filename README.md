# Developing a method to study the effect of germline variants and somatic mutations in breast cancer patients on mRNA splicing

## Description
<p>
    This is a pipeline developed to analyze the effect of genetic variants in breast cancer patients on mRNA splicing. It can be divided into 5 steps: </p>

1. Identifying potential Casette Exons in GENCODE's annotated exons.
2.   Extracting reads from patients alignment files, spanning these exons. Calculate Percent Spliced In (PSI) scores per potential Casette Exon.
3.   Extracting variant/locationinformation from the same patient samples' .vcf files. 
4.   Fill gaps to create a complete genotype table.
5.   Identify the variants laying within the coordinates of potential Casette Exons.

The scripts for these steps, as well as the script detecting novel splice sites, have been written in Python v3.8.13.
Statistical testing has been conducted on an example gene using R v4.1.1. 

### Requirements

#### Input Files:
 * GENCODE .tsv file, containing gene ID's and their corresponding transcript ID's from UCSC Table Browser https://genome.ucsc.edu/cgi-bin/hgTables.
 * Human Genome .bed file, containing annotated exons from UCSC Table Browser https://genome.ucsc.edu/cgi-bin/hgTables.
 * sample set of alignment files (.bam/.bai) from the SCAN-B initiative.
 * sample set of variant calling files (.vcf) from the SCAN-B initiative.
 * sample set of gene information files (gene.tsv) from the SCAN-B initiative.

#### Python Packages:
* pysam (v0.19.0), has to be installed separately via conda.
```
conda install -c bioconda pysam 
```
* argparse
* glob
* time
* re
* gzip
* OrderedDict, from collections

## 1. and 2. Identifying potential Casette Exons and calculating PSI
```
usage: Find GENECODE events -s INPUT-FOLDER -o OUTPUT [-c] "chrX:XXXXXX-XXXXXX" -db BEDFILE  -gc GENCODE_TSV

Returns potential Casette Exons and their PSI scores, based on sample bam files.

optional arguments:
    -h, --help                                    show this help message and exit
    --samples SAMPLES, -s SAMPLES                 folder containing sample folders containing among others, the bam files.
    --gencode GENCODE, -gc GENCODE                tsv file containing gene_ID and transcript IDs.
    --database DATABASE, -db DATABASE             BED file containing already known/annotated splice junctions from the ucsc
    --out OUT, -o OUT                             Output txt file, where PSI scores are printed to.
    --coordinates COORDINATES, -c COORDINATES     Start and stop coordinates of region of interest, as well as chromosome. Format: chr[]:start-stop 
```
This was done using the script *events_gencode.py*. It finds potential CE in gencode bed file, and calculates the PSI for these events, based on the .bam files from the SCAN-B database. These contain mRNA sequences, aligned to a human reference genome.Further required inputs are GENCODE's .tsv and .bed file. The script can be run with coordinates, to look at a specific section of the genome
```
#Example Coordinates for Estrogen Receptor, ESR1.
python events_gencode.py -s . -gc geneID_hg38_GENCODE39.tsv -db hg38_GENCODE39.bed -o CE_ESR1.txt -c "chr6:151600000-152150000"
```
or without coordinates.
```
python events_gencode.py -s . -gc geneID_hg38_GENCODE39.tsv -db hg38_GENCODE39.bed -o CE_PSI_all_out.txt  
```
The output is a tab-delimited table containing PSI scores, having the different samples as columns and the different exons as rows. For every new gene, there is a header line, separating its exons from the previous gene's exons. 

## 3. Parse .vcf to create a variant table
```
usage: VCF Parser -s INPUT-FOLDER -o OUTPUT [-c] "chrX:XXXXXX-XXXXXX"

Creates a variant table out of several samples. Containing location x sample.

optional arguments:
    -h, --help                                   show this help message and exit
    --samples SAMPLES, -s SAMPLES                folder containing sample folders containing among others, the vcf files.
    --out OUT, -o OUT                            Output txt file, containing genotype table, positions with no variant in a sample but variants in other sample(s) are marked with NAN.
    --coordinates COORDINATES, -c COORDINATES    Start and stop coordinates of region of interest, as well as chromosome. Format: chr[]:start-stop
```

With the script *VCF_parser.py*, variants are extracted from the .vcf file for each sample. The variants are filtered and then their location and reference extracted. The script can again be run with coordinates, to only inspect a region of the genome.
```
#Example Coordinates for Estrogen Receptor, ESR1.
python VCF_parser.py -s . -o vcf_parser_ESR1.txt -c "chr6:151690496-152103274"
```
Note that the coordinates are selected much more accurately here, as we want to avoid extracting too many variants around the gene of interest, which would be excluded in the next step anyways. The script can also run for the whole genome.
```
python VCF_parser.py -s . -o vcf_parser_all.txt
```
The output is a tab delimited table containing variants. The first and second column contain location and reference respectively, whereas the remaining columns contain the samples genotype. Since most locations only have a variant in a few samples, there is a lot of gaps in the table.

## 4. Fill the gaps in the genotype table
```
usage: Genotype Table -s INPUT-FOLDER -o OUTPUT -t VCFPARSER TABLE

Fills gaps in VCF_parser.py's output table.

optional arguments:

    -h, --help                      show this help message and exit
    --samples SAMPLES, -s SAMPLES   folder containing sample folders containing among others, their gene.tsv files.
    --out OUT, -o OUT               Output txt file, containing updated genotype table, no more NAN.
    --table TABLE, -t TABLE         output from vcf parser. Table containing variants, empty spots held by NAN.
```
This was done using the script *genotype_table.py*. The table from the previous step was read, and for each *NAN* in the table, it was checked whether the variant was within a gene using the *gene.tsv* files.
```
python genotype_table.py -s . -o genotype_out.txt -t vcf_parser_all.txt
```
For variants in a gene, we assume from here on that the gene is expressed for this sample, but not having a variant at this spot. Thus the *NAN* is replaced with *HMZ* for homozygote. If the variant is not in a gene, then we assume it is not expressed *NE*, as we are working with mRNA sequences.
## 5. Extract variants from genotype table, which lay within a potential Casette Exon
```
usage: Variants in Exons [-h] --out OUT --genotype GENOTYPE --PSI_table PSI_TABLE

Filters out variants not in exons from genotype table.

optional arguments:
    -h, --help                             show this help message and exit
    --out OUT, -o OUT                      Output txt file, containing only varians in exons, includes a column for exon coordinates.
    --genotype GENOTYPE, -g GENOTYPE       txt file containing genotype per sample/location.
    --PSI_table PSI_TABLE, -p PSI_TABLE    output from events_gencode. Table containing PSI scores per sample and exon. 
```
Lastly, using the script *variants_in_exons.py*, variants not within potential Casette Exons get filtered out, and for those kept, a new column with the exon coordinates is added.

```
python variants_in_exons.py -g genotype_out.txt -p CE_ESR1_all.txt -o genotypes_in_exons.txt
```
The output is a genotype table, containing information on what position, and in which exon a variant is laying. As well as which genotype all samples have at this position.

## Statistical testing on output tables for ESR1
A set of one-way ANOVA's were conducted on the tables, so that for each exon-location pair, it was tested whether there was a significant difference between genotypes. The p-values of these tests where then corrected for multiple testing. The script used for this was statistics_ESR1.R. Note that to run this script yourself, the working directory and the table names have to be changed in the script. Both the results and the script can be found in the Statistical_Testing subfolder.

## Further steps
Using the PSI table and the final genotype table, one can analyze the potential correlation between the expression of an exon and the genetic variants found in it.


## Identification of Novel Splice Junctions
```
usage: Find Splicesites -s SAM-INPUT -o OUTPUT

Returns coordinates of unique splicesites in SAM file.

optional arguments:
    -h, --help                                   show this help message and exit 
    --bamfile BAMFILE, -b BAMFILE                BAM file, containing alignment information.
    --database DATABASE, -db DATABASE            BED file containing already known/annotated splice junctions from the ucsc.
    --out OUT, -o OUT                            Output bed file, where comparison data should be printed to.
    --coordinates COORDINATES, -c COORDINATES    Start and stop coordinates of region of interest, as well as chromosome. Format: chr[]:start-stop 
```
A final goal of this project would be to be able to apply this pipeline not only to annotated splice junctions and exons, but also to novel splice junctions found in the .bam files of the SCAN-B initiative.
To identify such splice junctions, a script was written. It reads an alignment file and extracts coordinates of reads, spanning splice junctions. They are then compared to the gencode database, to see which splice junctions are found in both, .bam and GENCODE .bed, and which in either. The output .bed file is colorcoded.
```
#with coordinates for Estrogen Receptor
python splice_junctions.py -b alignment.bam -db hg38_GENCODE39.bed -o test_out.txt -c "chr6:151690496-152103274"

#no coordinates
python splice_junctions.py -b alignment.bam -db hg38_GENCODE39.bed -o test_all.txt"
```

