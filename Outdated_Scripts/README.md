# Developing a method to study the effect of germline variants and somatic mutations in breast cancer patients on mRNA splicing

## Description
<p>
    This is a pipeline developed to analyze the effect of genetic variants in breast cancer patients on mRNA splicing. It can be divided into 5 steps: </p>

1. Identifying potential Casette Exons in GENCODE's annotated exons.
2.   Extracting reads from patients alignment files, spanning these exons. Calculate Percent Spliced In (PSI) scores per potential Casette Exon.
3.   Extracting variant/locationinformation from the same patient samples' .vcf files. Assign genotypes where not extracted from variant table.
4.   Identify the variants laying within the coordinates of potential Casette Exons.

The scripts for these steps, as well as the script detecting novel splice sites, have been written in Python v3.8.13.
Statistical testing has been conducted on an example gene using R v4.1.1. 

### Requirements

#### Input Files:
 * Human Genome .tsv file GENCODE, containing annotated exons from UCSC Table Browser https://genome.ucsc.edu/cgi-bin/hgTables.
 * Human Genome .tsv file RefSeq, containing annotated exons from UCSC Table Browser https://genome.ucsc.edu/cgi-bin/hgTable
 * sample set of alignment files (.bam/.bai) from the SCAN-B initiative.
 * sample set of variant calling files (.vcf) from the SCAN-B initiative.
 * sample set of gene expression information files (gene.tsv) from the SCAN-B initiative.

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

## 1. and 2. Identifying potential Casette Exons and calculating PSI
```
usage: Find Alternative Splicing events -s INPUT-FOLDER -o OUTPUT -g GENCODE-TSV -r REFSEQ-TSV [-c] "chrX:XXXXXX-XXXXXX" [-n] "ABCDE"                                                                                                                                                                                           Returns potential Casette Exons based on GENCODE39 and RefSeq and their PSI scores, based on sample bam files.                                                   

optional arguments:
-h, --help                                  show this help message and exit
--samples SAMPLES, -s SAMPLES               folder containing sample folders containing among others, the bam files.
--gencode GENCODE, -g GENCODE               tsv file containing bed file information on annotated exons from GENCODE39 as well as gene names.
--refseq REFSEQ, -r REFSEQ                   tsv file containing bed file information on annotated exons from RefSeq as well as gene names.
--out OUT, -o OUT                           Output tsv file, where PSI scores are printed to.
--coordinates COORDINATES, -c COORDINATES   Start and stop coordinates of region of interest, as well as chromosome. Format: chr[]:start-                                                                   stop
--name NAME, -n NAME                        Symbol of gene of interest. f.e. ESR1. (Symbol by HUGO Gene Nomenclature Committee).
--AS AS, -as AS                             Which type of alternative splicing event we are interested in. "CE" for Casette Exons, "AA" for alternative                                                     acceptors, "AD" for alternative donors, "IR" for intron retention and "ALL" for all of the types. Several seperated                                             by ,.
```
This was done using the script *AS_events.py*. It finds potential CE in both RefSeq and GENCODE tsv files, and calculates the PSI for these events, based on the .bam files from the SCAN-B database. These contain mRNA sequences, aligned to the human reference genome. Further required inputs are GENCODE's and RefSeq .tsv. The script can be run with coordinates, to look at a specific section of the genome.
```
#Example Coordinates for Estrogen Receptor, ESR1.
python AS_events.py -s ../Sample_Data/ -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -o CE_ESR1_cor.tsv -c "chr6:151650000-152130000" -as CE
```
or without coordinates.
```
python AS_events.py -s ../Sample_Data/ -g Database/hg38_GENCODE39_all.tsv -r Database/hg38_NCBI_all.tsv -o CE_ESR1_sym.tsv -n ESR1 -as CE
```
The output is a tab-delimited table containing PSI scores, having the different samples as columns and the different exons as rows. For every new gene, there is a header line, separating its exons from the previous gene's exons. 

## 3. Parse .vcf to assign genotypes per variant location and sample.
```
usage: VCF Parser -s INPUT-FOLDER -o OUTPUT [-c] "chrX:XXXXXX-XXXXXX"                                                                                           

Creates a genotype table out of several samples. Containing location x sample.                                                                                   

optional arguments:
-h, --help                                  show this help message and exit
--samples SAMPLES, -s SAMPLES               folder containing sample folders containing among others, the vcf files.
--out OUT, -o OUT                           Output tsv file, containing genotype table.
--coordinates COORDINATES, -c COORDINATES   Start and stop coordinates of region of interest, as well as chromosome. Format: chr[]:start-stop
```

With the script *AS_events.py*, variants are extracted from the .vcf file for each sample. The variants are filtered and then their genotype HETZ or HMZA (Heterozygous Allele or Homozygous Alternative Allele) extracted. Resulting genotype table is filled by assigning HMZR (Homozygous Reference Allele) or NOEXP (Not expressed). The script can again be run with coordinates, to only inspect a region of the genome.
```
#Example Coordinates for Estrogen Receptor, ESR1.
python VCF_Genotypes.py -s . -o genotype_ESR1.tsv -c "chr6:151690496-152103274"
```
Note that the coordinates are selected much more accurately here, as we want to avoid extracting too many variants around the gene of interest, which would be excluded in the next step anyways. The script can also run for the whole genome.
```
python VCF_Genotypes.py -s . -o genoytpe_wholegenome.tsv
```
The output is a tab delimited table containing genotypes. The first and second column contain location and reference respectively, whereas the remaining columns contain the samples genotype.

## 4. Extract variants from genotype table, which lay within a potential Casette Exon
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

