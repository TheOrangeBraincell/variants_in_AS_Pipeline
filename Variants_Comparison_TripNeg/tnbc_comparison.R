# tnbc_comparison.R
# Mirjam Karlsson-Müller
# 29.09.22
#
# Description: Comparison of variants found based on DNA samples of 252 patients,
# compared to variants found on RNA samples of the same 252 patients.
#
#

library(tidyverse)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\TNBC_comparison")

dna_variants<- read_tsv("tnbc_variants.tsv")
rna_variants_filtered<-read_tsv("rna_variants_filtered.txt")
rna_variants_unfiltered<-read_tsv("rna_variants.txt")



# How many of the variants in the dna files did we find with each versin of the rna variants?
dna_variants
#genotype 0/1 and 1/0 have no clear distance. so we set them all to 0/1

rna_variants_filtered$genotype<-str_replace(rna_variants_filtered$genotype, "1/0", "0/1")
rna_variants_unfiltered$genotype<-str_replace(rna_variants_unfiltered$genotype, "1/0", "0/1")


#Remove samples that are only in one cohort.
rna_variants_filtered %>% 
  select(Sample) %>% 
  distinct() %>% 
  pull(Sample)-> rna_samples

dna_variants %>% 
  select(Sample) %>% 
  distinct() %>% 
  pull(Sample)->dna_samples


rna_variants_filtered %>% 
  filter(Sample %in% dna_samples)->rna_vfu #rna variants filtered updated

rna_variants_unfiltered %>% 
  filter(Sample %in% dna_samples)->rna_vu

dna_variants %>% 
  filter(Sample %in% rna_samples)->dna_vu #dna variants updated


#Only focus on single base substitutions at this point:
# so exclude all lines with alt or ref longer than 1 char.
rna_vfu

rna_vfu %>% 
  filter(nchar(alt)==1 & nchar(ref)==1)->rna_vfu

rna_vu %>% 
  filter(nchar(alt)==1 & nchar(ref)==1)->rna_vu

dna_vu %>% 
  filter(nchar(alt)==1 & nchar(ref)==1)->dna_vu


#Now i should basically be able to count identical rows.

common_variants<-inner_join(dna_vu, rna_vfu)
common_variants
#6603

common_variants %>% 
  group_by(genotype) %>% 
  count()

# Lets only compare variants in exons.
dna_exons<- dna_vu %>% 
  filter(Location=="Exon")

dna_exons #ca 55000

rnaf_exons<-rna_vfu %>% 
  filter(Location=="Exon")
rnaf_exons #most of these are in exons. 4.8-> 4.1 Mio.

rna_exons<-rna_vu %>% 
  filter(Location=="Exon")
rna_exons # same. 10 Mio->8.9 Mio

common_exon_variants<-inner_join(dna_exons, rnaf_exons)
common_exon_variants
#5886 filtered, 8945 unfiltered 

#How many of the RNA variants are germline?
rnaf_exons %>% 
  group_by(dbsnp) %>% 
  count()
#Most. Only 117 624 are not germline variants.

#common number between no germline and dna?
rnaf_exons_no_germ<-rnaf_exons %>% 
  filter(dbsnp=="No")

rnaf_exons_no_germ #117 624

common_wo_germline_exon<-inner_join(dna_exons, rnaf_exons_no_germ)
common_wo_germline_exon
#4417
#So most matches are not germline variants. which makes sense given that
# germline variants are excluded from the dna data set (mostly.)

rna_exons_no_germ<-rna_exons %>% 
  filter(dbsnp=="No")

rna_exons_no_germ #117 624

common_wo_germline_exon<-inner_join(dna_exons, rna_exons_no_germ)
common_wo_germline_exon


#Relationship to Expression level?


