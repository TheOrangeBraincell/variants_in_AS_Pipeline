#!/usr/bin/env Rscript 
#27.10.22
# Sequence Depth Script

library(tidyverse)

#setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\TNBC_comparison")

dna_variants<- read_tsv("tnbc_variants_fpkm.tsv")
rna_variants_filtered<-read_tsv("rna_variants_filtered.txt")
rna_variants_unfiltered<-read_tsv("rna_variants.txt")

#Adjust heterozygous genotype format
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
  filter(Sample %in% dna_samples) %>% 
  filter(nchar(alt)==1 & nchar(ref)==1) %>% 
  filter(Location=="Exon") %>% 
  mutate(Match ="filtered")->rna_subs_f #rna substitutions filtered

rm(rna_variants_filtered)

rna_variants_unfiltered %>% 
  filter(Sample %in% dna_samples) %>% 
  filter(nchar(alt)==1 & nchar(ref)==1) %>% 
  filter(Location=="Exon")->rna_subs

rm(rna_variants_unfiltered)

dna_variants %>% 
  filter(Sample %in% rna_samples) %>% 
  filter(nchar(alt)==1 & nchar(ref)==1) %>% 
  filter(Location=="Exon")->dna_subs #dna substitutions

rm(dna_variants)

tibble(full_Path=list.files(full.names=TRUE, recursive=TRUE)) %>%
  filter(str_detect(full_Path, pattern="alignment.bam")) %>%
  mutate(Sample=substring(full_Path,3,9))-> paths

print(paths)

rna_subs_f%>% 
  mutate(Match="Filtered") %>% 
  right_join(rna_subs) %>% 
  mutate_at(vars(Match), ~replace_na(., "unfiltered")) %>% 
  right_join(dna_subs) %>% 
  mutate_at(vars(Match), ~replace_na(.,"absent")) %>% 
  filter(FPKM>=10) %>% 
  inner_join(paths, by="Sample") %>% 
  mutate(Start= position-1) %>%  #Because same format as bed file, i.e. 0 based.
  mutate(Stop= position) %>% 
  mutate(cmd_args = paste0("depth -r ", chrom, ":", Start, "-", Stop, " ", full_Path )) %>% 
  mutate(depth = map(cmd_args, ~system2(command = "samtools", args = .x, stdout = T))) %>%
  unnest(depth) %>%
  separate(depth, into = c("chr", "nt", "depth"), extra = "merge", sep = "\t")-> substitution_depths

         
write.table(substitution_depths, "subs_depth.tsv", dec=".", col.names=TRUE, append=FALSE)
#Do not do write table next time... do write_tsv instead.