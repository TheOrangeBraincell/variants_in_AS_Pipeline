# tnbc_comparison.R
# Mirjam Karlsson-Müller
# 29.09.22
#
# Description: Comparison of variants found based on DNA samples of 252 patients,
# compared to variants found on RNA samples of the same 252 patients.
#
#

library(tidyverse)
library(RColorBrewer)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\TNBC_comparison")

dna_variants<- read_tsv("tnbc_variants_fpkm.tsv")
rna_variants_filtered<-read_tsv("rna_variants_filtered.txt")
rna_variants_unfiltered<-read_tsv("rna_variants.txt")




# How many of the variants in the dna files did we find with each versin of the rna variants?
#dna_variants
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

rm(rna_variants_filtered)

rna_variants_unfiltered %>% 
  filter(Sample %in% dna_samples)->rna_vu
rm(rna_variants_unfiltered)

dna_variants %>% 
  filter(Sample %in% rna_samples)->dna_vu #dna variants updated

rm(dna_variants)

#-----------------------------SUBSTITUTIONS--------------------------------------
# so exclude all lines with alt or ref longer than 1 char.
rna_vfu

rna_vfu %>% 
  filter(nchar(alt)==1 & nchar(ref)==1)->rna_vfu

rna_vu %>% 
  filter(nchar(alt)==1 & nchar(ref)==1)->rna_vu

dna_vu %>% 
  filter(nchar(alt)==1 & nchar(ref)==1)->dna_vu

#Final Plot attempt--------------------
#Taking helenas comments into consideration, we redo this into one plot!
#library(forcats)
dna_exons<- dna_vu %>% 
  filter(Location=="Exon")

rnaf_exons<-rna_vfu %>% 
  filter(Location=="Exon")

rna_exons<-rna_vu %>% 
  filter(Location=="Exon")

#No exclusion of germline variants this time!

#smaller bins
#dna_exons$bins<-cut(dna_exons$FPKM, breaks=seq(from=0, to=100, by=10), include.lowest = TRUE,)
#replace bins=NA with >100
dna_exons

dna_exons %>%
  mutate(FPKM=round(FPKM, digits=0)) %>% 
  mutate(bins=case_when(FPKM>=0 & FPKM<=1 ~ "<1",
                        FPKM>1 & FPKM<=10 ~ "1-10",
                        FPKM>10 & FPKM<=20 ~ "11-20",
                        FPKM>20 & FPKM<=30 ~ "21-30",
                        FPKM>30 & FPKM<=40 ~ "31-40",
                        FPKM>40 & FPKM<=50 ~ "41-50",
                        FPKM>50 & FPKM<=60 ~ "51-60",
                        FPKM>60 & FPKM<=70 ~ "61-70",
                        FPKM>70 & FPKM<=80 ~ "71-80",
                        FPKM>80 & FPKM<=90 ~ "81-90",
                        FPKM>90 & FPKM<=100 ~ "91-100",
                        FPKM>100 ~ ">100"))  %>% 
  #mutate_at("bins", ~replace_na(., ">100")) %>%
  #select(bins) %>% 
  #arrange() %>% 
  #distinct()
  mutate(bins=as.factor(bins))->dna_exons

#Add both filtered and unfiltered with different Match columns.
rnaf_exons %>%
  mutate(Match="filtered")->rnaf_exons

rna_exons %>% 
  left_join(rnaf_exons) %>%
  mutate_at(vars(Match), ~replace_na(., "unfiltered")) %>% 
  mutate(FPKM=round(FPKM, digits=0))->rna_exons

dna_exons %>% 
  left_join(rna_exons) %>% 
  mutate_at(vars(Match), ~replace_na(.,"absent")) %>% 
  mutate(FPKM=as.numeric(FPKM)) %>% 
  drop_na(FPKM) %>% 
  mutate(bins=factor(bins,levels = c("<1","1-10", "11-20", "21-30", "31-40","41-50","51-60","61-70","71-80","81-90","91-100", ">100")))-> dna_exons

dna_exons %>% 
  filter(is.na(bins))

ggplot(dna_exons, aes(fill=factor(Match, levels=c("absent", "unfiltered", "filtered")), x=bins))+
  theme_classic()+
  #scale_fill_viridis(discrete = T, name = "Fraction of Variants found in", labels = c("DNA", "DNA + RNA unfiltered", "DNA + RNA filtered"))+
  xlab("FPKM values")+
  ylab("%Variants found")+
  scale_color_brewer(palette = "BuPu")+
  scale_fill_brewer(name = "Variants found in", labels = c("DNA", "DNA + RNA unfiltered", "DNA + RNA filtered"))+
  ggtitle("Variants found in DNA data, that are also found in RNA data")+
  geom_bar(position="fill")+
  theme(text = element_text(size = 15), axis.text=element_text(size=15))+
  #geom_text(aes(label = ..count..), size=5, stat = "count", position="fill", vjust =1.5)+
  scale_x_discrete(guide=guide_axis(n.dodge=2))

# Another. But for bins from 1-10
dna_exons<- dna_vu %>% 
  filter(Location=="Exon")

dna_exons %>%
  mutate(FPKM=round(FPKM, digits=0)) %>% 
  mutate(bins=case_when(FPKM>=0 & FPKM<=1 ~ "<1",
                        FPKM>1 & FPKM<=2 ~ "1-2",
                        FPKM>2 & FPKM<=3 ~ "2-3",
                        FPKM>3 & FPKM<=4 ~ "3-4",
                        FPKM>4 & FPKM<=5 ~ "4-5",
                        FPKM>5 & FPKM<=6 ~ "5-6",
                        FPKM>6 & FPKM<=7 ~ "6-7",
                        FPKM>7 & FPKM<=8 ~ "7-8",
                        FPKM>8 & FPKM<=9 ~ "8-9",
                        FPKM>9 & FPKM<=10 ~ "9-10",
                        FPKM>10 ~ ">10"))  %>% 
  #mutate_at("bins", ~replace_na(., ">100")) %>%
  #select(bins) %>% 
  #arrange() %>% 
  #distinct()
  mutate(bins=as.factor(bins))->dna_exons

#Add both filtered and unfiltered with different Match columns.
rnaf_exons %>%
  mutate(Match="filtered")->rnaf_exons

rna_exons %>% 
  left_join(rnaf_exons) %>%
  mutate_at(vars(Match), ~replace_na(., "unfiltered")) %>% 
  mutate(FPKM=round(FPKM, digits=0))->rna_exons

dna_exons %>% 
  left_join(rna_exons) %>% 
  mutate_at(vars(Match), ~replace_na(.,"absent")) %>% 
  mutate(FPKM=as.numeric(FPKM)) %>% 
  drop_na(FPKM) %>% 
  mutate(bins=factor(bins,levels = c("<1","1-2", "2-3", "3-4", "4-5","5-6","6-7","7-8","8-9","9-10",">10")))-> dna_exons

dna_exons %>% 
  filter(is.na(bins))

ggplot(dna_exons, aes(fill=factor(Match, levels=c("absent", "unfiltered", "filtered")), x=bins))+
  theme_classic()+
  xlab("FPKM values")+
  ylab("%Variants found")+
  scale_color_brewer(palette = "BuPu")+
  scale_fill_brewer(name = "Variants found in", labels = c("DNA", "DNA + RNA unfiltered", "DNA + RNA filtered"))+
  ggtitle("Variants found in DNA data, that are also found in RNA data")+
  geom_bar(position="fill")+
  theme(text = element_text(size = 15), axis.text=element_text(size=15))+
  #geom_text(aes(label = ..count..), stat = "count", position="fill", vjust =1.5 )+
  scale_x_discrete(guide=guide_axis(n.dodge=2))
#General numbers of overlap!

#How many germline variants in dna_vu?
rna_vfu %>% 
  filter(dbsnp=="Yes") %>% 
  mutate(Match="Filtered") %>% 
  right_join(dna_vu) %>% 
  group_by(Match) %>% 
  count()

#Make a table with overlap numbers

rna_vfu%>% 
  mutate(Match="Filtered") %>% 
  right_join(rna_vu) %>% 
  mutate_at(vars(Match), ~replace_na(., "unfiltered")) %>% 
  #filter(dbsnp=="No") %>%
  left_join(dna_vu) %>% 
  mutate_at(vars(Match), ~replace_na(.,"absent")) %>% 
  #filter(Location=="Intergenic") %>% 
  group_by(Match) %>% 
  count()

#unique rna variants:
dna_vu %>% 
  mutate(Match="True") %>% 
  right_join(rna_vfu) %>% 
  mutate_at(vars(Match), ~replace_na(.,"False")) %>%
  filter(dbsnp=="No") %>%
  #filter(Location=="Intergenic") %>%
  group_by(Match) %>% 
  count()


# Filter Criteria Plot ------------------------------

library(tidyverse)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\TNBC_comparison")

dna_variants<- read_tsv("tnbc_variants_fpkm.tsv")
criteria<-read.csv("filter_criteria_vcf.tsv", header=TRUE, sep="\t")
#criteria<-read_tsv("filter_criteria_vcf.tsv")

#Remove samples not present in dna
criteria %>% 
  select(Sample) %>% 
  distinct() %>% 
  pull(Sample)-> rna_samples

dna_variants %>% 
  select(Sample) %>% 
  distinct() %>% 
  pull(Sample)->dna_samples


dna_variants %>% 
  filter(Sample %in% rna_samples)->dna_vu #dna variants updated

rm(dna_variants)

criteria %>% 
  filter(Sample %in% dna_samples)->rna_criteria 

rm(criteria)

#Format genotype
rna_criteria$genotype<-str_replace(rna_criteria$genotype, "1/0", "0/1")

#Substitutions only
rna_criteria %>% 
  filter(nchar(alt)==1 & nchar(ref)==1)->rna_criteria

head(rna_criteria)

#Prepare data for plotting.
rna_criteria %>%
  #mutate(Match="True") %>% 
  inner_join(dna_vu) %>%
  select(c(Sample, chrom, position, HMPOL, GC, VD, low_complexity,
           quality, ucsc_rep, MSI)) %>% 
  pivot_longer(cols=!c(Sample, chrom, position),
               names_to="Filter", values_to="Pass_Fail") ->dna_rna


ggplot(dna_rna, aes(fill=factor(Pass_Fail, levels=c("Pass", "Fail")), x=Filter))+
  theme_classic()+
  #scale_color_brewer(palette = "Pastel1")+
  xlab("Filter criteria to Fail")+
  ylab("Percentage of matching variants excluded by filter")+
  scale_fill_brewer(name = "Variants which", labels = c("Pass", "Fail"), palette="BuPu")+
  ggtitle("Matching variants between DNA and RNA data removed by filtering")+
  geom_bar(position="fill")+
  theme(text = element_text(size = 20), axis.text=element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(aes(label = ..count..), stat = "count", position="fill", vjust =1.5 )+
  scale_x_discrete(labels=c("GC-content>78%","HMPOL>5","low-complexity flag", "MSI>6", "quality score<55", "ucsc-repetitive flag", "Variant Depth <5"))

#Lets try and do arthurs plot.

dna_rna %>% 
  filter(Pass_Fail=="Fail") %>% 
  group_by(Sample, Filter, Pass_Fail) %>% 
  count()-> dna_rna_count

ggplot(dna_rna_count, aes(x=Filter, y=n ))+
  geom_boxplot()

dna_rna_count %>% 
  filter(n>=50)
#Its not really giving more information than the one we already made.


#How could we make a plot that includes matches, filter conditions and FPKM?

# a scatterplot- FPKM vs Filters, colors for Pass and Fail. Only with matches.
# Lets try.

rna_criteria %>%
  inner_join(dna_vu) %>%
  select(c(Sample, chrom, position, FPKM, HMPOL, GC, VD, low_complexity,
           quality, ucsc_rep, MSI)) %>% 
  pivot_longer(cols=!c(Sample, chrom, position, FPKM),
               names_to="Filter", values_to="Pass_Fail") %>% 
  drop_na(FPKM)-> dna_rna_fpkm
  #filter(FPKM<=100) %>% 
  #filter(FPKM>1)-> dna_rna_fpkm

ggplot(dna_rna_fpkm, aes(x=Filter, y=FPKM))+
  geom_boxplot(aes(fill=Pass_Fail))+xlab("Filter criteria to Fail")+
  scale_fill_discrete(name = "Fraction of Variants which", labels = c("Fail", "Pass"))+
  ggtitle("Do the matching variants who fail the filter criteria have higher FPKM values?")+
  ylim(0,100)+
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust=1))+
  scale_x_discrete(labels=c("GC-content>78%","HMPOL>5","low-complexity flag", "MSI>6", "quality score<55", "ucsc-repetitive flag", "Variant Depth <5"))



#How many of the entries failing the quality score, also fail other filters?

rna_criteria %>% 
  filter(quality=="Fail") %>% 
  filter(MSI=="Pass", HMPOL=="Pass", GC=="Pass", VD=="Pass", low_complexity=="Pass",
         ucsc_rep=="Pass") %>% 
  inner_join(dna_vu) %>% 
  count()






#-----------------------------INDELS---------------------------------------------
# We want to do the fpkm bin plots as well as tables with the counts for indels instead of substitutions.
# Just me trying to figure out why these things give different results.
# rna_vu
# 
# rna_vu %>% 
#   mutate(Type=case_when(
#     nchar(alt)==1 | nchar(ref)==1 ~"substitution",
#     nchar(alt)>1 | nchar(ref)>1 ~ "indel",
#     T ~ "other")  ) %>% 
#   group_by(Type) %>% 
#   count()
# 
# 
# rna_vu %>% 
#   mutate(Type=case_when(nchar(alt)==1 | nchar(ref)==1 ~"substitution",
#                         (nchar(alt)>1) | (nchar(ref)>1) ~ "indel")) %>%
#   filter(Type == "indel")
# 
# rna_vu %>% 
#   ungroup %>% 
#   filter(nchar(alt)>1 | nchar(ref)>1)->rna_indels
# 
# rna_indels
# 
# 
# rna_vu %>% 
#   mutate(across(c(ref, alt), ~nchar(.x), .names = "{.col}_n")) %>% 
#   filter(ref_n > 1 | alt_n > 1)
# 
# a
# 






# Create indel tables
rna_vfu %>% 
  filter(nchar(alt)>1 | nchar(ref)>1)->rna_vfu

rna_vu %>% 
  filter(nchar(alt)>1 | nchar(ref)>1)->rna_indels

rna_indels


#indels is a different table.
dna_indels <- read_tsv("tnbc_indels_fpkm.tsv")

dna_indels %>% 
  filter(Sample %in% rna_samples)->indels_vu

indels_exons<- indels_vu %>% 
  filter(Location=="Exon")

rnaf_exons<-rna_vfu %>% 
  filter(Location=="Exon")

rna_exons<-rna_vu %>% 
  filter(Location=="Exon")

indels_exons %>%
  mutate(FPKM=round(FPKM, digits=0)) %>% 
  mutate(bins=case_when(FPKM>=0 & FPKM<=1 ~ "<1",
                        FPKM>1 & FPKM<=10 ~ "1-10",
                        FPKM>10 & FPKM<=20 ~ "11-20",
                        FPKM>20 & FPKM<=30 ~ "21-30",
                        FPKM>30 & FPKM<=40 ~ "31-40",
                        FPKM>40 & FPKM<=50 ~ "41-50",
                        FPKM>50 & FPKM<=60 ~ "51-60",
                        FPKM>60 & FPKM<=70 ~ "61-70",
                        FPKM>70 & FPKM<=80 ~ "71-80",
                        FPKM>80 & FPKM<=90 ~ "81-90",
                        FPKM>90 & FPKM<=100 ~ "91-100",
                        FPKM>100 ~ ">100"))  %>% 
  #mutate_at("bins", ~replace_na(., ">100")) %>%
  #select(bins) %>% 
  #arrange() %>% 
  #distinct()
  mutate(bins=as.factor(bins))->binned_indels_exons

#Add both filtered and unfiltered with different Match columns.
rnaf_exons %>%
  mutate(Match="filtered")->rnaf_exons


rna_exons %>% 
  left_join(rnaf_exons) %>%
  mutate_at(vars(Match), ~replace_na(., "unfiltered")) %>% 
  mutate(FPKM=round(FPKM, digits=0))->rna_exons

binned_indels_exons %>% 
  left_join(rna_exons) %>% 
  mutate_at(vars(Match), ~replace_na(.,"absent")) %>% 
  mutate(FPKM=as.numeric(FPKM)) %>% 
  drop_na(FPKM) %>% 
  mutate(bins=factor(bins,levels = c("<1","1-10", "11-20", "21-30", "31-40","41-50","51-60","61-70","71-80","81-90","91-100", ">100")))-> binned_indels_exons


ggplot(binned_indels_exons, aes(fill=factor(Match, levels=c("absent", "unfiltered", "filtered")), x=bins))+
  theme_classic()+
  xlab("FPKM values")+
  ylab("Fraction of Variants found")+
  scale_fill_brewer(name = "Fraction of Variants found in", labels = c("DNA", "DNA + RNA unfiltered", "DNA + RNA filtered"), palette="BuPu")+
  ggtitle("Fraction of variants found in DNA data, that are also found in RNA data")+
  geom_bar(position="fill")+
  geom_text(aes(label = ..count..), stat = "count", position="fill", vjust =1.5 )+
  scale_x_discrete(guide=guide_axis(n.dodge=2))

# To make the tables.
rna_vfu%>% 
  mutate(Match="Filtered") %>% 
  right_join(rna_vu) %>% 
  mutate_at(vars(Match), ~replace_na(., "unfiltered")) %>% 
  filter(dbsnp=="No") %>%
  right_join(indels_vu) %>% 
  mutate_at(vars(Match), ~replace_na(.,"absent")) %>% 
  filter(Location=="Intergenic") %>% 
  group_by(Match) %>% 
  count()

#unique rna variants:
indels_vu %>% 
  mutate(Match="True") %>% 
  right_join(rna_vfu) %>% 
  mutate_at(vars(Match), ~replace_na(.,"False")) %>%
  filter(dbsnp=="No") %>%
  filter(Location=="Intergenic") %>%
  group_by(Match) %>% 
  count()


