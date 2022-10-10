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

dna_variants<- read_tsv("tnbc_variants_fpkm.tsv")
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

rm(rna_variants_filtered)

rna_variants_unfiltered %>% 
  filter(Sample %in% dna_samples)->rna_vu
rm(rna_variants_unfiltered)

dna_variants %>% 
  filter(Sample %in% rna_samples)->dna_vu #dna variants updated

rm(dna_variants)

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
  filter(dbsnp=="No") %>% 
  select(!dbsnp)

rnaf_exons_no_germ #117 624

common_wo_germline_exonf<-inner_join(dna_exons, rnaf_exons_no_germ)
common_wo_germline_exonf
#4417
#So most matches are not germline variants. which makes sense given that
# germline variants are excluded from the dna data set (mostly.)

rna_exons_no_germ<-rna_exons %>% 
  filter(dbsnp=="No")

rna_exons_no_germ #2 756 494

common_wo_germline_exon<-inner_join(dna_exons, rna_exons_no_germ)
common_wo_germline_exon


#Relationship to Expression level?
#using filtered data first
rnaf_exons_no_germ %>% 
  mutate(Match="T")->rnaf_exons_no_germ

dna_exons %>% 
  left_join(rnaf_exons_no_germ) %>%
  mutate(FPKM=as.numeric(FPKM)) %>% 
  mutate_at(vars(Match), ~replace_na(., "F")) %>% 
  drop_na(FPKM)->dna_exons_rnaf

dna_exons_rnaf

dna_exons_rnaf$Match<-as.factor(dna_exons_rnaf$Match)
#maximum_FPKM<-max(dna_exons$FPKM)
dna_exons_rnaf$bins<-cut(dna_exons_rnaf$FPKM, breaks=seq(from=0, to=ceiling(max(dna_exons_rnaf$FPKM)), by=500))

#dna_exons %>% 
#  mutate(bins=cut(FPKM, breaks=seq(from=0, to=ceiling(max(FPKM)), by=500)))

#Not percentages.
dna_exons_rnaf %>% 
  drop_na(FPKM) %>%
  drop_na(bins) %>% 
  ggplot(aes(fill=Match, x=bins))+
  xlab("ranges of FPKM values")+
  ylab("Number of variants")+
  scale_fill_discrete(name = "Variant found in", labels = c("DNA", "DNA + RNA"))+
  ggtitle("Variants found in DNA vs variants found in both DNA and RNA")+
  #xlim(-1,501)+
  geom_bar()


#Only look at counts for fpkm <500

dna_exons_rnaf_filtered<- dna_exons_rnaf %>%
  select(!bins) %>% 
  filter(FPKM<=500)

dna_exons_rnaf_filtered$bins<-cut(dna_exons_rnaf_filtered$FPKM, breaks=seq(from=0, to=500, by=50))
#dna_exons_filtered
dna_exons_rnaf_filtered %>% 
  drop_na(bins)->dna_exons_rnaf_filtered_no_na

ggplot(dna_exons_rnaf_filtered_no_na, aes(fill=Match, x=bins))+
  xlab("FPKM values")+
  ylab("Number of variants")+
  scale_fill_discrete(name = "Variant found in", labels = c("DNA", "DNA + RNA"))+
  ggtitle("Variants found in DNA vs variants found in both DNA and RNA (filtered)")+
  geom_bar(position="fill")+
  geom_text(aes(label = ..count..), stat = "count", position="fill", vjust =1.5 )+
  scale_x_discrete(guide=guide_axis(n.dodge=2))
  

#For unfiltered rna data.
dna_exons<- dna_vu %>% 
  filter(Location=="Exon")

rna_exons_no_germ %>% 
  mutate(Match="T")->rna_exons_no_germ

dna_exons %>% 
  left_join(rna_exons_no_germ) %>%
  mutate(FPKM=as.numeric(FPKM)) %>% 
  mutate_at(vars(Match), ~replace_na(., "F")) %>% 
  drop_na(FPKM)->dna_exons

dna_exons

dna_exons$Match<-as.factor(dna_exons$Match)
#maximum_FPKM<-max(dna_exons$FPKM)
dna_exons$bins<-cut(dna_exons$FPKM, breaks=seq(from=0, to=ceiling(max(dna_exons$FPKM)), by=500))

#dna_exons %>% 
#  mutate(bins=cut(FPKM, breaks=seq(from=0, to=ceiling(max(FPKM)), by=500)))

#Not percentages.
dna_exons %>% 
  drop_na(FPKM) %>%
  drop_na(bins) %>% 
  ggplot(aes(fill=Match, x=bins))+
  xlab("ranges of FPKM values")+
  ylab("Number of variants")+
  scale_fill_discrete(name = "Variant found in", labels = c("DNA", "DNA + RNA"))+
  ggtitle("Variants found in DNA vs variants found in both DNA and RNA")+
  #xlim(-1,501)+
  geom_bar()


#Only look at counts for fpkm <500

dna_exons_filtered<- dna_exons %>%
  select(!bins) %>% 
  filter(FPKM<=500)

dna_exons_filtered$bins<-cut(dna_exons_filtered$FPKM, breaks=seq(from=0, to=500, by=50))
#dna_exons_filtered
dna_exons_filtered %>% 
  drop_na(bins)->dna_exons_filtered_no_na

ggplot(dna_exons_filtered_no_na, aes(fill=Match, x=bins))+
  xlab("FPKM values")+
  ylab("Number of variants")+
  scale_fill_discrete(name = "Variant found in", labels = c("DNA", "DNA + RNA"))+
  ggtitle("Variants found in DNA vs variants found in both DNA and RNA (unfiltered)")+
  geom_bar(position="fill")+
  geom_text(aes(label = ..count..), stat = "count", position="fill", vjust =1.5 )+
  scale_x_discrete(guide=guide_axis(n.dodge=2))
#remove the position="fill" to not have percentage on y axis but counts.
