#!/usr/bin/env Rscript
# Exploring Genotype Output tables
# 23.02.23
# Mirjam Karlsson-Müller
# Genotype.R
#
#
# Description:
#



library(tidyverse)

#setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\Genotypes_Unfiltered")
#setwd("\\proj\\sens2019019\\nobackup\\wharf\\mimuka\\mimuka-sens2019019\\Outputs\\Genotype_Tables")

#read in chromosome genotype tables separately. to not overload things 

#Vlundur recommends 100 lines per file, then write code for those and then do it in interactive bianca session when r code is functioning. 
#Lets do that. tomorrow.

#Skipped Y and M. By recommendation of Helena.

genotype<- c(1:22, "X") %>%
  as_tibble() %>%
  #mutate(path = paste0("/proj/sens2019019/nobackup/wharf/mimuka/mimuka-sens2019019/Outputs/Filtered_Genotype_Tables/chr", value, "_genotype_table.tsv")) %>%
  mutate(path = paste0("chr", value, "_genotype_table.tsv")) %>%
  mutate(skip_n= ifelse((value =="1" | value=="2"), 9, 8 )) %>%
  #mutate(df = map2(path, skip_n, ~read_tsv(.x, show_col_types = F, skip = .y, n_max = 100)))
  mutate(df = map2(path, skip_n, ~read_tsv(.x, show_col_types = F, skip = .y)))

genotype %>% 
  unnest(df) %>% 
  select(!c(value, path, skip_n))-> genotype

#Now most of these are useless, as we can only do comparisons for locations where at least 2 genotypes are found. Especially reference is of interest. 
#So we make the same plot again but with only locations that have homozygous reference at at least one sample.
genotype %>% 
  mutate(hmza=rowSums(genotype=="1/1")) %>% 
  mutate(hetz=rowSums(genotype=="0/1")) %>% 
  mutate(hmzr=rowSums(genotype=="0/0")) %>% 
  #Filter rows where no other genotype but NE, ND is found
  filter((hmzr+hmza+hetz)!=0) %>% 
  #Remove entries without any hmzr. No thats skewing the results  
  #filter(hmzr!=0) %>% 
  mutate(percentage=(hmzr+hmza+hetz)/3455*100) %>% 
  mutate(bin=case_when(percentage<1 ~ "<1",
                       percentage>=1 & percentage<=10 ~ "1-10",
                       percentage>10 & percentage<=20 ~ "11-20",
                       percentage>20 & percentage<=30 ~ "21-30",
                       percentage>30 & percentage<=40 ~ "31-40",
                       percentage>40 & percentage<=50 ~ "41-50",
                       percentage>50 & percentage<=60 ~ "51-60",
                       percentage>60 & percentage<=70 ~ "61-70",
                       percentage>70 & percentage<=80 ~ "71-80",
                       percentage>80 & percentage<=90 ~ "81-90",
                       percentage>90 & percentage<=100 ~ "91-100")) %>% 
  mutate(bin=factor(bin,levels = c("<1","1-10", "11-20", "21-30", "31-40","41-50","51-60","61-70","71-80","81-90","91-100", ">100"))) %>% 
  select(!c(hmzr, hmza, hetz)) %>% 
  pivot_longer(cols=!c(Location, percentage, bin), names_to="Samples", values_to="Genotype") %>% 
  filter((Genotype=="0/0" | Genotype=="0/1" | Genotype=="1/1"))-> histo_data

ggplot(histo_data, aes(fill=factor(Genotype, levels=c("1/1", "0/1", "0/0")), x=bin))+
  xlab("Location found in percentage of samples")+
  ylab("Number of genotype per type")+
  ggtitle("Which genotypes are identified for locations found in percentage of samples")+
  geom_bar(position="stack")+
  theme_classic()+
  scale_color_brewer(palette = "BuPu")+
  scale_fill_brewer(name = "Genotype", labels = c("1/1", "0/1", "0/0"))+
  theme(text = element_text(size = 20))+
  #geom_text(aes(label = ..count..), stat = "count", position="stack", vjust =1.5 )+
  scale_x_discrete(guide=guide_axis(n.dodge=2))


ggsave("Genotype.png",
       device = "png",
       width = 175,
       height = 200,
       units = "mm")

rm(histo_data)

#Do same plots but only with variants that have a common dnsnp id.

dbsnp_vcf<-read_tsv("..\\Database\\no_header.vcf", name_repair="universal")

#str(dbsnp_vcf)

dbsnp_vcf %>%
  filter(nchar(ALT)==1 & nchar(REF)==1) %>%
  mutate(CHROM=.CHROM) %>%
  select(!.CHROM) %>% 
  mutate(Location=paste(paste("chr",CHROM, sep=""), as.character(POS), REF, paste("(", ALT, ")", sep=""), sep="_")) %>% 
  select(c(Location, ID))-> dbsnp

rm(dbsnp_vcf)

#head(dbsnp)
#head(genotype$Location)

genotype %>% 
  inner_join(dbsnp, by="Location")->variants_dbsnp


#Make same plot
variants_dbsnp %>% 
  mutate(hmza=rowSums(variants_dbsnp=="1/1")) %>% 
  mutate(hetz=rowSums(variants_dbsnp=="0/1")) %>% 
  mutate(hmzr=rowSums(variants_dbsnp=="0/0")) %>% 
  #Filter rows where no other genotype but NE, ND is found
  filter((hmzr+hmza+hetz)!=0) %>% 
  #Remove entries without any hmzr
  #filter(hmzr!=0) %>% 
  mutate(percentage=(hmzr+hmza+hetz)/3455*100) %>% 
  mutate(bin=case_when(percentage<1 ~ "<1",
                       percentage>=1 & percentage<=10 ~ "1-10",
                       percentage>10 & percentage<=20 ~ "11-20",
                       percentage>20 & percentage<=30 ~ "21-30",
                       percentage>30 & percentage<=40 ~ "31-40",
                       percentage>40 & percentage<=50 ~ "41-50",
                       percentage>50 & percentage<=60 ~ "51-60",
                       percentage>60 & percentage<=70 ~ "61-70",
                       percentage>70 & percentage<=80 ~ "71-80",
                       percentage>80 & percentage<=90 ~ "81-90",
                       percentage>90 & percentage<=100 ~ "91-100")) %>% 
  mutate(bin=factor(bin,levels = c("<1","1-10", "11-20", "21-30", "31-40","41-50","51-60","61-70","71-80","81-90","91-100", ">100"))) %>% 
  select(!c(hmzr, hmza, hetz)) %>% 
  pivot_longer(cols=!c(Location, percentage, bin, ID), names_to="Samples", values_to="Genotype") %>% 
  filter((Genotype=="0/0" | Genotype=="0/1" | Genotype=="1/1"))-> histo_data

ggplot(histo_data, aes(fill=factor(Genotype, levels=c("1/1", "0/1", "0/0")), x=bin))+
  xlab("Location found in percentage of samples")+
  ylab("Number of genotype per type")+
  scale_fill_discrete(name = "Genotype", labels = c("1/1", "0/1", "0/0"))+
  ggtitle("Which genotypes are identified for locations found in percentage of samples")+
  geom_bar(position="stack")+
  #geom_text(aes(label = ..count..), stat = "count", position="stack", vjust =1.5 )+
  theme_classic()+
  scale_color_brewer(palette = "BuPu")+
  scale_fill_brewer(name = "Genotype", labels = c("1/1", "0/1", "0/0"))+
  theme(text = element_text(size = 20))+
  scale_x_discrete(guide=guide_axis(n.dodge=2))

ggsave("dbsnp_Genotype.png",
       device = "png",
       width = 175,
       height = 200,
       units = "mm")

#Summary statistics: Where are the variants located? In exons? Introns?
#For both unfiltered and filtered genotypes.



