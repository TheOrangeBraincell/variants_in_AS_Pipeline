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

#setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\Download220223")
#setwd("\\proj\\sens2019019\\nobackup\\wharf\\mimuka\\mimuka-sens2019019\\Outputs\\Genotype_Tables")

#read in chromosome genotype tables separately. to not overload things 

#Vlundur recommends 100 lines per file, then write code for those and then do it in interactive bianca session when r code is functioning. 
#Lets do that. tomorrow.

#Skipped Y and M. By recommendation of Helena.

genotype<- c(1:22, "X") %>%
  as_tibble() %>%
  mutate(path = paste0("chr", value, "_genotype_table.tsv")) %>%
  mutate(skip_n= ifelse((value =="1" | value=="2"), 9, 8 )) %>%
  mutate(df = map2(path, skip_n, ~read_tsv(.x, show_col_types = F, skip = .y, n_max = 100)))
  #mutate(df = map2(path, skip_n, ~read_tsv(.x, show_col_types = F, skip = .y)))

genotype %>% 
  unnest(df) %>% 
  select(!c(value, path, skip_n))-> genotype

genotype %>% 
  mutate(hmza=rowSums(genotype=="1/1")) %>% 
  mutate(hetz=rowSums(genotype=="0/1")) %>% 
  mutate(hmzr=rowSums(genotype=="0/0")) %>% 
  #Filter rows where no other genotype but NE, ND is found
  filter((hmzr+hmza+hetz)!=0) %>% 
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
  #select(c(percentage, bin)) %>% 
  pivot_longer(cols=!c(Location, hmzr, hetz, hmza, percentage, bin), names_to="Samples", values_to="Genotype") %>% 
  filter((Genotype=="0/0" | Genotype=="0/1" | Genotype=="1/1"))-> histo_data

png("../all_genotypes.png")

ggplot(histo_data, aes(fill=factor(Genotype, levels=c("0/0", "0/1", "1/1")), x=bin))+
  xlab("Location found in percentage of samples")+
  ylab("Number of genotype per type")+
  scale_fill_discrete(name = "Genotype", labels = c("0/0", "0/1", "1/1"))+
  ggtitle("Which genotypes are identified for locations found in percentage of samples")+
  geom_bar(position="fill")+
  geom_text(aes(label = ..count..), stat = "count", position="fill", vjust =1.5 )+
  scale_x_discrete(guide=guide_axis(n.dodge=2))

dev.off()

#Now most of these are useless, as we can only do comparisons for locations where at least 2 genotypes are found. Especially reference is of interest. 
#So we make the same plot again but with only locations that have homozygous reference at at least one sample.
genotype %>% 
  mutate(hmza=rowSums(genotype=="1/1")) %>% 
  mutate(hetz=rowSums(genotype=="0/1")) %>% 
  mutate(hmzr=rowSums(genotype=="0/0")) %>% 
  #Filter rows where no other genotype but NE, ND is found
  filter((hmzr+hmza+hetz)!=0) %>% 
  #Remove entries without any hmzr
  filter(hmzr!=0) %>% 
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
  #select(c(percentage, bin)) %>% 
  pivot_longer(cols=!c(Location, hmzr, hetz, hmza, percentage, bin), names_to="Samples", values_to="Genotype") %>% 
  filter((Genotype=="0/0" | Genotype=="0/1" | Genotype=="1/1"))-> histo_data

png("../only_hmzr.png")
ggplot(histo_data, aes(fill=factor(Genotype, levels=c("1/1", "0/1", "0/0")), x=bin))+
  xlab("Location found in percentage of samples")+
  ylab("Number of genotype per type")+
  scale_fill_discrete(name = "Genotype", labels = c("1/1", "0/1", "0/0"))+
  ggtitle("Which genotypes are identified for locations found in percentage of samples")+
  geom_bar(position="stack")+
  geom_text(aes(label = ..count..), stat = "count", position="stack", vjust =1.5 )+
  scale_x_discrete(guide=guide_axis(n.dodge=2))

dev.off()

#How are the genotypes distributed over the samples? How many are common?









