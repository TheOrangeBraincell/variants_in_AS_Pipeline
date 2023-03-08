# FPKM_Plot.R
# Mirjam Karlsson-Müller
# 30.01.23
#
# Description: 
#
#

library(tidyverse)
setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs")

fpkm_data<-read_tsv("fpkm_table.tsv")
genes<-read_tsv("gene_ranges.tsv")

head(fpkm_data)
head(genes)
# Pivot data for plotting

fpkm_pivot<-fpkm_data %>% 
  pivot_longer(cols=!c(Location), names_to="Samples", values_to="FPKM")


#How about a plot of counts instead?
fpkm_pivot %>% 
  mutate(classifier=case_when(FPKM>10 ~ TRUE,
                        FPKM<=10 ~ FALSE)) %>% 
  group_by(Location) %>% 
  summarize(n=length(which(classifier))) %>% 
  mutate(percentage= (n/3455)*100) %>% 
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
  mutate(bin=factor(bin,levels = c("<1","1-10", "11-20", "21-30", "31-40","41-50","51-60","61-70","71-80","81-90","91-100", ">100")))-> histo_data

#Write into output file for filtering psi scores.
histo_data %>% 
  filter(bin!="<1") %>% 
  select("Location")->output



#Merge with gene ranges to find coordinates of these genes.
output%>% 
  left_join(genes, by=c("Location"="Gene")) %>% 
  group_by(Chrom) %>% 
  arrange(Chrom, start) %>% 
  drop_na()-> genes_of_interest

#Now we have a table that tells us what ranges to run the PSIs for.
write.table(genes_of_interest, "expressed_genes.tsv", append =FALSE, sep="\t", dec=".", col.names=TRUE, row.names=FALSE,quote=FALSE)

ggplot(histo_data, aes(x=bin))+
  geom_bar()+
  xlab("Percentage of samples with FPKM >10")+
  ylab("Number of genes")+
  ggtitle("Number of genes with an FPKM >10 in certain percentage of the samples")+
  geom_text(aes(label = ..count..), stat = "count", position="fill", vjust =1.5 )+
  scale_x_discrete(labels=c("0-1%","1-10%","10-20%","20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"), guide=guide_axis(n.dodge=2))
  
