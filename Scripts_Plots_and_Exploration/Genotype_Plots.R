# Genotype Tables Plot
# 07-04-23
# Mirjam Karlsson-Müller


library(tidyverse)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\Genotype_Plots")
#Read in all the count tables

genotypes<- read_tsv("gene_list.txt", col_names="Gene") %>% 
  mutate(path=paste0("Genotype_Counts/", Gene, "_GT_Counts.tsv")) %>% 
  mutate(df = map(path, ~read_tsv(.x, show_col_types=F, col_types = cols("0/0" = col_character(),"0/X" = col_character(),"X/X" = col_character()))))

 genotypes%>% 
  unnest(c(df)) %>% 
  select(!path) ->GT_DF

 #I gave it invalid names in..... so ill redo column names.
colnames(GT_DF)<- c("Gene", "Location", "A", "B", "C")
GT_DF %>% 
  mutate(total=as.integer(A)+as.integer(B)+as.integer(C)) %>% 
  mutate(hmzr=as.integer(A)/total) %>% 
  mutate(hetz=as.integer(B)/total) %>% 
  mutate(hmza=as.integer(C)/total) %>% 
  select(!c(A, B, C, Location, total)) %>% 
  arrange(desc(hmzr, hetz, hmza)) %>% 
  mutate(n=seq.int(nrow(.))) %>% 
  pivot_longer(cols=!c(Gene, n), names_to = "Genotype", values_to = "Fraction")-> GT_fractions

#Subsample cause too many points.

uniq_locations <- GT_fractions %>% 
  select(n) %>% 
  distinct()

subsample<-uniq_locations[sample(nrow(uniq_locations), 1500),]
subsample %>% 
  arrange(n) %>% 
  mutate(m=seq.int(nrow(.))) %>% 
  left_join(GT_fractions, by=c("n")) -> subset_fractions


ggplot(subset_fractions, aes(x=m, y=Fraction, group=Genotype, color=Genotype))+
  geom_point(size=2)+
  theme_classic()+
  xlab("Random subset of variants sorted by frequency of 0/0")+
  theme(text = element_text(size = 20), axis.text=element_text(size=20))+
  ylab("Fractions of genotype")+
  labs(color = "Genotype") +
  scale_color_brewer(labels = c("1/1", "0/1", "0/0"), palette ="Accent")
 