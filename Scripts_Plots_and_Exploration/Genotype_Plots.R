# Genotype Tables Plot
# 07-04-23
# Mirjam Karlsson-Müller
#
# Description: Plots output of Genotyping tables. Prepped with genotype_Counts.py, 
# whose output files we are reading in here. Shows the fractions of genotypes determined.
# 



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
  theme(text = element_text(size = 25), axis.text=element_text(size=25))+
  ylab("Fractions of genotypes assigned")+
  labs(color = "Genotype") +
  scale_color_brewer(labels = c("1/1", "0/1", "0/0"), palette ="Accent")
 

#Lets try another. like the OG of genotype plots I had made. (That showed us it was wrong ._.)

GT_DF %>% 
  mutate(total=as.integer(A)+as.integer(B)+as.integer(C)) %>% 
  mutate(percentage=total/3455 *100) %>% 
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
                       percentage>90 & percentage<=100 ~ "91-100")) -> binned_data
  #The fraction would have to be over the total amount in that location bin...
binned_data %>% 
  select(!c(total, percentage)) %>%
  pivot_longer(cols=!c(Gene, Location, bin), names_to="Genotype", values_to="Counts") %>% 
  select(!c(Gene, Location)) %>% 
  group_by(bin, Genotype) %>% 
  mutate(n=as.integer(Counts)) %>% 
  select(!Counts) %>% 
  summarise(sum_counts=sum(n)) %>% 
  pivot_wider(names_from=Genotype, values_from=sum_counts) %>% 
  mutate(total=A+B+C) %>% 
  mutate(hmzr=A/total) %>% 
  mutate(hetz=B/total) %>% 
  mutate(hmza=C/total) %>% 
  select(!c(A, B, C, total)) %>% 
  pivot_longer(cols=!c(bin), names_to="Genotype", values_to="Fraction")-> GT_histo


GT_histo


ggplot(GT_histo, aes(fill=factor(Genotype, levels=c("hmzr", "hetz", "hmza")), y=Fraction, x=bin)) + 
  theme_classic()+
  geom_bar(position="stack", stat="identity") +
  scale_color_brewer(palette = "BuPu")+
  scale_fill_brewer(name = "Genotype", labels = c("0/0", "0/1", "1/1")) +
  xlab("Location found in percentage of samples")+
  ylab("Fractions of genotypes assigned")+
  theme(text = element_text(size = 25), axis.text=element_text(size=25))+
  #geom_text(aes(label = ..count..), stat = "count", position="stack", vjust =1.5 )+
  scale_x_discrete(guide=guide_axis(n.dodge=2))


