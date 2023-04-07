# Plot PSI Scores
# 07-04-23  
# Mirjam Karlsson-Müller
#
# Create Overview Plots of the PSI Score distribution per AS Type. 
#

library(tidyverse)


setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\PSI_Results")

#Read in all count tables
PSI<- read_tsv("gene_list.txt", col_names="Gene") %>% 
  mutate(path=paste0("PSI_Counts/PSI_Counts/", Gene, "_Counts.tsv")) %>% 
  mutate(df = map(path, ~read_tsv(.x, show_col_types=F, col_types = cols(n = col_character()))))


PSI %>% 
  unnest(c(df)) %>% 
  select(!path) %>% 
  mutate(n=as.integer(n))-> PSI_DF

PSI_DF %>%
  filter(event=="IR") %>% 
ggplot(aes(x=bin, y=n))+
theme_classic()+
geom_col(fill="darkgoldenrod")+
xlab("PSI score interval")+
  ylab("Number of scores")+
  theme(text = element_text(size = 20), axis.text=element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))
  #geom_text(aes(label = ..count..), position="fill", vjust =1.5 )


PSI_DF %>%
  filter(event=="AA") %>% 
  ggplot(aes(x=bin, y=n))+
  theme_classic()+
  geom_col(fill="darkred")+
  xlab("PSI score interval")+
  ylab("Number of scores")+
  theme(text = element_text(size = 20), axis.text=element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))
#geom_text(aes(label = ..count..), position="fill", vjust =1.5 )


PSI_DF %>%
  filter(event=="AD") %>% 
  ggplot(aes(x=bin, y=n))+
  theme_classic()+
  geom_col(fill="darkolivegreen4")+
  xlab("PSI score interval")+
  ylab("Number of scores")+
  theme(text = element_text(size = 20), axis.text=element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))
#geom_text(aes(label = ..count..), position="fill", vjust =1.5 )


PSI_DF %>%
  filter(event=="CE") %>% 
  ggplot(aes(x=bin, y=n))+
  theme_classic()+
  geom_col(fill="deepskyblue3")+
  xlab("PSI score interval")+
  ylab("Number of scores")+
  theme(text = element_text(size = 20), axis.text=element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))
#geom_text(aes(label = ..count..), position="fill", vjust =1.5 )
