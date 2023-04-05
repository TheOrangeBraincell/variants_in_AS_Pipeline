# Bin PSI counts per AS event
# Bin_PSI.R
# Mirjam Karlsson-Müller
# 04-04-23


library(tidyverse)

args=commandArgs(trailingOnly=TRUE)

PSI_file<-args[1]
Output<-args[2]
gene <- args[3]


psi<-read_tsv(PSI_file) %>%
  #mutate(Gene=gene) %>% 
  pivot_longer(cols=!c(Location), names_to="Sample", values_to="PSI") %>% 
  filter(PSI!="NaN") %>% 
  #(event=="CE") %>% 
  separate(col=Location, c("event", "A", "B", "C", "D"), sep="_") %>% 
  mutate(event= ifelse(startsWith(A, "A"), A, event)) %>% 
  select(event, PSI) %>%
  mutate(bin=case_when(PSI<0.1 ~ "[0,0.1)",
                       PSI>=0.1 & PSI<=0.2 ~ "[0.1,0.2)",
                       PSI>0.2 & PSI<=0.3 ~ "[0.2,0.3)",
                       PSI>0.3 & PSI<=0.4 ~ "[0.3,0.4)",
                       PSI>0.4 & PSI<=0.5 ~ "[0.4,0.5)",
                       PSI>0.5 & PSI<=0.6 ~ "[0.5,0.6)",
                       PSI>0.6 & PSI<=0.7 ~ "[0.6,0.7)",
                       PSI>0.7 & PSI<=0.8 ~ "[0.7,0.8)",
                       PSI>0.8 & PSI<=0.9 ~ "[0.8,0.9)",
                       PSI>0.9 & PSI<=1 ~ "[0.9-1]")) %>% 
  group_by(event, bin) %>% 
  count() %>% 
  write_tsv(Output, na = "NA",
            quote = "none",
            escape = c("double", "backslash", "none"),
            eol = "\n",
            num_threads = readr_threads(),
            progress = show_progress())