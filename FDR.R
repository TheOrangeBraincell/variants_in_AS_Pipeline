# 29.03.23
# FDR correction for KW outputs
# Mirjam Karlsson-Müller
# 
# Description: Takes a tsv file as input, where one column is p.value. In case of this pipeline,
# the output of Statistics.R is used as an input for this script. Also requires numbers of tests to correct for.
# Will return a new tsv file with an additional column q.values which contains corrected p values.
# 
# Useage: 
#   Rscript FDR.R 18677433 ESR1_KW.tsv
# 

library(tidyverse)

args=commandArgs(trailingOnly=TRUE)

number_of_tests<-as.integer(args[1])
file_to_correct<- args[2]

#Read in file
KW_scores<-read_tsv(file_to_correct, show_col_types=FALSE)

#Create new column with fdr corrected p values
KW_scores$q.values=p.adjust(KW_scores$p.value, method="fdr", n=number_of_tests)

#Write new output file
KW_scores %>% 
  arrange(q.values) %>%
  filter(q.values<0.05) %>% 
  write_tsv(
  paste0(substr(file_to_correct, 1, nchar(file_to_correct)-4), "_fdr.tsv"),
  na = "NA",
  quote = "none",
  escape = c("double", "backslash", "none"),
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress())