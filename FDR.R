# 29.03.23
# FDR correction for KW outputs
# Mirjam Karlsson-Müller

library(tidyverse)

args=commandArgs(trailingOnly=TRUE)

number_of_tests<-as.integer(args[1])
file_to_correct<- args[2]

KW_scores<-read_tsv(file_to_correct, show_col_types=FALSE)


KW_scores$q.values=p.adjust(KW_scores$p.value, method="fdr", n=number_of_tests)

KW_scores %>% 
  arrange(q.values) %>% 
  write_tsv(
  paste0(substr(file_to_correct, 1, nchar(file_to_correct)-4), "_fdr.tsv"),
  na = "NA",
  quote = "none",
  escape = c("double", "backslash", "none"),
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress())