#Statistics Script per gene
#27.03.23
#Mirjam Karlsson-Müller


library(tidyverse)

args=commandArgs(trailingOnly=TRUE)

PSI_file<-args[1]
genotype_file<-args[2]
gene<-args[3]


#prepped Genotype file has layout: Location, Gene, S1, S2 etc, skip 8 lines for header
#prepped PSI file has layout: event, Gene, S1, S2 etc. Skip 1 line for header.


genotype<- read_tsv(genotype_file, skip=8, show_col_types = FALSE)
PSI<-read_tsv(PSI_file, show_col_types=FALSE)

PSI %>% 
  mutate("event"=Location) %>% 
  select(!Location) %>% 
  pivot_longer(cols=!c(event, Gene), names_to="Sample", values_to="PSI") %>% 
  filter(PSI!="NAN")-> PSI_long

rm(PSI)

genotype %>% 
  mutate(variant=Location) %>% 
  select(!Location) %>% 
  pivot_longer(cols=!c(variant, Gene), names_to="Sample", values_to="Genotype") %>% 
  filter(Genotype!="NE") %>% 
  filter(Genotype!="ND") %>% 
  left_join(by=c("Gene", "Sample")) %>% 
  #Remove rows that either have no PSI or no genotype in that sample
  drop_na() -> comparison_data

rm(genotype, PSI_long)

#But now there might be only 1 distinct PSI value or genotype at a location! So there's no need to test
#So we need to find those and remove them.
comparison_data %>% 
  select(!c(PSI, Sample)) %>%
  distinct() %>% 
  group_by(Gene, variant, event) %>% 
  count() %>% 
  filter(n<2) %>%
  select(!n) %>% 
  left_join(comparison_data, by=c("Gene", "variant", "event")) %>% 
  ungroup() %>% 
  select(!c(PSI, Genotype)) %>% 
  distinct() %>% 
  group_by(Gene, variant, event) %>% 
  count() %>% 
  filter(n<2) %>%
  select(!n) %>% 
  left_join(comparison_data, by=c("Gene", "variant", "event"))-> full_table
  
  
full_table %>% 
  group_by(variant, event, Gene) %>% 
  nest() %>% 
  mutate(kw=map(data, ~kruskal.test(PSI ~ Genotype, data=.x))) %>% 
  #select(-data) %>% 
  mutate(kw_tidy=map(kw, broom::tidy)) %>% 
  unnest(kw_tidy) -> kw_outcomes
  
kw_outcomes %>% 
  select(-c(data, kw, statistic, parameter, method)) -> kw_save

kw_save %>% 
  arrange(p.value)

write_tsv(
  kw_save,
  paste0(gene, "_KW.tsv"),
  na = "NA",
  quote = "none",
  escape = c("double", "backslash", "none"),
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)

  
  
