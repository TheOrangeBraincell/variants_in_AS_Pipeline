# 16.03.23
# Statistic Testing for ESR1, i.e. proof of concept.
# Mirjam Karlsson-Müller


library(tidyverse)
library(GenomicRanges)
library(rstatix)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\Statistical_Test\\ESR1")

#read in prepped PSI file

PSI<- read_tsv("ESR1_PSI_prepped.tsv")

#This is prefiltered, so there should be no NaN rows. But lets check
PSI %>% 
  mutate(s=rowSums(is.na(PSI))==ncol(PSI)-1) %>% 
  select(c(Location, s)) %>% 
  filter(s==FALSE)

#ok

#Read in prepped Genotype File
genotype<- read_tsv("prepped_ESR1_genotype_table.tsv")

genotype

#also no filtering necessary, as its only variants with more than one genotype at every location.

#Read in AS events
AS<-read_tsv("AS_events_ESR1.tsv")
#This one is needed to add gene information to every PSI scored event
AS %>% 
  drop_na() %>% 
  #The infostring has different format to compare to PSI table, it needs to have AS included, not in separate column
  mutate(merged=paste(AS_Type, sep="_", Location)) %>% 
  select(c(merged, Gene_Name)) %>% 
  #Merge with the PSI scores, to give gene information to each PSI score!
  right_join(PSI, by=c("merged"="Location")) %>% 
  pivot_longer(cols=!c(merged, Gene_Name), names_to="Sample", values_to="PSI")-> PSI_gene


#To match variant locations to the PSI scores in the same gene, we need to know what gene they are in
ranges<- read_tsv("gene_ranges.tsv")
ranges

genotype %>% 
  #Split infostring into separate columns.
  separate(col=Location, c("other", "alternative"), sep="\\(") %>% 
  separate(col=other, c("chrom", "position", "ref", "empty"), sep="_") %>% 
  select(!empty) %>% 
  mutate(alternative=str_replace(alternative, "\\)", "")) %>% 
  mutate(position=as.integer(position)) %>% 
  mutate(start=position) %>% 
  #gene ranges needs a start and a stop position, the stop is exclusive, so we add plus one to get a "range"
  mutate(stop=position+1)-> genotypes_expanded

gr<-makeGRangesFromDataFrame(ranges)
gl<-makeGRangesFromDataFrame(genotypes_expanded)


#find overlap between the two genomic ranges
coordinates<-as_tibble(findOverlaps(gl, gr))

#Add in the genotype information
coordinates %>% 
  mutate(ranges[subjectHits,1]) %>% 
  mutate(genotypes_expanded[queryHits,1]) %>% 
  mutate(genotypes_expanded[queryHits,2]) %>% 
  mutate(genotypes_expanded[queryHits,3]) %>% 
  mutate(genotypes_expanded[queryHits,4]) %>% 
  left_join(genotypes_expanded) %>% 
  select(!c(queryHits, subjectHits, start, stop)) %>% 
  group_by(Gene) %>% 
  pivot_longer(cols=!c(Gene, position, ref, alternative, chrom), names_to = "Sample", values_to="Genotype") -> gt_gene

rm(AS)
rm(PSI)
rm(coordinates)
rm(gl)
rm(gr)
rm(ranges)
rm(genotype)
rm(genotypes_expanded)

#Its too big here to process the next step. But its looking good! 

#Join dataframes based on genes!
gt_gene %>% 
  left_join(PSI_gene, by=c("Gene"="Gene_Name", "Sample")) %>% 
  #Remove variants that couldnt be given to an event
  filter(merged!="NA") %>% 
  filter(!(Genotype=="NE" | Genotype=="ND")) %>% 
  filter(!(is.nan(PSI))) -> filtered_correlation



#We now have entries with non numeric PSI value, which we cant use for statistical testing.
#We also cant do stats for genotype NE and ND. So all of those rows need to be removed

correlation_table %>% 
  filter(!(Genotype=="NE" | Genotype=="ND")) %>% 
  filter(!(is.nan(PSI))) -> filtered_correlation

#Now We removed a lot of rows. So there might be locations with only one genotype left that has numerical PSI values.
#So we need to check that there is still 2 groups of PSI values to compare!

filtered_correlation %>%
  select(c(Gene, chrom, position, ref, alternative, Genotype, merged)) %>% 
  distinct() %>% 
  group_by(Gene, chrom, position, ref, alternative, merged) %>% 
  count() %>% 
  filter(n>=2) %>% 
  select(!n) %>% 
  left_join(filtered_correlation)-> testable_pairs


testable_pairs %>% 
  ungroup() %>% 
  mutate(Genotype=as.factor(Genotype)) %>% 
  mutate(event=as.factor(merged)) %>% 
  select(!merged) %>% 
  mutate(variant=paste(chrom, sep="_", position)) %>% 
  select(!c(chrom, position, ref, alternative)) %>% 
  group_by(variant, event, Gene) %>% 
  drop_na() %>% 
  nest() %>% 
  mutate(kw=map(data, ~kruskal.test(PSI ~ Genotype, data=.x))) %>% 
  #select(-data) %>% 
  mutate(kw_tidy=map(kw, broom::tidy)) %>% 
  unnest(kw_tidy) -> kw_outcomes

#Correct for multiple testing
kw_outcomes$q.values=p.adjust(kw_outcomes$p.value, method="fdr", n=length(kw_outcomes$p.value))

#format kw_outcomes for potential save
kw_outcomes %>% 
  select(-c(data, kw, statistic, parameter, method)) -> kw_save

kw_save %>% 
  arrange(q.values)-> kw_save

write_tsv(
  kw_save,
  "output_KW_ESR1.tsv",
  na = "NA",
  quote = "none",
  escape = c("double", "backslash", "none"),
  eol = "\n",
  num_threads = readr_threads(),
  progress = show_progress()
)

kw_outcomes %>% 
  unnest(cols=data) %>%
  ungroup() %>% 
  select(-c(kw, statistic, parameter, method)) %>% 
  #select(Gene, event, variant, p.value, q.values) %>% 
  arrange(q.values) %>% 
  group_by(event, q.values) %>%
  count()
  #Filter for significant q values
  filter(q.values<0.05) %>% 
  select(-c(p.value, q.values)) %>% 
  group_by(variant, event, Gene, Genotype) %>% 
  select(-Sample) %>% 
  summarize(PSI_v=list(PSI))-> prepare

  prepare %>%
    rename(Genotype="Genotype2") %>% 
    rename(PSI_v="PSI_v2") %>% 
    inner_join(prepare, by=c("Gene", "event", "variant")) %>% 
    #Remove redundancies
    filter(Genotype !=Genotype2) %>% 
    mutate(Genotype=as.character(Genotype)) %>% 
    mutate(Genotype2=as.character(Genotype2)) %>% 
    mutate(event=as.character(event)) %>% 
    #make filter column
    rowwise() %>% 
    #mutate(ID=paste(str_sort(c(Genotype, Genotype2)),collapse="_"))
    mutate(ID=paste(Gene,event, variant,paste(str_sort(c(Genotype, Genotype2)),collapse="_"), sep="_")) %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    select(-ID) %>% 
    #Once more remove pairs that have genotype NE or ND, because theres no comparison there.
    filter(Genotype!="NE" & Genotype2 !="NE") %>% 
    #Now test. No need to state pairwise wilcox, because we manually put it into pairwise format.
    mutate(wilcox_p=wilcox.test(PSI_v, PSI_v2)$p.value) %>% 
    #Because this test is done separately  by line, no correction for multiple testing has been done yet. Do manually.
    mutate(wilcox_q=p.adjust(wilcox_p, method="fdr", n=nrow(.))) -> output_Wilcox
  
  write_tsv(
    output_Wilcox,
    "output_wilcox_ESR1.tsv",
    na = "NA",
    quote = "none",
    escape = c("double", "backslash", "none"),
    eol = "\n",
    num_threads = readr_threads(),
    progress = show_progress()
  )




