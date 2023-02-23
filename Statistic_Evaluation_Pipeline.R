# 23-02-23
# Mirjam Karlsson-Müller
# Statistic_Evaluation_Pipeline.R
# 
# Description:
#   Reads in PSI tables and Genotype tables. Then it matches the PSI scores of events to the variants in the same gene and looks for significant correlations
#   between the PSI value for different genotypes (Kruskal Wallis Test). This test was chosen, as the PSI scores most commonly are very close to 1 or very
#   close to 0, thus do not follow a normal distribution.
#   For the significant correlations between genotypes and PSI values, a paired Wilcox test was conducted and then corrected to see which genotype has the 
#   significant correlation.
#   
  

library(tidyverse)
library(GenomicRanges)
library(rstatix)
library(nplyr)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\Statistical_Test\\Statistical_Test_Practice_Data")

#Need inputs matching the final results from the pipeline, but much smaller. 
genotypes<-read_tsv("genotype_table_ESR1_new.tsv", skip =8)

PSI<-read_tsv("PSI_ESR1_new.tsv")

# Plan of approach:
#
# Three types of tests. 
# Filter outputs first, then run tests for each event against all variants in same gene.
# Correct for multiple testing.



#1 Filter: Remove lines in PSI table, where all entries are NA

PSI %>% 
  mutate(s=rowSums(is.na(PSI))==ncol(PSI)-1) %>% 
  filter(s==FALSE)%>% 
  select(- s) %>% 
  #Also remove rows where all PSI scores are 1 or 0. We do that by checking if all are the same. Without psi variation, theres no interest.
  mutate(s= if_all(starts_with('S0'), ~ S000001==.x)) %>% 
  filter(s==FALSE) %>% 
  select(-s)->PSI_filtered

#2 Filter: Remove lines in genotype table, where no entry has genotype homozygous reference 0/0
# Note that in the practice data, there is still rows with only NE and ND. These are removed in the real data.

genotypes %>% 
  mutate(hmzr=rowSums(genotypes=="0/0")) %>% 
  #For practice data. Remove ND/NE rows
  mutate(r=rowSums(genotypes=="NE")) %>% 
  mutate(t=rowSums(genotypes=="ND")) %>% 
  mutate(unclassified=r+t) %>% 
  #We also want to remove lines, that have only one genotype that passes filtering. (Again, wont be needed for final data.)
  mutate(hmza=rowSums(genotypes=="1/1")) %>% 
  mutate(hetz=rowSums(genotypes=="0/1")) %>% 
  #For there to be 2 genotypes, out of x, y, s at least two have to be not 0.
  mutate(filter_cond=ifelse((hmzr==0 | unclassified==ncol(genotypes)-1 |(hmzr==0 & hmza==0) | (hmzr==0 & hetz==0)| (hmza==0 & hetz==0)), T, F)) %>% 
  filter(filter_cond==FALSE) %>% 
  select(- c(hmzr,r,t,unclassified,hetz,hmza, filter_cond))-> genotypes_filtered
  

# remove unfiltered dataframes to prevent using too much space.
rm(PSI)
rm(genotypes)


# 3 Pair events with variants. 
# Every event gets paired with variants in the same gene.
# Means I need gene information ._.
# AS_event table contains gene information for each exon.

AS<- read_tsv("AS_events_ESR1.tsv")

#The gene headers are a problem. Remove them.To join them with PSI, requires AS type to be in location
AS %>% 
  drop_na() %>% 
  mutate(merged=paste(AS_Type, sep="_", Location)) %>% 
  select(c(merged, Gene_Name)) %>% 
  #Merge with the PSI scores, to give gene information to each PSI score!
  right_join(PSI_filtered, by=c("merged"="Location")) %>% 
  pivot_longer(cols=!c(merged, Gene_Name), names_to="Sample", values_to="PSI")-> PSI_gene

rm(AS)
# Now I need to link the variants to the same gene names. But the AS_events does not contain gene ranges.
# However, gene_ranges does. So we import that.

ranges<- read_tsv("gene_ranges.tsv")

genotypes_filtered %>% 
  separate(col=Location, c("other", "alternative"), sep="\\(") %>% 
  separate(col=other, c("chrom", "position", "ref", "empty"), sep="_") %>% 
  select(!empty) %>% 
  mutate(alternative=str_replace(alternative, "\\)", "")) %>% 
  mutate(position=as.integer(position)) %>% 
  mutate(start=position) %>% 
  mutate(stop=position+1)-> genotypes_expanded

rm(genotypes_filtered)
#Create genomic ranges objects for genotypes and genetic ranges

gr<-makeGRangesFromDataFrame(ranges)
gl<-makeGRangesFromDataFrame(genotypes_expanded)


#find overlaps

coordinates<-as_tibble(findOverlaps(gl, gr))

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

rm(coordinates, gr, gl)
#Now we should be able to join dataframes based on genes!
gt_gene %>% 
  left_join(PSI_gene, by=c("Gene"="Gene_Name", "Sample")) %>% 
  #Remove variants that couldnt be given to an event
  filter(merged!="NA") -> correlation_table

rm(gt_gene)
#HOORAAAAYYY the table is in the right format.

#1. Logistic Regression: Transform psi values
# He suggested id use a logit transformation. But that doesnt work for PSI=1 or PSI=0. Which is the majority of the values... So nvm that. lets try linear first.

#Lets try on one event first.
# correlation_table %>% 
#   filter(merged=="AD_2_chr6_+_151807509") %>% 
#   filter(Genotype==c("1/1","0/0", "0/1")) %>% 
#   drop_na() %>% 
#   mutate(position=as.factor(position))-> test
# test

# x<-lm(PSI ~ Genotype, data=test)
# summary(x)

#Also these arent normally distributed probably. Lets try ANOVA.
# anova<-aov(PSI ~ Genotype, data=test)
# summary(anova)

# ... its worse isnt it. I need to test for each event x variant.
# which means double loop.
# This is how to get the p-value
# summary(anova)[[1]][["Pr(>F)"]][1]
# empty dataframe tibble
# p_values<-as_tibble(matrix(nrow = 0, ncol = 1), .name_repair=~ c("p-values"))



#This is an attempt to do it tidy. But its not working ^^
# Nested is a good attempt, but basically we want each "nest" to have its own test. 
#which means id need nests, within nests, within nests. Because I need to test gene x position x event
correlation_table %>% 
  mutate(Genotype=as.factor(Genotype)) %>% 
  mutate(event=as.factor(merged)) %>% 
  select(!merged) %>% 
  mutate(variant=paste(chrom, sep="_", position)) %>% 
  select(!c(chrom, position, ref, alternative)) %>% 
  group_by(variant, event, Gene) %>% 
  drop_na() -> full_table
  #So after dropping the NA, there might be only numeric PSI values for one genotype... So we need to check again.
  #Split it into two dataframes... Filter, rejoin.


full_table %>% 
  select(!c(Genotype)) -> full_PSI

rm(correlation_table)

# ANOVA IMPLEMENTATION, BUT WE PROBABLY DO NOT FULFILL REQUIREMENT OF NORMALITY. HENCE COMMENTED OUT.

# full_table %>% 
#   select(!PSI) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from=Sample, values_from="Genotype") %>%  
#   #Check if each row has more than 1 type of genotype
#   mutate(s= ifelse(if_all(starts_with('S0'), ~ S000001==.x), T, F)) %>% 
#   mutate(s= ifelse(is.na(s), T, F)) %>% 
#   filter(s==FALSE) %>% 
#   select(!s) %>% 
#   pivot_longer(cols=!c(Gene, event, variant), names_to="Sample", values_to="Genotype") %>% 
#   drop_na() %>% 
#   left_join(full_PSI) %>% 
#   group_by(variant, event, Gene) %>% 
#   nest() %>% 
#   mutate(aov=map(data, ~aov(PSI ~ Genotype, data=.x))) %>% 
#   select(-data) %>% 
#   mutate(aov_tidy=map(aov, broom::tidy)) %>% 
#   unnest(aov_tidy)-> aov_outcomes
# 
# 
# 
# #Multiple Testing Correction with FDR. yay.
# aov_outcomes$q.values=p.adjust(aov_outcomes$p.value, method="fdr", n=length(aov_outcomes$p.value))
# 
# #In case we want this saved
# aov_outcomes %>% 
#   select(-c(aov, term, df, sumsq, meansq, statistic))

#However, we can really only use this, if the residuals are normally distributed. Which i have no idea if they are. But my mind says its unlikely. When we have the whole dataset, we have to check!
# Actually the nature of PSI scores suggests that they are always gonna be closest to 1 or 0 and barely any in the middle, which is the opposite to normal distribution.
#Furthermore, being this large of a data set means a shapiro test (f.e.) would be oversensitive and be significant regardless of how the data looks. 
#Do i have enough data points to say that by the CLT normality is given?
#Does this even work when its the residuals that are supposed to be normal?

#aov_outcomes %>% 
#  mutate(residuals=aov$residuals)
#  mutate(norm_res= shapiro.test(aov$residuals))

#In case they are not, we prepare by also doing a "non-parametric ANOVA", i.e. Kruskal Wallis test.
#This is from the R exercises last year, as reference: file "R ref script ex 7 ht21.R"
  # Carry out non-parametric Kruskal-Wallis test 
  # kruskal.test(len~treat, data=peas)

  # Non-parametric posthoc test: pairwise Wilcoxon signed-ranks test with FDR correction
  # pairwise.wilcox.test(peas$len, peas$treat, p.adjust.method="fdr")

full_table %>% 
  select(!PSI) %>% 
  ungroup() %>% 
  pivot_wider(names_from=Sample, values_from="Genotype") %>%  
  #Check if each row has more than 1 type of genotype
  mutate(s= ifelse(if_all(starts_with('S0'), ~ S000001==.x), T, F)) %>% 
  mutate(s= ifelse(is.na(s), T, F)) %>% 
  filter(s==FALSE) %>% 
  select(!s) %>% 
  pivot_longer(cols=!c(Gene, event, variant), names_to="Sample", values_to="Genotype") %>% 
  drop_na() %>% 
  left_join(full_PSI) %>% 
  group_by(variant, event, Gene) %>% 
  nest() %>% 
  mutate(kw=map(data, ~kruskal.test(PSI ~ Genotype, data=.x))) %>% 
  #select(-data) %>% 
  mutate(kw_tidy=map(kw, broom::tidy)) %>% 
  unnest(kw_tidy) -> kw_outcomes

#remove previous dataframe
rm(full_table)

#Correct for multiple testing
kw_outcomes$q.values=p.adjust(kw_outcomes$p.value, method="fdr", n=length(kw_outcomes$p.value))

#format kw_outcomes for potential save
kw_outcomes %>% 
  select(-c(data, kw, statistic, parameter, method)) -> kw_save

#Sort by q values
kw_outcomes %>% 
  unnest(cols=data) %>%
  ungroup() %>% 
  select(-c(kw, statistic, parameter, method)) %>% 
  #select(Gene, event, variant, p.value, q.values) %>% 
  arrange(q.values) %>% 
  #Filter for significant q values
  filter(q.values<0.05) %>% 
  select(-c(p.value, q.values)) %>% 
  group_by(variant, event, Gene, Genotype) %>% 
  select(-Sample) %>% 
  summarize(PSI_v=list(PSI))-> prepare

#This last one throws an error if there is no significant results for Kruskal Wallis!

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
  mutate(wilcox_q=p.adjust(wilcox_p, method="fdr", n=nrow(.)))



