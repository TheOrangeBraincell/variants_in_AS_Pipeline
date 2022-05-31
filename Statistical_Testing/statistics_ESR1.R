#
# name: statistics_ESR1.R
# author: Mirjam Karlsson-Müller
# date: 24.05.22
# 
# description: Conducts a set of one way ANOVA tests, based on the output tables
# containing genotype and PSI score, to see if there is a 
# significant difference in p-values between genotypes within a potential 
# cassette exon.
# 

#Import library
library(tidyverse)


#Change working directory based on where you have saved your output files.
setwd("C:\\Users\\mirja\\Documents\\University\\ResearchProject\\Statistics")

# 1. Read in data -------
#Psi score table
psi_raw <- read_tsv("CE_ESR1_all.txt", comment = "#")

psi<-psi_raw %>% 
  pivot_longer(cols = -Location, names_to = "Sample", values_to = "psi")

psi<-psi[,c("Sample", "Location", "psi")]
psi[order(psi$Sample),]
head(psi)
psi<-psi %>%
  rename(exon=Location)
#psi

#genotype table
genotype_raw<-read_tsv("genotypes_in_exons.txt", comment="#")
#genotype_raw

genotype<-genotype_raw %>% 
  pivot_longer(cols = -c(Location,exon), names_to = "Sample" , values_to = "Genotype")

#genotype

genotype<-genotype[,c("Sample", "Location","exon", "Genotype")]
genotype[order(genotype$Sample),]

#Join Dataframes

data<-left_join(genotype, psi, by=c("Sample","exon"))
#remove rows with NaN PSI value.
data=drop_na(data)
data<-data[,c("Sample", "Location", "exon", "psi", "Genotype")]
data<-data[order(data$Sample),]

data$Genotype=factor(data$Genotype)


#Statistical testing -----
"One way anova for each pair of exons and location"

#initializing result data frame
p_frame<-data.frame(matrix(ncol=3, nrow=0))
colnames(p_frame)<-c("Exon", "Location", "pvalue")
# Extract exons, remove duplicates
exons<-c(data$exon)[!duplicated(data$exon)]
# loop through exons
for (e in exons){
  # make data subset for each exon
  data_exon<-data[data$exon==e,]
  # make variant list for each exon
  variants<-c(data_exon$Location)[!duplicated(data_exon$Location)]
  # loop through variants
  for (v in variants){
    # make data subset for each variant from previous subset for exon
    data_variant<-data_exon[data_exon$Location==v,]
    data_variant$Genotype<-factor(data_variant$Genotype)
    #If there is only one genotype left for this variant, then there is an
    # error. skip to next.
    tryCatch({
      # make linear model for subset 2, add result into new dataframe.
      fit.data<-lm(psi~Genotype, data=data_variant)
      summary_data<-summary(fit.data)
      p<-unname(pf(summary_data$fstatistic[1],summary_data$fstatistic[2], summary_data$fstatistic[3],
                   lower.tail=FALSE))
      df<-data.frame(e,v,p)
      colnames(df)<-c("Exon", "Location", "pvalue")
      p_frame<-rbind(p_frame, df)
    }, error=function(e){})
  }  
}

p_frame
#none of them are significant. Even less so after multiple testing correction

#Correct for multiple testing
p_frame$p_adjust<-p.adjust(c(p_frame$pvalue), method="bonferroni", n=length(c(p_frame$pvalue)))
#no significant result ._.
p_frame
#all 1's.

write.table(p_frame, file='p_values_ESR1.tsv', quote=FALSE, sep='\t')
