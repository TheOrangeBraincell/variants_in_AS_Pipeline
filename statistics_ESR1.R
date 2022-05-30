#!/usr/bin/env Rscript
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
# useage: Can be run in command line
#     Rscript --vanilla statistics_ESR1.R CE_ESR1_all.txt genotypes_in_exons.txt

# 0. Import library and set up inputs----
library(tidyverse)
library(optparse)

#Console input
option_list = list(
  make_option(c("-p", "--psi"), type="character", default=NULL, 
              help="PSI score table file name.", metavar="character"),
  make_option(c("-g", "--genotype"), type="character", default=NULL, 
              help="genotype table file name.", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="p_value.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Check that both input files are given
if (is.null(opt$psi)){
  print_help(opt_parser)
  stop("The psi score table must be provided as input with -p.", call.=FALSE)
}
if (is.null(opt$genotype)){
  print_help(opt_parser)
  stop("The genotype table must be provided as input with -g.", call.=FALSE)
}

# 1. Read in data -------
#Psi score table
psi_raw <- read_tsv(opt$psi, comment = "#")

psi<-psi_raw %>% 
  pivot_longer(cols = -Location, names_to = "Sample", values_to = "psi")

psi<-psi[,c("Sample", "Location", "psi")]
psi[order(psi$Sample),]
head(psi)
psi<-psi %>%
  rename(exon=Location)
#psi

#genotype table
genotype_raw<-read_tsv(opt$genotype, comment="#")
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

#write p-values and adjusted p-values in output file
sink(opt$out)
cat(p_frame)
sink()