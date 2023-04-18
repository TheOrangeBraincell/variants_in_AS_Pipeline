# 16.03.23
# Statistic Testing for ESR1, i.e. proof of concept.
# Mirjam Karlsson-Mller


library(GenomicRanges)
library(rstatix)
library(tidyverse)
#library(patchwork)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\Statistical_Test\\ESR1")

#read in prepped PSI file

PSI<- read_tsv("ESR1_PSI_prepped.tsv")

#This is prefiltered, so there should be no NaN rows. But lets check
#PSI %>%
#  mutate(s=rowSums(is.na(PSI))==ncol(PSI)-1) %>%
#  select(c(Location, s)) %>%
#  filter(s==FALSE)

#ok

#Read in prepped Genotype File
genotype<- read_tsv("prepped_ESR1_genotype_table.tsv")

#genotype
#-------------------------------
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
ranges<- read_tsv("../Database/gene_ranges.tsv")

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

#correlation_table %>%
#  filter(!(Genotype=="NE" | Genotype=="ND")) %>%
#  filter(!(is.nan(PSI))) -> filtered_correlation


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

#Correct for multiple testing: i dont know how many tests ill do in the end.
#But if we assume ESR1 as an average gene. And we had 14000 tests for it. And we have ca. 10000 genes. And correct for that?
#That should give us an idea.
kw_outcomes$q.values=p.adjust(kw_outcomes$p.value, method="fdr", n=140000000)

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

#---------------------------------
#For potential posthoc. But i think we should check interesting examples first. so wait with rest of script.
# 
# kw_outcomes %>%
#   unnest(cols=data) %>%
#   ungroup() %>%
#   select(-c(kw, statistic, parameter, method)) %>%
#   #select(Gene, event, variant, p.value, q.values) %>%
#   arrange(q.values) %>%
#   #group_by(event, q.values) %>%
#   #count()
#   #Filter for significant q values
#   filter(q.values<0.05) %>%
#   select(-c(p.value, q.values)) %>%
#   group_by(variant, event, Gene, Genotype) %>%
#   select(-Sample) %>%
#   summarize(PSI_v=list(PSI))-> prepare
# 
# prepare %>%
#   rename(Genotype="Genotype2") %>%
#   rename(PSI_v="PSI_v2") %>%
#   inner_join(prepare, by=c("Gene", "event", "variant")) %>%
#   #Remove redundancies
#   filter(Genotype !=Genotype2) %>%
#   mutate(Genotype=as.character(Genotype)) %>%
#   mutate(Genotype2=as.character(Genotype2)) %>%
#   mutate(event=as.character(event)) %>%
#   #make filter column
#   rowwise() %>%
#   #mutate(ID=paste(str_sort(c(Genotype, Genotype2)),collapse="_"))
#   mutate(ID=paste(Gene,event, variant,paste(str_sort(c(Genotype, Genotype2)),collapse="_"), sep="_")) %>%
#   distinct(ID, .keep_all = TRUE) %>%
#   select(-ID) %>%
#   #Once more remove pairs that have genotype NE or ND, because theres no comparison there.
#   filter(Genotype!="NE" & Genotype2 !="NE") %>%
#   #Now test. No need to state pairwise wilcox, because we manually put it into pairwise format.
#   mutate(wilcox_p=wilcox.test(PSI_v, PSI_v2)$p.value) %>%
#   #Because this test is done separately  by line, no correction for multiple testing has been done yet. Do manually.
#   mutate(wilcox_q=p.adjust(wilcox_p, method="fdr", n=nrow(.))) -> output_Wilcox
# 
# write_tsv(
#   output_Wilcox,
#   "output_wilcox_ESR1.tsv",
#   na = "NA",
#   quote = "none",
#   escape = c("double", "backslash", "none"),
#   eol = "\n",
#   num_threads = readr_threads(),
#   progress = show_progress()
# )


#--------------------------------
#Plots as well! Genotypes first

genotype<- read_tsv("ESR1_genotype_table.tsv", skip=8)

#variants in exons only!

rna_variants_filtered
genotype
  

genotype %>% 
  mutate(hmza=rowSums(genotype=="1/1"| genotype=="2/2" | genotype =="3/3")) %>% 
  mutate(hetz=rowSums(genotype=="0/1"| genotype=="1/2"| genotype=="0/2"| genotype=="2/3" | genotype=="0/3")) %>% 
  mutate(hmzr=rowSums(genotype=="0/0")) %>% 
  mutate(total=hmza+hmzr+hetz) %>% 
  #select(hmza, hetz, hmzr, total)
  mutate("Ref"=(hmzr/total)) %>% 
  mutate("Het"=hetz/total) %>% 
  mutate("Alt"=hmza/total) %>% 
  select(c(Location, Ref, Het, Alt)) %>% 
  arrange(desc(Ref)) %>% 
  mutate(n=seq.int(nrow(.))) %>% 
  pivot_longer(cols=!c(Location,n), names_to = "Genotype", values_to = "Fraction")-> genotype_fractions


ggplot(genotype_fractions, aes(x=n, y=Fraction, group=Genotype, color=Genotype))+
  geom_line(size=0.75)+
  theme_classic()+
  xlab("Variant Locations sorted by frequency of genotype homozygous reference")+
  ylab("Fractions of genotype occurance at location")+
  labs(color = "Genotype") +
  scale_color_brewer(labels = c("1/1", "0/1", "0/0"), palette ="Accent")



#So most of these variants are superrare. What if we want to plot those that arent?
genotype %>% 
  mutate(hmza=rowSums(genotype=="1/1"| genotype=="2/2" | genotype =="3/3")) %>% 
  mutate(hetz=rowSums(genotype=="0/1"| genotype=="1/2"| genotype=="0/2"| genotype=="2/3" | genotype=="0/3")) %>% 
  mutate(hmzr=rowSums(genotype=="0/0")) %>% 
  mutate(total=hmza+hmzr+hetz) %>% 
  filter(hetz+hmza>5) %>% 
  mutate("Ref"=(hmzr/total)) %>% 
  mutate("Het"=hetz/total) %>% 
  mutate("Alt"=hmza/total) %>% 
  select(c(Location, Ref, Het, Alt)) %>% 
  arrange(desc(Ref)) %>% 
  mutate(n=seq.int(nrow(.))) %>% 
  pivot_longer(cols=!c(Location,n), names_to = "Genotype", values_to = "Fraction")-> genotype_fractions


ggplot(genotype_fractions, aes(x=n, y=Fraction, group=Genotype, color=Genotype))+
  geom_line(size=0.75)+
  theme_classic()+
  xlab("Variant Locations sorted by frequency of 0/0")+
  ylab("Fractions of genotype occurance at location")+
  labs(color = "Genotype") +
  theme(text = element_text(size = 20), axis.text=element_text(size=20))+
  scale_color_brewer(labels = c("1/1", "0/1", "0/0"), palette ="Accent")

#Now is this what we expect?
setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\Statistical_Test\\ESR1")
swegen<- read_tsv("SweGen_ESR1_all_filtered_variants_phyloP_TargetScan.txt")

#We need the coordinates to lift them to 38
swegen %>% 
  select(c(Variant, Chrom, `Allele Frequency`)) %>% 
  mutate(Observed_Expected="expected") %>% 
  rename("Allele Frequency"="Allele_Frequency") %>% 
  mutate(chrom=paste0("chr", Chrom)) %>% 
  relocate(chrom, .before=Variant) %>% 
  select(!Chrom) %>% 
  separate(Variant, c("temp", "position", "ref", "alt"), sep="-") %>%
  mutate(position=as.integer(position)) %>% 
  select(!temp) %>% 
  #only point mutations no indels
  filter(nchar(alt)==1 & nchar(ref)==1) %>% 
  mutate(chromStop=position+1) %>% 
  mutate(chromStart=position) %>% 
  select(!position) %>% 
  relocate(chromStop, .after=chromStart) %>% 
  unite("name", c(ref, alt, Allele_Frequency, Observed_Expected), sep="_") %>% 
  relocate(name, .after=chromStop)->swegen_lift


write_tsv(swegen_lift, "swegen_lift_hg37.bed")

swegen38<- read_tsv("hglft_genome_27da8_d7c440.bed", col_names=FALSE)
  
swegen38 %>% 
  rename("chrom"="X1") %>% 
  rename("position"="X2") %>% 
  select(!c(X3, X5)) %>% 
  separate(X4, c("ref", "alt", "Allele_Frequency_Expected", "Observed_Expected"), sep="_") %>% 
  select(!Observed_Expected) %>% 
  mutate(Allele_Frequency_Expected=as.double(Allele_Frequency_Expected))->swegen_alt

#Now when we expand genotypes, check that its only locations in exons.

genotype %>%
  #Split infostring into separate columns.
  separate(col=Location, c("other", "alt"), sep="\\(") %>%
  separate(col=other, c("chrom", "position", "ref", "empty"), sep="_") %>%
  select(!empty) %>%
  mutate(alt=str_replace(alt, "\\)", "")) %>%
  mutate(position=as.integer(position)) %>% 
  mutate(hmza=rowSums(genotype=="1/1"| genotype=="2/2" | genotype =="3/3")) %>% 
  mutate(hetz=rowSums(genotype=="0/1"| genotype=="1/2"| genotype=="0/2"| genotype=="2/3" | genotype=="0/3")) %>% 
  mutate(hmzr=rowSums(genotype=="0/0")) %>% 
  mutate(total=hmza+hmzr+hetz) %>% 
  mutate("Ref"=(2*hmzr/(2*total))) %>% 
  mutate("Het"=hetz/(2*total)) %>% 
  mutate("Alt"=2*hmza/(2*total)) %>% 
  mutate('Allele_Frequency_observed'=(2*hmza + hetz)/(2*total)) %>% 
  select(c(chrom, position, ref, alt, Allele_Frequency_observed))-> genotypes_expanded


#Check which variants are in exons
setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\Database")
annotation<-read_tsv("gencode.v39.annotation.gff3", skip=7, col_names=FALSE)

annotation %>% 
  select(c(X1, X3, X4, X5)) %>% 
  filter(X1=="chr6") %>% 
  rename("chrom"="X1") %>% 
  rename("start"="X4") %>% 
  rename("stop"="X5") %>% 
  rename("Location"="X3") %>% 
  filter(Location=="exon")-> ann_interest

ann_interest      
rm(annotation)

gr_ann=makeGRangesFromDataFrame(ann_interest)
genotypes_expanded %>% 
  select(!c(ref, alt, Allele_Frequency_observed)) %>% 
  mutate(start=position) %>% 
  mutate(stop=position+1) %>% 
  makeGRangesFromDataFrame()-> gr_geno

coordinates<-as_tibble(findOverlaps(gr_geno, gr_ann))

genotypes_expanded

coordinates %>% 
  #mutate(ann_interest[subjectHits,1]) %>% 
  mutate(ann_interest[subjectHits,2]) %>% 
  mutate(genotypes_expanded[queryHits,1]) %>%
  mutate(genotypes_expanded[queryHits,2]) %>%
  mutate(genotypes_expanded[queryHits,3]) %>%
  mutate(genotypes_expanded[queryHits,4]) %>% 
  left_join(genotypes_expanded) %>% 
  select(!c(queryHits, subjectHits)) %>% 
  distinct() %>% 
  left_join(swegen_alt, by=c("alt", "ref", "chrom", "position")) %>% 
  #select(!c(chrom,Location)) %>% 
  drop_na() ->comparison_data
  

#pivot_longer(cols=!c(Location,n), names_to = "Genotype", values_to = "Fraction")
  #pivot_longer(cols=!c(position, ref, alt), names_to="flag", values_to="Allele_Frequency") %>% 
  #mutate("Observed_Expected"=ifelse(flag=="Allele_Frequency_observed", "observed", "expected")) %>% 
  #select(!flag) %>% 
  #mutate("Observed_Expected"=as.factor(Observed_Expected))-> comparison_data
setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\Statistical_Test\\ESR1")
write_tsv(comparison_data, "AF_Comparison_Table.tsv")

ggplot(comparison_data, aes(x=Allele_Frequency_observed, y=Allele_Frequency_Expected)) +
  theme_classic()+
  geom_point(size=2)+
  xlab("Allele Frequency observed in SCAN-B Analysis")+
  ylab("Allele Frequency observed in SweGen Analysis")+
  theme(text = element_text(size = 20), axis.text=element_text(size=20))
  #geom_segment(aes(x=position, y=Allele_Frequency,xend=position, yend=Allele_Frequency))


#What if the genotype frequency looks more reasonable too?

genotype %>% 
  mutate(hmza=rowSums(genotype=="1/1"| genotype=="2/2" | genotype =="3/3")) %>% 
  mutate(hetz=rowSums(genotype=="0/1"| genotype=="1/2"| genotype=="0/2"| genotype=="2/3" | genotype=="0/3")) %>% 
  mutate(hmzr=rowSums(genotype=="0/0")) %>% 
  mutate(total=hmza+hmzr+hetz) %>% 
  filter(hetz+hmza>5) %>% 
  mutate("Ref"=(hmzr/total)) %>% 
  mutate("Het"=hetz/total) %>% 
  mutate("Alt"=hmza/total) %>% 
  select(c(Location, Ref, Het, Alt))-> genotype_fractions


genotype_fractions %>% 
  separate(col=Location, c("other", "alt"), sep="\\(") %>%
  separate(col=other, c("chrom", "position", "ref", "empty"), sep="_") %>%
  select(!empty) %>%
  mutate(alt=str_replace(alt, "\\)", "")) %>%
  mutate(position=as.integer(position)) %>% 
  select(!c(ref, alt, Ref, Het, Alt)) %>% 
  mutate(start=position) %>% 
  mutate(stop=position+1) %>% 
  makeGRangesFromDataFrame()-> gr_geno

coordinates<-as_tibble(findOverlaps(gr_geno, gr_ann))


coordinates %>% 
  #mutate(ann_interest[subjectHits,1]) %>% 
  mutate(ann_interest[subjectHits,2]) %>% 
  mutate(genotype_fractions[queryHits,1]) %>%
  mutate(genotype_fractions[queryHits,2]) %>%
  mutate(genotype_fractions[queryHits,3]) %>%
  mutate(genotype_fractions[queryHits,4]) %>% 
  left_join(genotype_fractions) %>% 
  select(!c(queryHits, subjectHits)) %>% 
  distinct() %>% 
  arrange(desc(Ref)) %>%  
  mutate(n=seq.int(nrow(.))) %>% 
  pivot_longer(cols=!c(Location, n), names_to="Genotype", values_to="Fraction")-> genotype_freq

  

genotype_fractions

ggplot(genotype_freq, aes(x=n, y=Fraction, group=Genotype, color=Genotype))+
  geom_line(size=0.75)+
  theme_classic()+
  xlab("Variant Locations sorted by frequency of 0/0")+
  ylab("Fractions of genotype occurance at location")+
  labs(color = "Genotype") +
  theme(text = element_text(size = 20), axis.text=element_text(size=20))+
  scale_color_brewer(labels = c("1/1", "0/1", "0/0"), palette ="Accent")
  
  
#-----------------------------------
#PSI PLOTS

PSI<-read_tsv("ESR1_PSI.tsv")

PSI %>% 
  separate(col=Location, c("event", "A", "B", "C", "D"), sep="_") %>%  
  #There might be several headers. they will create NA values
  #drop_na() %>% 
  mutate(chrom=ifelse(startsWith(A, "chr"), A, B)) %>% 
  mutate(start=ifelse(startsWith(event, "A"), as.integer(D), as.integer(C))) %>% 
  mutate(stop=ifelse(startsWith(event, "A"), start+1, as.integer(D))) %>% 
  mutate(strand=ifelse(startsWith(event, "A"), C, B)) %>% 
  select(!c(B,C,D)) %>% 
  relocate(c(chrom, strand, start, stop), .after=event) ->PSI_expanded


PSI_expanded %>% 
  pivot_longer(cols=!c(event, chrom, strand, start, stop, A), names_to="Samples", values_to="PSI") %>% 
  drop_na() %>% 
  filter(event=="CE") %>% 
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
  mutate(bin=factor(bin,levels = c("[0,0.1)","[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)", "[0.4,0.5)","[0.5,0.6)","[0.6,0.7)","[0.7,0.8)","[0.8,0.9)","[0.9-1]")))-> histo_data

histo_data

ggplot(histo_data, aes(x=bin))+
  theme_classic()+
  geom_bar(fill="deepskyblue3")+
  xlab("PSI score interval")+
  ylab("Number of scores")+
  theme(text = element_text(size = 20), axis.text=element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(aes(label = ..count..), stat = "count", position="fill", vjust =1.5 )
  #ggtitle("Number of genes with an FPKM >10 in certain percentage of the samples")+
  #geom_text(aes(label = ..count..), stat = "count", position="fill", vjust =1.5 )+
  #scale_x_discrete(labels=c("0-1%","1-10%","10-20%","20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%"), guide=guide_axis(n.dodge=2))


ggplot(event_PSI, aes(x=PSI)) +
  theme_classic()+
  ylim(c(0,3))+
  geom_density(fill="darkolivegreen4", color="darkolivegreen4", alpha=0.8) +
  xlab("PSI score")+
  ylab("Density")


#---------------------------
#Plot for statistics results

library(GGally)

KW<-read_tsv("KW_ESR1_170323.tsv", skip=1)


ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()
