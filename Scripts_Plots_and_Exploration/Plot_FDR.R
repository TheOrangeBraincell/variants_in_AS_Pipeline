# Plot_FDR.R
#Mirjam Karlsson-Müller
# 04.04.23

library(tidyverse)
library(GenomicRanges)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\FDR_Outputs")


fdr<-read_tsv("merged_FDR_table.tsv")

fdr %>% 
  arrange(q.values)

fdr %>% 
  separate(col=event, c("A","B", "C", "D", "E"), sep="_") %>% 
  mutate(AS_Type=ifelse(startsWith(B, "A"), B, A)) %>% 
  #mutate(AS_Type= substr(event, 1, 2)) %>% 
  group_by(AS_Type) %>% 
  count()

# A tibble: 4 x 2
# # Groups:   AS_Type [4]
# AS_Type      n
# <chr>    <int>
# 1 AA       82369
# 2 AD        7070
# 3 CE      165177
# 4 IR       11335

fdr %>% 
  separate(col=event, c("A","B", "C", "D", "E"), sep="_") %>% 
  mutate(AS_Type=ifelse(startsWith(B, "A"), B, A)) %>% 
  filter(AS_Type!="IR")

#254'616 significant without IR events.

#How many interactions in same exon? How many of those have the variant in the splice site? (4 bases before intron)

#For exons need a gff file.

gencode<- read_tsv("..\\Database\\gencode.v39.annotation.gff3",skip=7, col_names=FALSE) %>% 
  select(c(X1, X3, X4, X5)) %>% 
  filter(X3=="exon") %>% 
  mutate("chrom"=X1) %>% 
  mutate("start"=X4) %>% 
  mutate("stop"=X5) %>% 
  mutate("Location"=X3) %>% 
  select(!c(X1, X3, X4, X5))

#Match variants to gencode coordinates to see if they are in exons.

fdr %>% 
  filter(!startsWith(event, "IR")) %>% 
  select(Gene, variant) %>% 
  separate(col=variant, c("other", "alt"), sep="\\(") %>%
  separate(col=other, c("chrom", "position", "ref", "empty"), sep="_") %>% 
  select(!empty) %>% 
  mutate(alt=str_replace(alt, "\\)", "")) %>%
  mutate(start=as.integer(position)) %>% 
  mutate(stop=start+1) -> fdr_genotypes

fdr_genotypes %>%   
  makeGRangesFromDataFrame()-> gr_geno

gr_ann=makeGRangesFromDataFrame(gencode)  

coordinates<-as_tibble(findOverlaps(gr_geno, gr_ann))

coordinates %>% 
  mutate(fdr_genotypes[queryHits,1]) %>% 
  mutate(fdr_genotypes[queryHits,2]) %>% 
  mutate(fdr_genotypes[queryHits,3]) %>% 
  mutate(fdr_genotypes[queryHits,4]) %>% 
  mutate(fdr_genotypes[queryHits,5]) %>% 
  mutate(gencode[subjectHits, 2]) %>% 
  mutate(gencode[subjectHits, 3]) %>% 
  mutate(exon_start_v=start) %>% 
  mutate(exon_stop_v=stop) %>% 
  select(!c(start, stop, queryHits, subjectHits)) -> variants_in_exons

variants_in_exons %>% 
  mutate(Alt_Base=paste0("(", alt, ")")) %>% 
  mutate(variant=paste(chrom, position, ref,Alt_Base, sep="_")) %>% 
  left_join(fdr) %>% 
  distinct() %>% 
  select(!c(chrom, position, ref, alt, Alt_Base)) %>% 
  group_by(Gene, variant, event, p.value, q.values) %>% 
  count()

#221'861 variants are in exons. 

#Now how many are in the same exon as their correlated splice event?

variants_in_exons %>% 
  mutate(Alt_Base=paste0("(", alt, ")")) %>% 
  mutate(variant=paste(chrom, position, ref,Alt_Base, sep="_")) %>% 
  left_join(fdr) %>% 
  distinct() %>% 
  select(!c(chrom, position, ref, alt, Alt_Base)) -> fdr_variants_in_exons

#So i can merge them later, without size problem
fdr_variants_in_exons %>% 
  select(!c(p.value, q.values)) -> fdr_merge

fdr_variants_in_exons %>% 
  select(c(Gene, event)) %>% 
  separate(event, c("A", "B", "C", "D", "E"), sep="_") %>% 
  mutate(AS_Type=ifelse(startsWith(B, "A"), B, A)) %>% 
  mutate(chrom=ifelse(startsWith(AS_Type, "A"), C, B)) %>% 
  mutate(strand=ifelse(startsWith(AS_Type, "A"), D, C)) %>% 
  mutate(start=ifelse(startsWith(AS_Type, "A"), as.integer(E), as.integer(D))) %>% 
  mutate(stop=ifelse(startsWith(AS_Type, "A"), start+1, as.integer(E))) %>% 
  unite("event",c(A, B, C, D ,E)) -> fdr_PSI

#fdr_PSI

fdr_PSI %>% 
  makeGRangesFromDataFrame()-> gr_psi

coordinates<-as_tibble(findOverlaps(gr_psi, gr_ann))

coordinates %>% 
  mutate(fdr_PSI[queryHits,1]) %>% 
  mutate(fdr_PSI[queryHits,2]) %>% 
  #mutate(fdr_PSI[queryHits,3]) %>% 
  #mutate(fdr_PSI[queryHits,5]) %>% 
  #mutate(fdr_PSI[queryHits,6]) %>% 
  mutate(gencode[subjectHits, 2]) %>% 
  mutate(gencode[subjectHits, 3]) %>%
  mutate(exon_start_e=start) %>% 
  mutate(exon_stop_e=stop) %>% 
  select(!c(start, stop, queryHits, subjectHits))-> as_in_exons

rm(coordinates)
#rm(fdr)
rm(gencode)
rm(gr_ann)
rm(gr_psi)
rm(gr_geno)
rm(variants_in_exons)
rm(fdr_genotypes)
rm(fdr_PSI)

as_in_exons %>% 
  distinct() %>% 
  left_join(fdr_merge, by=c("exon_start_e"="exon_start_v", "exon_stop_e"="exon_stop_v", "Gene", "event")) %>% 
  drop_na() %>% 
  relocate(variant, .after=event) %>% 
  distinct()-> same_exon
 
write_tsv(same_exon,
    paste0("Significant_q_same_exon.tsv"),
    na = "NA",
    quote = "none",
    escape = c("double", "backslash", "none"),
    eol = "\n",
    num_threads = readr_threads(),
    progress = show_progress())

#Now we want variants in splice region. Which means for CE in beginning or end of exon.
#For AA and AD its only interesting if its on the side of the event...
#So its strand specific -.-

#same_exon %>% 
#  filter(Gene=="ESR1")

#1 pair in same exon. 
same_exon %>% 
  select(!c(exon_start_e, exon_stop_e)) %>% 
  distinct()

same_exon %>% 
  separate(event, c("A", "B", "C", "D", "E"), sep="_") %>% 
  mutate(AS_Type=ifelse(startsWith(B, "A"), B, A)) %>% 
  mutate(chrom=ifelse(startsWith(AS_Type, "A"), C, B)) %>% 
  mutate(strand=ifelse(startsWith(AS_Type, "A"), D, C)) %>% 
  mutate(start=ifelse(startsWith(AS_Type, "A"), as.integer(E), as.integer(D))) %>% 
  mutate(stop=ifelse(startsWith(AS_Type, "A"), start+1, as.integer(E))) %>% 
  unite("event",c(A, B, C, D ,E)) %>% 
  separate(col=variant, c("other", "alt"), sep="\\(") %>%
  separate(col=other, c("chrom", "position", "ref", "empty"), sep="_") %>% 
  select(!empty) %>% 
  mutate(alt=str_replace(alt, "\\)", "")) %>%
  mutate(position=as.integer(position)) %>% 
  mutate(Alt_Base=paste0("(", alt, ")")) %>% 
  mutate(variant=paste(chrom, as.character(position), ref,Alt_Base, sep="_")) %>% 
  select(!Alt_Base) %>% 
  relocate(variant, .after=event) %>% 
  mutate(splice_end= exon_stop_e-position) %>% 
  mutate(splice_start=position-exon_start_e) %>% 
  #arrange(splice_end) %>% 
  filter((splice_end<4) | (splice_start<4)) %>% 
  #Now for AA and AD they need to be on the right side.
  filter((AS_Type=="AA" & strand=="-" & splice_end<4) | (AS_Type=="AD" & strand=="+" & splice_end<4) | AS_Type=="CE" | (AS_Type=="AD" & strand=="-" & splice_start<4) | (AS_Type=="AA"& strand=="+"& splice_start<4) ) -> splice_sites_expanded
  
  
#Write into output table
splice_sites_expanded

splice_sites_expanded %>% 
  #relocate(chrom, .after=Gene) %>% 
  select(!c(chrom, position, ref, alt, AS_Type, strand, start, stop, splice_end, splice_start)) %>% 
  left_join(fdr) %>%
  distinct() %>% 
  arrange(q.values)-> splice_sites 

write_tsv(splice_sites, "Splice_Site_pairs.tsv")
  
splice_sites %>% 
  #The event/variant coordinates got linked to several different exons.
  #That makes for duplicate pairs in the list, eventhough it is just one test. So we remove those.
  select(!c(exon_start_e, exon_stop_e)) %>% 
  distinct() %>% 
  arrange(q.values)%>% 
  write_tsv("Splice_Site_pairs_noexon.tsv")

#841 significantly correlated events with variant in splice region.

#Looking at examples

splice_sites %>% 
  filter(event=="CE_chr17_-_7673700_7673837") %>% 
  filter(variant=="chr17_7673700_C_(A_T_G)")


PSI<-read_tsv("Interesting_PSI.tsv") %>% 
  pivot_longer(cols=!Location, names_to="Sample", values_to="PSI") %>% 
  filter(!is.nan(PSI))

genotypes<-read_tsv("Interesting_Genotypes.tsv")

splice_sites %>% 
  right_join(genotypes, by=c("variant"="Location", "Gene")) %>% 
  select(!c(exon_stop_e, exon_start_e, p.value, q.values)) %>% 
  pivot_longer(cols=!c(Gene, event, variant), names_to="Sample", values_to="Genotype") %>% 
  filter(Genotype!="NE") %>% 
  filter(Genotype!="ND") %>% 
  right_join(PSI, by=c("event"="Location", "Sample")) -> concat_table


#lets see which one of them would be nice. One with several genotypes perhaps?

concat_table %>% 
  group_by(Gene, event, variant, Genotype) %>% 
  count() %>% 
  filter(Gene=="PAK1")

# A tibble: 3 x 5
# Groups:   Gene, event, variant, Genotype [3]
# Gene  event                        variant              Genotype     n
# <chr> <chr>                        <chr>                <chr>    <int>
#   1 PAK1  CE_chr11_-_77379893_77379994 chr11_77379894_C_(T) 0/0       3121
# 2 PAK1  CE_chr11_-_77379893_77379994 chr11_77379894_C_(T) 0/1        108
# 3 PAK1  CE_chr11_-_77379893_77379994 chr11_77379894_C_(T) 1/1          9

#Its the only one who has more than 1 alternative allele.... -.-

concat_table %>% 
  filter(Gene=="TP53") %>% 
  mutate(Genotype=as.factor(Genotype)) %>% 
  ggplot(aes(x=Genotype, y=PSI))+
  theme_classic()+
  geom_boxplot()+
  geom_point()+
  ylim(0.8, 1)

#This ones interesting cause its associated with like half the cancers out there.
#So having a mutation spot, especially one where every nucleotide is found, is cool

#So this one is actually a known variant, causing intron retention in the space after. 
#which ofc we dont have reliable psi scores for -.-
#But it could be proof at least that we find it?


#Fishing for other examples (note that the files are overwritten...)
PSI<-read_tsv("Interesting_PSI.tsv") %>% 
  distinct() %>% 
  pivot_longer(cols=!Location, names_to="Sample", values_to="PSI") %>% 
  filter(!is.nan(PSI))

genotypes<-read_tsv("Interesting_Genotypes.tsv") %>% 
  distinct() %>% 
  pivot_longer(cols=!c(Location, Gene), names_to="Sample", values_to="Genotype") %>% 
  filter(Genotype!="NE") %>% 
  filter(Genotype!="ND")

splice_sites %>% 
  select(!c(exon_start_e, exon_stop_e)) %>% 
  distinct() %>% 
  left_join(genotypes, by=c("Gene", "variant"="Location")) %>% 
  left_join(PSI, by=c("event"="Location", "Sample")) -> concat_table
  
concat_table %>% 
  select(!c(p.value,q.values)) %>% 
  group_by(Gene, event, variant, Genotype) %>% 
  count() %>%
  filter(n>10) %>% 
  filter(n<3444) %>%
  #arrange(desc(n)) %>% 
  group_by(Gene, event, variant) %>% 
  count() %>% 
  filter(n>1) -> candidates

candidates %>% 
  ungroup() %>% 
  select(Gene) %>% 
  left_join(concat_table) %>% 
  select(!c())
  group_by(Gene, event, variant, Genotype) %>% 
  count()-> overview
  
concat_table %>% 
  filter(Gene=="APIP") %>% 
  mutate(Genotype=as.factor(Genotype)) %>% 
  ggplot(aes(x=Genotype, y=PSI))+
  theme_classic()+
  geom_jitter(aes(col=Genotype))+
  geom_boxplot(alpha = 0.8)+
  theme(text = element_text(size = 20), axis.text=element_text(size=20))

concat_table %>% 
  filter(Gene=="APIP") %>% 
  mutate(Genotype=as.factor(Genotype)) %>% 
  group_by(Gene, event, variant, Genotype) %>% 
  count()

#rs61734605 from C to T

# # A tibble: 3 x 5
# # Groups:   Gene, event, variant, Genotype [3]
# Gene  event                        variant                Genotype     n
# <chr> <chr>                        <chr>                  <fct>    <int>
#   1 APIP  CE_chr11_-_34895009_34895110 chr11_34895110_C_(T_A) 0/0        583
# 2 APIP  CE_chr11_-_34895009_34895110 chr11_34895110_C_(T_A) 0/1       1408
# 3 APIP  CE_chr11_-_34895009_34895110 chr11_34895110_C_(T_A) 1/1        395

splice_sites %>% 
  filter(Gene=="APIP")
