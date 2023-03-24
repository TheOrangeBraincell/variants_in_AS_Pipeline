# 17.03.23
# Mirjam Karlsson-Müller
# Plotting ESR1 Stats outputs


library(tidyverse)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\Statistical_Test\\ESR1")


KW<-read_tsv("KW_ESR1_170323.tsv", skip=1)

head(KW)

#Because the correction for multiple testing might change, we extract the top 10 significant. 
#With correction for all genes, those would still be significant (as there was 17)
KW %>% 
  filter(p.value<0.00005) %>% 
  filter(!startsWith(event, "IR")) %>% 
  head(10)->sign_KW
  

sign_KW

#For those 10 we would like to extract PSI and Genotypes to see if having a variant makes the PSI scores smaller or bigger.


PSI<-read_tsv("ESR1_PSI_prepped.tsv")

sign_KW %>% 
  left_join(PSI, by=c("event"="Location")) %>% 
  select(!c(Gene, p.value)) %>% 
  pivot_longer(cols=!c(event, variant, q.values), names_to="Sample", values_to="PSI") %>% 
  drop_na()->PSI_KW

#rm(PSI)
PSI_KW

genotype<-read_tsv("prepped_ESR1_genotype_table.tsv")


genotype %>% 
  separate(col=Location, c("other", "alternative"), sep="\\(") %>% 
  separate(col=other, c("chrom", "position", "ref", "empty"), sep="_") %>% 
  mutate(variant=paste(chrom, sep="_", position)) %>% 
  select(!c(empty, chrom, position, ref, alternative)) %>% 
  pivot_longer(cols=!(variant), names_to="Sample", values_to="Genotype") %>% 
  filter(!(Genotype=="NE" | Genotype=="ND")) %>% 
  right_join(PSI_KW) %>% 
  drop_na()-> KW_data

# rm(genotype)
# rm(sign_KW)
# rm(KW)
# rm(PSI_KW)

#has all the information needed for plotting.
#Now i need to find out for each of the 17 tests, which genotypes affect the PSI scores how, to make a flag that decides the color of my plot dots later.

#is it possible to do a correlation plot with lil boxplots? Since we only have 17 significant tests.


  #11 variants, 8 events
  #so it would be 8 x 11 plot. that might be hard to see...
  #So maybe color is better anyways... But do I do a boxplot for each to see what direction it goes?


#3 events, 7 variants.

KW_data %>% 
  group_by(event, variant) %>% 
  count()

KW_data %>% 
  filter(event=="CE_chr6_+_152011655_152011794") %>% 
  filter(variant=="chr6_152099237") %>% 
  #drop_na() %>% 
  mutate(Genotype=as.factor(Genotype)) %>% 
  ggplot(aes(x=Genotype, y=PSI))+
  theme_classic()+
  geom_boxplot()

#All variants seem to lower the PSI scores...
#so there is no point colour coding them.
#instead make dot size or smth for p value?

sign_KW %>% 
  separate(variant, c("chrom", "position"), sep="_") %>% 
  separate(event, c("event", "chrom1", "strand", "start", "stop"), sep="_") %>% 
  select(c(event, position, start, stop, p.value)) %>% 
  mutate(event=paste(event, paste(start, stop, sep="-"), sep=":")) %>% 
  select(!c(start, stop)) %>% 
  mutate(bubble= case_when(p.value<1e-13 ~ "1e-14",
                           p.value>=1e-13 & p.value<=1e-12 ~ "1e-13",
                           p.value>1e-12 & p.value<=1e-11 ~ "1e-12",
                           p.value>1e-11 & p.value<=1e-10 ~ "1e-11",
                           p.value>1e-10 ~ "1e-10"))->sign_alt

library(viridis)

ggplot(sign_alt, aes(event, position, color=bubble))+
  geom_point(size=10)+
  theme_classic()+
  scale_color_viridis_d("P-value")+
  #scale_fill_viridis(discrete=TRUE)+
  xlab("Alternative Splicing Event")+
  ylab("Variant")+
  #scale_colour_discrete("P-value")+
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  theme(text = element_text(size = 15), axis.text=element_text(size=15))  


#the p-value scale is completely off....
sign_alt
