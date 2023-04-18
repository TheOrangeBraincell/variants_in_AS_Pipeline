# Variants in ESR1

library(tidyverse)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs")

variants<- read_tsv("genotype_ESR1_Variants_AS.txt")


#Check how many non unique variants:

filtered_variants <- variants %>% 
  pivot_longer(cols = !Location, names_to = "Patient", values_to = "genotype") %>% 
  filter(!genotype %in% c("ND", "NE")) %>% 
  mutate(interesting = ifelse(genotype == "0/0", F, T)) %>% 
  filter(interesting) %>% 
  count(Location, interesting, sort=T) %>% 
  filter(n>=5)

#25 locations with more than 5 samples showing a variant. 

#Import ESR1 coordinates

ESR1<-read_tsv("rewritten_ESR1_coor.txt")
spec(ESR1)
ESR1

ESR1 %>% 
  select(exonStarts, exonEnds) %>% 
  distinct()


x<-filtered_variants %>% 
  mutate(coordinates=as.numeric(str_remove(Location, "chr\\d_")))

in_exon_vector<-vector()
for (i in x$coordinates){
  in_exon="in intron"
  for (row in 1:nrow(ESR1)){
    if (i>= ESR1[row, "exonStarts"] & i<=ESR1[row, "exonEnds"]){
      in_exon="in exon"
    } 
  }
  in_exon_vector=cbind(in_exon_vector, in_exon)
}


x<-cbind(x, matrix(in_exon_vector))


#Adjusting code from internet for chromosome plot?

vector_x<-pull(x, coordinates)


ggplot() + geom_vline(aes(xintercept=x)) + 
  xlim(c(151690496, 152103274))+
  theme_bw()

#Put in exons? Though coordinates not 100 % clear
#Lollipop plot?


ggplot(x, aes(x=coordinates, y=n, colour=matrix(in_exon_vector)))+
  geom_point(size=3)+
  geom_segment(aes(x=coordinates, xend=coordinates, y=0, yend=n))+
  xlim(c(151690496, 152103274))+
  labs(title="Variants found in ESR1 for more than 5% of samples", x="Coordinates", y="Number of samples with variant", color="Variant Position")+
  scale_color_manual(values=c("coral", "cornflowerblue"))



#ggplot(x, aes(x=coordinates, y=n))+ geom_bar(stat="identity")
