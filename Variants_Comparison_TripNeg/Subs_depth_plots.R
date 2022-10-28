#28.10.22
#Evaluation of sequencing depth table
#TNBC comparison

library(tidyverse)

setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis\\R and outputs\\TNBC_comparison")


subs_depth_df<- read.delim("subs_depth.tsv", sep=" ")
subs_depth_df


subs_depth_df %>% 
  select(-cmd_args, -chr, -Start, -Stop, -full_Path, -Location, -row_id) -> subs_depth

#Boxplot FPKM>10
subs_depth %>%
  #(depth<=75) %>%
  #filter(FPKM>10) %>% 
  ggplot(aes(x=Match, y=depth, fill=Match))+
  geom_boxplot()+
  #geom_jitter(color="black", size=0.4, alpha=0.9)+
  xlab("Matches between DNA data and RNA data.")+
  ylab("Read depth")+
  ylim(0,75)+
  ggtitle("Sequencing depth of substitutions with FPKM >10 in RNA data")+
  scale_fill_discrete(name = "Substitutions found in", labels = c("DNA", "DNA + RNA filtered", "DNA + RNA unfiltered"))+
  theme(axis.text.x = element_blank())
  

#Scatterplot

subs_depth %>% 
  ggplot(aes(x=FPKM, y=depth, color=Match))+
  geom_point()+
  xlim(0,100)+
  ylim(0,100)
#Pretty but way too many points. ._.

  

