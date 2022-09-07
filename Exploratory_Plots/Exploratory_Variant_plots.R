#Exploratory plots
#To investigate the reliability of our genotype predictions.

#Want to make a heatmap to see the density of data points for correlation between
#Allele fraction and total number of reads.

library(tidyverse)
setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis")

variants<-read_tsv("AF_TR_out.txt")

#Plotting Allele Frequency------------------------
setwd("C:\\Users\\mirja\\Documents\\University\\Master Thesis")

variants %>% 
  ggplot(aes(x=DP, y=AF) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
  scale_fill_continuous(trans = 'log', low="#0000FF", high="#FFFF00", limits=c(2e-15, 1), na.value="#0000FF")

#not working.

#Lets restrict to one chromosome

variants %>% 
  ggplot(aes(x=DP, y=AF))+ geom_point()

#Can we do a heatmap of that?

#custom_breaks <- function(limits){
#  print(limits)
#  breaks<-abs(log10(limits))
#  return(breaks)
#}


variants %>%
  filter(chrom=="chr6") %>% 
  ggplot(aes(x=DP, y=AF) ) +
  geom_density_2d_filled()#+
  #scale_fill_gradient(trans="log")
  #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  #scale_fill_manual(limits=c(0,0.12), breaks=function(lims){seq(abs(log10(lims[1])), abs(log10(lims[2])))})
  #scale_fill_viridis_c()


#if we exclude the top right corner, we sort of can. 
#Maybe with a 2d histogram

variants %>% 
  #filter(chrom=="chr6") %>% 
  ggplot(aes(x=DP, y=AF))+ geom_bin2d( bins=100)+
  scale_fill_gradient(trans = 'log', low="#0000FF", high="#FFFF00")

#Close enough.


#Lets plot just AF
variants %>% 
  ggplot(aes(x=AF) ) +
  geom_histogram(bins=50)

# Plotting Variant Counts --------------------------

#Only for S00001
gene<-read_tsv("Data_Plots/gene.tsv")
counts<-read_tsv("Data_Plots/S00001_out.tsv")

gene %>% 
  full_join(counts, by=c("Gene Name"="Gene_ID")) %>% 
  mutate(gene_length=End-Start) %>% 
  mutate(counts_kb=Counts*1000/gene_length) %>% 
  filter(FPKM>0) %>% 
  filter(counts_kb<250) %>% 
  ggplot(aes(x=log2(FPKM+1), y=counts_kb))+geom_bin2d(bins=100)+
  scale_fill_gradient(low="#0000FF", high="#FFFF00")

#Maybe requires some zoom in.

gene %>% 
  full_join(counts, by=c("Gene Name"="Gene_ID")) %>% 
  mutate(gene_length=End-Start) %>% 
  mutate(counts_kb=Counts*1000/gene_length) %>% 
  filter(FPKM>0) %>% 
  ggplot(aes(x=log2(FPKM+1), y=counts_kb))+geom_point()

#Outliers need to be explained but we already looked into that.

gene %>% 
  full_join(counts, by=c("Gene Name"="Gene_ID")) %>% 
  mutate(gene_length=End-Start) %>% 
  mutate(counts_kb=Counts*1000/gene_length) %>% 
  ggplot(aes(x=Coverage, y=counts_kb))+geom_point()

gene %>% 
  full_join(counts, by=c("Gene Name"="Gene_ID")) %>% 
  mutate(gene_length=End-Start) %>% 
  mutate(counts_kb=Counts*1000/gene_length) %>% 
  ggplot(aes(x=FPKM, y=counts_kb))+geom_point()

#Those two basically produce the same plot. ._.


#Only for ESR1

ESR1<-read_tsv("Data_Plots/ESR1_out.txt")
str(ESR1)

ESR1 %>% 
  ggplot(aes(x=log2(FPKM+1), y=Count))+geom_point()
#This could be linear, but its definitely monotone, which means we can run a statistical test!

library(lmtest)
library(lmodel2)

#Lets check the assumptions for a Pearson correlation test:

#Normality of data:

#Make sure normal curve is calculated based on correct number of bins in both plots.
x <- ESR1$Count
h <- hist(x, col="#FF99FF")
#making a normal curve with values ranging in data values, to see how it matches.
xfit <- seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE),length=100)
yfit <- dnorm(xfit,mean=mean(x, na.rm=TRUE),sd=sd(x, na.rm=TRUE))
yfit <- yfit*max(h$counts)/max(yfit)
lines(xfit, yfit, lwd=2)
#matches quite well!

shapiro.test(ESR1$Count)
#significant

#kolmogorov:
x <- ESR1$Count
y <- pnorm(summary(x), mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
ks.test(x, y)
#significant, very close to 0.05

qqnorm(ESR1$Count)
qqline(ESR1$Count)
#looks good. Has a tail that leaves the line, but thats common.

#Make sure normal curve is calculated based on correct number of bins in both plots.
x <- ESR1$FPKM
h <- hist(x, col="#FF99FF")
#making a normal curve with values ranging in data values, to see how it matches.
xfit <- seq(min(x, na.rm=TRUE),max(x, na.rm=TRUE),length=100)
yfit <- dnorm(xfit,mean=mean(x, na.rm=TRUE),sd=sd(x, na.rm=TRUE))
yfit <- yfit*max(h$counts)/max(yfit)
lines(xfit, yfit, lwd=2)
#not great but ok.

shapiro.test(ESR1$FPKM)
#significant

#kolmogorov:
x <- ESR1$FPKM
y <- pnorm(summary(x), mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
ks.test(x, y)
#significant

qqnorm(ESR1$FPKM)
qqline(ESR1$FPKM)
#also not great. but we try. Has a tail that leaves the line, but thats common.



cor.test(ESR1$Count, ESR1$FPKM)
#significant: means true correlation is not equal to 0. (based on p value)