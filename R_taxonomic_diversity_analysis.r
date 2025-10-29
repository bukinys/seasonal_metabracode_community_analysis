
library(vegan)

setwd("G:\\Taxonomic_diversity")

#reading the genotype counte table
data<-read.csv("ASV_16S.txt", header=TRUE, sep="\t", check.names = FALSE)

data_count<-data[, -ncol(data)]
rownames(data_count)<-data_count[,1]
data_count<-data_count[,-1]

d_f<-as.data.frame(t(data_count))

#plotting rarefaction curves
rc<-rarecurve(d_f, step = 20, sample=70, col = "blue", cex=0.6)


#calculation of the expected number of genotypes - Chao1 and ASE, calculation of the Shannon and Simpson index
diversiti<-t(estimateR(d_f))
shannon<-diversity(d_f, index = "shannon")
simpson<-diversity(d_f, index = "simpson")
diversiti<-cbind(sample_id=rownames(diversiti), reads_number=rowSums(d_f),diversiti,  shannon=shannon, simpson=simpson)
diversiti<-as.data.frame(diversiti)

diversiti #visualization of a table with taxonomic diversity indicators


