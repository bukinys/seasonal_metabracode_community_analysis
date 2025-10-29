
library(compositions)
library(igraph)
library(qgraph)

setwd("G:\\Networks total")

#reconstruction of a network with positive correlation coefficients
#reading and preparing data
data<-read.table("16S_top_90_ASV.txt",header=TRUE,sep="\t") #16S_top_90_ASV.txt - for bacteria, 18S_top_90_ASV.txt - for microeukaryotes
rownames(data)<-data[,1]
data<-data[,-1]

ASV<-data

#normalization - conversion to relative values
ASV_relative<-ASV
S<-colSums(ASV_relative)
for(i in 1:ncol(ASV_relative)) ASV_relative[,i]<-ASV_relative[,i]/S[i]
colSums(ASV_relative)
Psize<-rowMeans(ASV_relative)

#clr transformation
ASV_clr<-ASV_relative
for(i in 1:ncol(ASV_clr)) ASV_clr[,i]<-clr(ASV_clr[,i])
ASV_relative_total<-as.matrix(ASV_relative)
ASV_clr<-as.matrix(ASV_clr)
N_ASV<-nrow(ASV_clr)

ASV_clr<-rbind(ASV_clr)

#construction of a correlation matrix

DD=matrix(nrow=nrow(ASV_clr), ncol=nrow(ASV_clr))
rownames(DD)<-rownames(ASV_clr)
colnames(DD)<-rownames(ASV_clr)

DP<-DD

for(i in 1:nrow(ASV_clr))
 {
   for(j in 1:nrow(ASV_clr))
    {
     R=cor.test(ASV_clr[i,], ASV_clr[j,], method ="spearman")
     DD[i,j]=R$estimate
	 DP[i,j]<-R$p.value
	 DP[i,j]<-p.adjust(R$p.value, method="BH", n=nrow(ASV_clr))
	 if(sqrt(DD[i,j]^2)<=0.65) DD[i,j]=0
     if(i==j) DD[i,j]=0
    }
 }

for(i in 1:N_ASV)
 {
   for(j in 1:N_ASV)
    {
     if(DD[i,j]<0) DD[i,j]=0
    }
 }

DDV<-sqrt(DD^2)
S<-rowSums(DDV)
n<-which(S==0)

lab<-c(rownames(ASV_clr))

if(length(n)>=1)
 {
  DD<-DD[-n, -n]
  Psize<-Psize[-n]
  N_ASV<-N_ASV-length(n)
  lab<-lab[-n]
 } 

#writing the matrix of pairwise correlation to disk
write.table(cbind(ASV_id=rownames(DD), DD), paste0("16S_cor_matrix_0.65_positive", ".tsv"), quote=F, sep="\t", row.names=F)

group<-list(c(1:N_ASV)) 
col=c("white") 
sh<-c(rep("circle",N_ASV)) 
Ps<-c(sqrt(1000*unname(c(Psize))))
Ps

#reconstruction and display of the network of pairwise correlation relationships
net<-qgraph(DD, minimum=0, layout="spring", label.cex=3.0, vsize=Ps, labels=lab, groups=group, color=col, shape=sh, border.width=2.5, fade=F, legend=F, esize=3.5)


#reconstruction of a network with negative correlation coefficients
#reading and preparing data
data<-read.table("16S_top_90_ASV.txt",header=TRUE,sep="\t") #16S_top_90_ASV.txt - for bacteria, 18S_top_90_ASV.txt - for microeukaryotes
rownames(data)<-data[,1]
data<-data[,-1]

ASV<-data

#normalization - conversion to relative values
ASV_relative<-ASV
S<-colSums(ASV_relative)
for(i in 1:ncol(ASV_relative)) ASV_relative[,i]<-ASV_relative[,i]/S[i]
colSums(ASV_relative)
Psize<-rowMeans(ASV_relative)

#clr transformation
ASV_clr<-ASV_relative
for(i in 1:ncol(ASV_clr)) ASV_clr[,i]<-clr(ASV_clr[,i])
ASV_relative_total<-as.matrix(ASV_relative)
ASV_clr<-as.matrix(ASV_clr)
N_ASV<-nrow(ASV_clr)

ASV_clr<-rbind(ASV_clr)

#construction of a correlation matrix

DD=matrix(nrow=nrow(ASV_clr), ncol=nrow(ASV_clr))
rownames(DD)<-rownames(ASV_clr)
colnames(DD)<-rownames(ASV_clr)

DP<-DD

for(i in 1:nrow(ASV_clr))
 {
   for(j in 1:nrow(ASV_clr))
    {
     R=cor.test(ASV_clr[i,], ASV_clr[j,], method ="spearman")
     DD[i,j]=R$estimate
	 DP[i,j]<-R$p.value
	 DP[i,j]<-p.adjust(R$p.value, method="BH", n=nrow(ASV_clr))
	 if(sqrt(DD[i,j]^2)<=0.65) DD[i,j]=0
     if(i==j) DD[i,j]=0
    }
 }

for(i in 1:N_ASV)
 {
   for(j in 1:N_ASV)
    {
     if(DD[i,j]>0) DD[i,j]=0
    }
 }

DDV<-sqrt(DD^2)
S<-rowSums(DDV)
n<-which(S==0)

lab<-c(rownames(ASV_clr))

if(length(n)>=1)
 {
  DD<-DD[-n, -n]
  Psize<-Psize[-n]
  N_ASV<-N_ASV-length(n)
  lab<-lab[-n]
 } 
 
#writing the matrix of pairwise correlation to disk
write.table(cbind(ASV_id=rownames(DD), DD), paste0("16S_cor_matrix_0.65_negative", ".tsv"), quote=F, sep="\t", row.names=F)

group<-list(c(1:N_ASV)) 
col=c("white") 
sh<-c(rep("circle",N_ASV)) 
Ps<-c(sqrt(1000*unname(c(Psize))))
Ps

#reconstruction and display of the network of pairwise correlation relationships
net<-qgraph(DD, minimum=0, layout="spring", label.cex=3.0, vsize=Ps, labels=lab, groups=group, color=col, shape=sh, border.width=2.5, fade=F, legend=F, esize=3.5)


#//////////////////////////////////////////////////////////
#reconstruction of a network correlation relationships with environmental variables

#reading and preparing data
data<-read.table("16S_top_90_ASV.txt",header=TRUE,sep="\t") #16S_aligned.fasta - for bacteria, 18S_aligned.fasta - for microeukaryotes
rownames(data)<-data[,1]
data<-data[,-1]

ASV<-data

#normalization - conversion to relative values
ASV_relative<-ASV
S<-colSums(ASV_relative)
for(i in 1:ncol(ASV_relative)) ASV_relative[,i]<-ASV_relative[,i]/S[i]
colSums(ASV_relative)
Psize<-rowMeans(ASV_relative)

#clr transformation
ASV_clr<-ASV_relative
for(i in 1:ncol(ASV_clr)) ASV_clr[,i]<-clr(ASV_clr[,i])
ASV_relative_total<-as.matrix(ASV_relative)
ASV_clr<-as.matrix(ASV_clr)
N_ASV<-nrow(ASV_clr)

#reading data about environment variables
data<-read.table("Expl.tsv",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-1]
Expl<-t(data)

#obtaining a complete data array
ASV_clr<-rbind(ASV_clr, Expl)

#construction of a correlation matrix

DD=matrix(nrow=nrow(ASV_clr), ncol=nrow(ASV_clr))
rownames(DD)<-rownames(ASV_clr)
colnames(DD)<-rownames(ASV_clr)
DP<-DD
for(i in 1:nrow(ASV_clr))
 {
   for(j in 1:nrow(ASV_clr))
    {
     R=cor.test(ASV_clr[i,], ASV_clr[j,], method ="spearman")
     DD[i,j]=R$estimate
	 DP[i,j]<-p.adjust(R$p.value, method="BH", n=nrow(ASV_clr))
	 if(sqrt(DD[i,j]^2)<=0.65) DD[i,j]=0
     if(i==j) DD[i,j]=0
    }
 }

for(i in 1:N_ASV)
 {
   for(j in 1:N_ASV)
    {
     DD[i,j]=0
    }
 }


DDV<-sqrt(DD^2)
S<-rowSums(DDV)
n<-which(S==0)

lab<-c(rownames(ASV_clr))

if(length(n)>=1)
 {
  DD<-DD[-n, -n]
  Psize<-Psize[-n]
  N_ASV<-N_ASV-length(n)
  lab<-lab[-n]
 } 

#writing the matrix of pairwise correlation for the cluster to disk
write.table(cbind(ASV_id=rownames(DD), DD), paste0("16S-Env_cor_matrix_0.65_environment", ".tsv"), quote=F, sep="\t", row.names=F)


group<-list(c(1:N_ASV), c((N_ASV+1):(N_ASV+7)))
col=c("white", "gray88") 
sh<-c(rep("circle",N_ASV), rep("square",7))
Ps<-c(sqrt(1000*unname(c(Psize))), rep(3, 7))
Ps

#reconstruction and display of the network of pairwise correlation relationships
net<-qgraph(DD, minimum=0, layout="spring", label.cex=3.0, vsize=Ps, labels=lab, groups=group, color=col, shape=sh, border.width=2.5, fade=F, legend=F, esize=3.5)
