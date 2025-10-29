
library(compositions)
library(igraph)
library(qgraph)
library(phangorn)

setwd("G:\\Networks based on genetic distance clusters")

#reading and preparing data
data<-read.table("16S_top_90_ASV.txt",header=TRUE,sep="\t") #16S_top_90_ASV.txt - for bacteria, 18S_top_90_ASV.txt - for microeukaryotes
data<-data[,-ncol(data)]
rownames(data)<-data[,1]
data<-data[,-1]

ASV<-data

#normalization - conversion to relative values
ASV_relative<-ASV
S<-colSums(ASV_relative)
for(i in 1:ncol(ASV_relative)) ASV_relative[,i]<-ASV_relative[,i]/S[i]
colSums(ASV_relative)

#clr transformation
ASV_clr<-ASV_relative
for(i in 1:ncol(ASV_clr)) ASV_clr[,i]<-clr(ASV_clr[,i])

ASV_relative_total<-as.matrix(ASV_relative)

ASV_clr_total<-as.matrix(ASV_clr)

#reading data on environmental parameter

data<-read.table("Expl.txt",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-1]
Expl<-t(data)



#clustering by synergistic distances, distribution of ASV into clusters with distance d<=0.1

pas<-read.dna("16S_aligned.fasta", format="fasta") #16S_aligned.fasta - for bacteria, 18S_aligned.fasta - for microeukaryotes
d<-dist.dna(pas, model = "raw")
fit<-hclust(d, method = "average")
k<-cutree(fit, h=0.10)
data_clust<-data.frame(ASV_id=names(k), cluster=unname(k))
data_clust<-data_clust[order(data_clust$cluster, decreasing=T),]
k<-table(data_clust$cluster)
k<-k[k>=6]
k<-as.numeric(names(k))
n<-which(data_clust$cluster %in% k)
data_clust<-data_clust[n,]

k<-unique(data_clust$cluster)
length(k)

nk=5 #enter cluster number

n<-which(rownames(ASV_clr_total) %in% data_clust[data_clust$cluster==k[nk],]$ASV_id)

ASV_clr<-ASV_clr_total[n,]

ASV_relative<-ASV_relative_total[n,]

Paize<-rowMeans(ASV_relative)

num<-nrow(ASV_clr)

#creation and recording to disk of a distance matrix for a cluster d<=0.1 

dm<-as.matrix(d)
n<-which(rownames(dm) %in% data_clust[data_clust$cluster==k[nk],]$ASV_id)
dm<-dm[n,n]
write.table(cbind(ASV_id=rownames(dm), dm), paste0("dist_matrix_10%_cluster_", nk, ".tsv"), quote=F, sep="\t", row.names=F)


#construction of a correlation matrix

ASV_clr<-rbind(ASV_clr, Expl)

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
	 if(DP[i,j]>0.05) DD[i,j]=0
     if(i==j) DD[i,j]=1
    }
 }

#writing the distance matrix for a cluster to disk
write.table(cbind(ASV_id=rownames(DD), DD), paste0("cor_matrix_10%_cluster_", nk, ".tsv"), quote=F, sep="\t", row.names=F)


group<-list(c(1:num), c((num+1):(num+6))) 
col=c("white", "gray88") 
sh<-c(rep("circle",num), rep("square",6))
Ps<-c(sqrt(2000*unname(Paize)), rep(3, 6))
Ps

lab<-c(rownames(ASV_clr))

#reconstruction and display of the network of pairwise correlation relationships
net<-qgraph(DD, minimum=0, layout="spring", label.cex=1.0, vsize=Ps, labels=lab, groups=group, color=col, shape=sh, border.width=2.5, fade=F, legend=F, esize=3.5)

#dissects the characteristics of the network graph of pairwise correlation relationships

DSIGN<-sign(DD) 
DDP<-sqrt(DD^2) 
net<-qgraph(DDP, minimum=0, layout="spring", label.cex=2.5, vsize=6, labels=lab) 
net<-as.igraph(net) 

nsv<-rowSums(abs(DSIGN))

DNEG<-DSIGN
DNEG[which(DNEG==1)]<-0
n_neg<-rowSums(abs(DNEG)) 

DPOS<-DSIGN
DPOS[which(DPOS==-1)]<-0
n_poz<-rowSums(abs(DPOS)) 

#visualization of a table with graph characteristics - with relationships
data.frame(var=lab, betweenness=round(betweenness(net, normalized=T), digits=3), nsv=nsv, n_poz=n_poz, n_neg=n_neg)


