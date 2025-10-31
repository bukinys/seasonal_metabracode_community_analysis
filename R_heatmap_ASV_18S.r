library(vegan)
library(tibble)
library(gplots)
library(pheatmap)


setwd("C:\\COR_heatmap")

#loading a table with environment variables

data<-read.table("Expl.txt",header=TRUE,sep="\t")
rownames(data)<-data[,1]
data<-data[,-1]
data<-data[, -c(1,2)]
Expl<-t(data)


#reading and preparing data
data<-read.table("18S_top_90_ASV.txt",header=TRUE,sep="\t")
Expl_row<-data[,c(1, ncol(data)-1, ncol(data))]
rownames(Expl_row)<-Expl_row[,1]
Expl_row<-Expl_row[,-1, drop=F]
data<-data[, -c(ncol(data)-1, ncol(data))]
rownames(data)<-data[,1]
data<-data[,-1]
Resp<-as.matrix(data)

#normalization - conversion to relative values
Resp_relative<-Resp
S<-colSums(Resp_relative)
for(i in 1:ncol(Resp_relative)) Resp_relative[,i]<-Resp_relative[,i]/S[i]
colSums(Resp_relative)
Resp<-Resp_relative


#construction of a correlation matrix
DD=matrix(nrow=nrow(Resp), ncol=nrow(Expl))

rownames(DD)<-rownames(Resp)
colnames(DD)<-rownames(Expl)

DP<-DD

for(i in 1:nrow(Resp))
 {
   for(j in 1:nrow(Expl))
    {
     R=cor.test(Resp[i,], Expl[j,], method ="spearman")
     DD[i,j]=R$estimate
	 DP[i,j]<-R$p.value
    }
 }

DD[sqrt(DD^2)<0.5]<-0

write.table(cbind(sample_id=rownames(DD), DD), "r_cor_18S.tsv", sep="\t")

#cluster column wiuth the Euclidean distance
d<-vegdist(t(DD), method="euclidean")
tree_col<-hclust(d, method="complete")
plot(tree_col)
#cluster row  wiuth the Euclidean distance
d<-vegdist(DD, method="euclidean")
tree_row<-hclust(d, method="complete")
plot(tree_row)


#Select a color to colorize the categorical affiliation of columns and bars on a heat map

col_row_Taxa=c(Ciliophora="midnightblue", Dinophyceae="dodgerblue3", Perkinsea="deepskyblue", Chrysophyceae="gold", Bacillariophyta="red", Dictyochophyceae="red4", Eustigmatophyceae="aquamarine", Synurophyceae="bisque", MAST2="darkviolet", Cryptophyta="mediumseagreen", Haptophyta="forestgreen", Katablepharidophyta="darkgoldenrod3", Telonemia="gray0", Fungi="violet", Choanoflagellida="violetred1", Chlorophyta="lawngreen", Cercozoa="gray47", Conosa="gray")      

col_row_Function=c(Autotrophs="green", Mixotrophs="blue", Heterotrophs="brown1", Parasits="gray36")

annoCol=list(Taxa=col_row_Taxa, Function=col_row_Function)

#building a heatmap
pheatmap(DD, cluster_rows=tree_row, cluster_cols=tree_col, fontsize=5, 
 annotation_row=Expl_row, annotation_colors=annoCol, cutree_rows=3, cutree_cols=2)
