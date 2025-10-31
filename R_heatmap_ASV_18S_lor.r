library(vegan)
library(tibble)
library(gplots)
library(pheatmap)


setwd("C:\\Heat map")


#Heat map with clustering of samples and ASV representation based on Bray-Curtis distances
#loading a table with environment variables

data<-read.table("Expl.txt",header=TRUE,sep="\t")
data<-data[,-ncol(data)]
rownames(data)<-data[,1]
data<-data[,-1]
data<-data[, c(1,2)]
Expl<-data


#reading and preparing data
data<-read.table("18S_top_90_ASV.txt",header=TRUE,sep="\t")
Expl_row<-data[,c(1, ncol(data)-1, ncol(data))]
rownames(Expl_row)<-Expl_row[,1]
Expl_row<-Expl_row[,-1, drop=F]
data<-data[, -c(ncol(data)-1, ncol(data))]
rownames(data)<-data[,1]
data<-data[,-1]

Resp<-data

#normalization - conversion to relative values and logarithmic transformation
Resp_p<-Resp
S<-colSums(Resp_p)
for(i in 1:ncol(Resp_p)) Resp_p[,i]<-Resp_p[,i]/S[i]
colSums(Resp_p)

Resp_log<-log10(Resp+1)


#cluster column wiuth the Bray - Curtis distance
d<-vegdist(t(Resp_log), method="bray")
tree_col<-hclust(d, method="complete")
plot(tree_col)


#cluster row  wiuth the Bray - Curtis distance
d<-vegdist(Resp_log, method="bray")
tree_row<-hclust(d, method="complete")
plot(tree_row)


#Select a color to colorize the categorical affiliation of columns and bars on a heat map

col_row_Taxa=c(Ciliophora="midnightblue", Dinophyceae="dodgerblue3", Perkinsea="deepskyblue", Chrysophyceae="gold", Bacillariophyta="red", Dictyochophyceae="red4", Eustigmatophyceae="aquamarine", Synurophyceae="bisque", MAST2="darkviolet", Cryptophyta="mediumseagreen", Haptophyta="forestgreen", Katablepharidophyta="darkgoldenrod3", Telonemia="gray0", Fungi="violet", Choanoflagellida="violetred1", Chlorophyta="lawngreen", Cercozoa="gray47", Conosa="gray")      

Seasonal=c(spring_under_ice="blue", spring_open_water="deepskyblue", summer="green", autumn="gold")

Basin=c(Middle="darkgray", South="paleturquoise1")

col_row_Function=c(Autotrophs="green", Mixotrophs="orange", Heterotrophs="cornflowerblue", Parasits="deeppink")

annoCol=list(Taxa=col_row_Taxa, Function=col_row_Function,Seasonal=Seasonal, Basin=Basin)


#building a heatmap
pheatmap(Resp_log, cluster_rows=tree_row, cluster_cols=tree_col, fontsize=5, 
 annotation_row=Expl_row, annotation_col=Expl, annotation_colors=annoCol, cutree_rows=1, cutree_cols=4)


 
#///////////////////////////////////////////////////////

#Heat map with clustering of samples based on Bray-Curtis distances and clustering ASV based on phylogenetic tree (IQTREE)


#cluster row   based on phylogenetic tree (IQTREE)
tree<-read.tree("18S_clad.tree") #reading a IQTREE phylogenetic tree from disk
tree_row<-as.hclust.phylo(tree)
plot(tree_row)

st<-tree_row$labels
Resp_log<-Resp_log[order(match(rownames(Resp_log),st)),]
Expl_row<-Expl_row[order(match(rownames(Expl_row),st)), ,drop=F]


#Select a color to colorize the categorical affiliation of columns and bars on a heat map

col_row_Taxa=c(Ciliophora="midnightblue", Dinophyceae="dodgerblue3", Perkinsea="deepskyblue", Chrysophyceae="gold", Bacillariophyta="red", Dictyochophyceae="red4", Eustigmatophyceae="aquamarine", Synurophyceae="bisque", MAST2="darkviolet", Cryptophyta="mediumseagreen", Haptophyta="forestgreen", Katablepharidophyta="darkgoldenrod3", Telonemia="gray0", Fungi="violet", Choanoflagellida="violetred1", Chlorophyta="lawngreen", Cercozoa="gray47", Conosa="gray")      

Seasonal=c(spring_under_ice="blue", spring_open_water="deepskyblue", summer="green", autumn="gold")

Basin=c(Middle="darkgray", South="paleturquoise1")

col_row_Function=c(Autotrophs="green", Mixotrophs="orange", Heterotrophs="cornflowerblue", Parasits="deeppink")

annoCol=list(Taxa=col_row_Taxa, Function=col_row_Function,Seasonal=Seasonal, Basin=Basin)

#building a heatmap
pheatmap(Resp_log, cluster_rows=tree_row, cluster_cols=tree_col, fontsize=5, 
 annotation_row=Expl_row, annotation_col=Expl, annotation_colors=annoCol, cutree_rows=1, cutree_cols=4)

