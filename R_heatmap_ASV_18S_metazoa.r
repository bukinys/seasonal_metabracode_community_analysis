library(vegan)
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
data<-read.table("18S_metazoa_ASV.txt",header=TRUE,sep="\t")
Expl_row<-data[,c(1, ncol(data))]
rownames(Expl_row)<-Expl_row[,1]
Expl_row<-Expl_row[,-1, drop=F]
data<-data[,-ncol(data)]
rownames(data)<-data[,1]
data<-data[,-1]

Resp<-data

#normalization - conversion to relative values and logarithmic transformation
Resp_norm<-log10(Resp+1)

#cluster column wiuth the Bray - Curtis distance
d<-vegdist(t(Resp_norm), method="bray")
tree_col<-hclust(d, method="complete")
plot(tree_col)


#cluster row  wiuth the Bray - Curtis distance
d<-vegdist(Resp_norm, method="bray")
tree_row<-hclust(d, method="complete")
plot(tree_row)


#Select a color to colorize the categorical affiliation of columns and bars on a heat map

Taxa=c(Copepoda="midnightblue", Amphipoda="aquamarine", Cladocera="gold", Rotifera="red", Nematoda="bisque", Annelida="mediumseagreen")

Seasonal=c(spring_under_ice="blue", spring_open_water="deepskyblue", summer="green", autumn="gold")

Basin=c(Middle="darkgray", South="paleturquoise1")

annoCol=list(Taxa=Taxa, Seasonal=Seasonal, Basin=Basin)

#building a heatmap
pheatmap(Resp_norm, cluster_rows=tree_row, cluster_cols=tree_col, fontsize=5, 
 annotation_row=Expl_row, annotation_col=Expl, annotation_colors=annoCol, cutree_rows=1, cutree_cols=1)



#///////////////////////////////////////////////////////

#Heat map with clustering of samples based on Bray-Curtis distances and clustering ASV based on phylogenetic tree (IQTREE)


#cluster row   based on phylogenetic tree (IQTREE)
tree<-read.tree("18S_metazoa.tree") #reading a IQTREE metazoa phylogenetic tree from disk
tree_row<-as.hclust.phylo(tree)
plot(tree_row)

st<-tree_row$labels
Resp_norm<-Resp_norm[order(match(rownames(Resp_norm),st)),]
Expl_row<-Expl_row[order(match(rownames(Expl_row),st)), ,drop=F]


#Select a color to colorize the categorical affiliation of columns and bars on a heat map

Taxa=c(Copepoda="midnightblue", Amphipoda="aquamarine", Cladocera="gold", Rotifera="red", Nematoda="bisque", Annelida="mediumseagreen")

Seasonal=c(spring_under_ice="blue", spring_open_water="deepskyblue", summer="green", autumn="gold")

Basin=c(Middle="darkgray", South="paleturquoise1")

annoCol=list(Taxa=Taxa, Seasonal=Seasonal, Basin=Basin)

#building a heatmap
pheatmap(Resp_norm, cluster_rows=tree_row, cluster_cols=tree_col, fontsize=5, 
 annotation_row=Expl_row, annotation_col=Expl, annotation_colors=annoCol, cutree_rows=1, cutree_cols=1)


