
library(vegan)

setwd("G:\\PERMANOVA_CCA")

#data preparation

data<-read.table("ASV_16S.txt",header=TRUE,sep="\t") #loading ASV table (ASV_16S.txt, ASV_18S.txt or 18S_metazoa_ASV.txt)
data<-data[,-ncol(data)]
rownames(data)<-data[,1]
data<-data[,-1]

S<-colSums(data)
Sm<-mean(S)
for(i in 1:ncol(data)) data[,i]<-ceiling(Sm*data[,i]/S[i])
colSums(data)

#writing a normalized ASV table to disk
write.table(cbind(ASV_id=rownames(data), data), "ASV_16S_norm.txt", row.names=F, col.names=TRUE, sep = "\t")

Resp<-as.data.frame(t(data))

data<-read.table("Expl.tsv",header=TRUE,sep="\t") #loading a table with environment variables
rownames(data)<-data[,1]
data<-data[,-1]
Expl<-data

rownames(Resp)<-rownames(Expl)

#PERMANOVA analysis

rez_permanova<-data.frame(factor=colnames(Expl), per_R2=rep(0, 11), per_P_value=rep(0, 11))

for(i in 1:11)
 {
  fit<-adonis(Resp ~ (Expl[, i, drop=T]), method="jaccard", permutations=1000)
  rez_permanova$per_R2[i]<-fit$aov.tab[,5][1]
  rez_permanova$per_P_value[i]<-fit$aov.tab[,6][1]
 } 

rez_permanova #output table with PERMANOVA results

#anova cca analysis

rez_anova_cca<-data.frame(factor=colnames(Expl), F_value=rep(0, 11), P_value=rep(0, 11))

for(i in 1:11)
 {
  ord<-cca(formula = Resp ~ (Expl[, i, drop=T]))
  fit<-anova(ord, by="term", permutations=1000)
  rez_anova_cca$F_value[i]<-fit[,3][1]
  rez_anova_cca$P_value[i]<-fit[,4][1]
 } 

rez_anova_cca #output table with anova cca results

#/////////////////////////////////////////
#Data normalization by ranking from 0 to 1

Resp<-decostand(Resp, method="rang", MARGIN=2)

#CCA analysis

ord<-cca(Resp ~ Basin, data=Expl)

#CCA plot

plot(ord)

plot(ord, typ="n", xlim=c(-1.0, 2), ylim=c(-1.0, 2.5))

points(ord, display="sites", pch=1, lwd=3.5, cex = 1.7, col="blue", select=which(Expl$Basin=="South" & Expl$Seasonal=="spring_under_ice"))
points(ord, display="sites", pch=1, lwd=3.5, cex = 1.7, col="green", select=which(Expl$Basin=="South" & Expl$Seasonal=="spring_open_water"))
points(ord, display="sites", pch=1, lwd=3.5, cex = 1.7, col="red", select=which(Expl$Basin=="South" & Expl$Seasonal=="summer"))
points(ord, display="sites", pch=1, lwd=3.5, cex = 1.7, col="black", select=which(Expl$Basin=="South" & Expl$Seasonal=="autumn"))

points(ord, display="sites", pch=0, lwd=3.5, cex = 1.7, col="blue", select=which(Expl$Basin=="Middle" & Expl$Seasonal=="spring_under_ice"))
points(ord, display="sites", pch=0, lwd=3.5, cex = 1.7, col="green", select=which(Expl$Basin=="Middle" & Expl$Seasonal=="spring_open_water"))
points(ord, display="sites", pch=0, lwd=3.5, cex = 1.7, col="red", select=which(Expl$Basin=="Middle" & Expl$Seasonal=="summer"))
points(ord, display="sites", pch=0, lwd=3.5, cex = 1.7, col="black", select=which(Expl$Basin=="Middle" & Expl$Seasonal=="autumn"))

text(ord, display = "sites")


st<-paste("ord ~ ", colnames(Expl)[1], sep="")
for(i in 2:ncol(Expl)) st<-paste(st,colnames(Expl)[i], sep="+")

sample<-rep("0", nrow(Expl))
for(i in 1:nrow(Expl))
 {
  sample[i]<-strsplit(rownames(Expl)[i], "_")[[1]][1]
 }

ord.fit<-envfit(as.formula(st), data=Expl, perm=1000)
ord.fit
plot(ord.fit)





















#///////////////анализ методом главынх компонент усреднение по факторам///////////////////

fit<-prcomp(Resp)

pca_factor<-Expl$Seasonal

fviz_pca_biplot(fit, habillage=pca_factor, addEllipses=TRUE, ellipse.level=0.90)

#//////////////////////////////////////////////////////////////////


d<-vegdist(Resp, method="euclidean")

fit<-cmdscale(d,eig=TRUE, k=2)

x<-fit$points[,1]
y<-fit$points[,2]
plot(x, y, xlim=c(-3, 9), ylim=c(-0.5, 1.2), xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS")
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS")
text(x, y, labels=pca_factor, cex=.7)


st<-paste("fit ~ ", colnames(Expl)[1], sep="")
for(i in 2:ncol(Expl)) st<-paste(st,colnames(Expl)[i], sep="+")

ord.fit<-envfit(as.formula(st), data=Expl, perm=10000)
ord.fit
plot(ord.fit)





ord<-metaMDS(Expl, distance = "euclidean")

st<-paste("ord ~ ", colnames(Expl)[1], sep="")
for(i in 2:ncol(Expl)) st<-paste(st,colnames(Expl)[i], sep="+")

sample<-rep("0", nrow(Expl))
for(i in 1:nrow(Expl))
 {
  sample[i]<-strsplit(rownames(Expl)[i], "_")[[1]][1]
 }

plot(ord, typ="n", xlim=c(-2.0,1), ylim=c(-1,1))

points(ord, display="sites", pch=1, lwd=3.5, cex = 1.7, col="black", select=which(sample=="Olenek"))
points(ord, display="sites", pch=0, lwd=3.5, cex = 1.7, col="red", select=which(sample=="Yana"))
points(ord, display="sites", pch=2, lwd=3.5, cex = 1.7, col="green", select=which(sample=="Indigirka"))
points(ord, display="sites", pch=6, lwd=3.5, cex = 1.7, col="blue", select=which(sample=="Kolyma"))

text(ord, display = "sites")

ord.fit<-envfit(as.formula(st), data=Expl, perm=10000)
ord.fit
plot(ord.fit)