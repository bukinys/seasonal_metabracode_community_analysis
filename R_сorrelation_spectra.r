library(compositions)
library(phangorn)


setwd("G:\\Correlation spectra for genetic clusters")

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

ASV_relative<-as.matrix(ASV_relative)

ASV_clr<-as.matrix(ASV_clr)


#analysis on a scale of relative values
#construction of a correlation matrix

DD=matrix(nrow=nrow(ASV_relative), ncol=nrow(ASV_relative))

rownames(DD)<-rownames(ASV_relative)
colnames(DD)<-rownames(ASV_relative)

DP<-DD

for(i in 1:nrow(ASV_relative))
 {
   for(j in 1:nrow(ASV_relative))
    {
     R=cor.test(ASV_relative[i,], ASV_relative[j,], method ="spearman")
     DD[i,j]=R$estimate
	 DP[i,j]<-p.adjust(R$p.value, method="BH", n=nrow(ASV_clr))
	 if(sqrt(DD[i,j]^2)<=0.53) DD[i,j]=0
     if(i==j) DD[i,j]=1
    }
 }


#Calculation of phylogenetic distances - the proportion of substitutions
pas<-read.dna("16S_aligned.fasta", format="fasta")  #16S_aligned.fasta - for bacteria, 18S_aligned.fasta - for microeukaryotes
DF<-dist.dna(pas, model = "raw")
DF<-as.matrix(DF)

n<-which(rownames(DF) %in% rownames(DD))

DF<-DF[order(match(rownames(DF),rownames(DD))),]
DF<-DF[, order(match(colnames(DF),colnames(DD)))]

rownames(DD)
rownames(DF)

ASV_row<-c()
ASV_col<-c()
ASV_dist<-c()
ASV_cor<-c()

for(i in 2:(nrow(ASV_relative)))
 {
   for(j in 1:(i-1))
    {
	 ASV_row<-c(ASV_row, rownames(DD)[i])
	 ASV_col<-c(ASV_col, colnames(DD)[j])
	 ASV_dist<-c(ASV_dist, DF[i,j])
	 ASV_cor<-c(ASV_cor, DD[i,j])
    }
 }

df_cor<-data.frame(ASV_row, ASV_col, ASV_dist, ASV_cor)

df_cor<-df_cor[order(df_cor$ASV_dist, decreasing=F),]

#recording on disk the table of genetic distances for the ASV and the corresponding correlation coefficients 
write.table(df_cor, "16S_cor_dist_data.tsv", sep="\t", row.names=F)


#cluster data extraction

DD<-as.dist(DD)
DF<-as.dist(DF)

plot(DF, DD)

DR<-as.vector(DD)

DF<-as.vector(DF)

df<-data.frame(DR, DF)

#obtaining a distance distribution table

df3<-df[df$DF<=0.03,]
rp<-nrow(df3[df3$DR>0,])
r0<-nrow(df3[df3$DR==0,])
rm<-nrow(df3[df3$DR<0,])
d_r<-data.frame(r0, rp, rm)

df5<-df[df$DF>0.03 & df$DF<=0.05,]
rp<-nrow(df5[df5$DR>0,])
r0<-nrow(df5[df5$DR==0,])
rm<-nrow(df5[df5$DR<0,])
d_r<-rbind(d_r, data.frame(r0, rp, rm))

df10<-df[df$DF>0.05 & df$DF<=0.1,]
rp<-nrow(df10[df10$DR>0,])
r0<-nrow(df10[df10$DR==0,])
rm<-nrow(df10[df10$DR<0,])
d_r<-rbind(d_r, data.frame(r0, rp, rm))

dft<-df[df$DF>0.1,]
rp<-nrow(dft[dft$DR>0,])
r0<-nrow(dft[dft$DR==0,])
rm<-nrow(dft[dft$DR<0,])
d_r<-rbind(d_r, data.frame(r0, rp, rm))

d_r<-t(d_r)

d_r #table of the number of negative, zero, and positive correlation coefficients

dist_value = c("d ≤ 3%", "3% < d ≤5%", "5% < d ≤ 10%", "d > 10%")

colnames(d_r)<-dist_value

S<-colSums(d_r)

for(i in 1:ncol(d_r)) d_r[,i]<-d_r[,i]/S[i]

d_r #table of proportions of negative, zero, and positive correlation coefficients

#plotting a barplot
colors = c("steelblue1", "red","green")
dist_value = c("d ≤ 3%", "3% < d ≤5%", "5% < d ≤ 10%", "d > 10%")
corvalue=c("r=0","r>0", "r<0")
barplot(as.matrix(d_r), main = "barplot", names.arg=dist_value,space=1, xlab = "dist", ylab = "proportion", col=colors)
legend("topleft", corvalue, cex = 1.3, fill = colors)


#analysis on a scale of clr values
#construction of a correlation matrix

DD=matrix(nrow=nrow(ASV_clr), ncol=nrow(ASV_clr))

rownames(DD)<-rownames(ASV_clr)
colnames(DD)<-rownames(ASV_clr)

DP<-DD

for(i in 1:nrow(ASV_relative))
 {
   for(j in 1:nrow(ASV_relative))
    {
     R=cor.test(ASV_clr[i,], ASV_clr[j,], method ="spearman")
     DD[i,j]=R$estimate
	 DP[i,j]<-p.adjust(R$p.value, method="BH", n=nrow(ASV_clr))
	 if(sqrt(DD[i,j]^2)<=0.53) DD[i,j]=0
     if(i==j) DD[i,j]=1
    }
 }

#Calculation of phylogenetic distances - the proportion of substitutions
pas<-read.dna("16S_aligned.fasta", format="fasta")  #16S_aligned.fasta - for bacteria, 18S_aligned.fasta - for microeukaryotes
DF<-dist.dna(pas, model = "raw")
DF<-as.matrix(DF)

n<-which(rownames(DF) %in% rownames(DD))

DF<-DF[order(match(rownames(DF),rownames(DD))),]
DF<-DF[, order(match(colnames(DF),colnames(DD)))]

rownames(DD)
rownames(DF)


ASV_row<-c()
ASV_col<-c()
ASV_dist<-c()
ASV_cor<-c()

for(i in 2:(nrow(ASV_relative)))
 {
   for(j in 1:(i-1))
    {
	 ASV_row<-c(ASV_row, rownames(DD)[i])
	 ASV_col<-c(ASV_col, colnames(DD)[j])
	 ASV_dist<-c(ASV_dist, DF[i,j])
	 ASV_cor<-c(ASV_cor, DD[i,j])
    }
 }

df_cor<-data.frame(ASV_row, ASV_col, ASV_dist, ASV_cor)

df_cor<-df_cor[order(df_cor$ASV_dist, decreasing=F),]

#recording on disk the table of genetic distances for the ASV and the corresponding correlation coefficients 
write.table(df_cor, "16S_cor_dist_data_clr.tsv", sep="\t", row.names=F)


#cluster data extraction

DD<-as.dist(DD)
DF<-as.dist(DF)

plot(DF, DD)

DR<-as.vector(DD)

DF<-as.vector(DF)

df<-data.frame(DR, DF)

#obtaining a distance distribution table

df3<-df[df$DF<=0.03,]
rp<-nrow(df3[df3$DR>0,])
r0<-nrow(df3[df3$DR==0,])
rm<-nrow(df3[df3$DR<0,])
d_r<-data.frame(r0, rp, rm)

df5<-df[df$DF>0.03 & df$DF<=0.05,]
rp<-nrow(df5[df5$DR>0,])
r0<-nrow(df5[df5$DR==0,])
rm<-nrow(df5[df5$DR<0,])
d_r<-rbind(d_r, data.frame(r0, rp, rm))

df10<-df[df$DF>0.05 & df$DF<=0.1,]
rp<-nrow(df10[df10$DR>0,])
r0<-nrow(df10[df10$DR==0,])
rm<-nrow(df10[df10$DR<0,])
d_r<-rbind(d_r, data.frame(r0, rp, rm))

dft<-df[df$DF>0.1,]
rp<-nrow(dft[dft$DR>0,])
r0<-nrow(dft[dft$DR==0,])
rm<-nrow(dft[dft$DR<0,])
d_r<-rbind(d_r, data.frame(r0, rp, rm))

d_r<-t(d_r)

d_r #table of the number of negative, zero, and positive correlation coefficients

dist_value = c("d ≤ 3%", "3% < d ≤5%", "5% < d ≤ 10%", "d > 10%")

colnames(d_r)<-dist_value

write.table(d_r, "16S_cor_clr_absolute_data.tsv", sep="\t", row.names=F)

S<-colSums(d_r)

for(i in 1:ncol(d_r)) d_r[,i]<-d_r[,i]/S[i]

d_r #table of proportions of negative, zero, and positive correlation coefficients


#plotting a barplot
colors = c("steelblue1", "red","green")
dist_value = c("d ≤ 3%", "3% < d ≤5%", "5% < d ≤ 10%", "d > 10%")
corvalue=c("r=0","r>0", "r<0")
barplot(as.matrix(d_r), main = "barplot", names.arg=dist_value,space=1, xlab = "dist", ylab = "proportion", col=colors)
legend("topleft", corvalue, cex = 1.3, fill = colors)

































#анализ дистанций 3%

df3<-df[df$DF<=0.03,]

rp<-nrow(df3[df3$DR>0,])/nrow(df3)
r0<-nrow(df3[df3$DR==0,])/nrow(df3)
rm<-nrow(df3[df3$DR<0,])/nrow(df3)



dbar<-data.frame(name=c("r>0", "r=0", "r<0"), value=c(rp, r0, rm))

barplot(height=dbar$value, names=dbar$name,  space=1, ylim=c(0, 0.8))

#анализ дистанций 5%

df5<-df[df$DF>0.03 & df$DF<=0.05,]

rp<-nrow(df5[df5$DR>0,])/nrow(df5)
r0<-nrow(df5[df5$DR==0,])/nrow(df5)
rm<-nrow(df5[df5$DR<0,])/nrow(df5)

dbar<-data.frame(name=c("r>0", "r=0", "r<0"), value=c(rp, r0, rm))

barplot(height=dbar$value, names=dbar$name,  space=1, ylim=c(0, 0.8))

#анализ дистанций 10%

df10<-df[df$DF>0.05 & df$DF<=0.1,]

rp<-nrow(df10[df10$DR>0,])/nrow(df10)
r0<-nrow(df10[df10$DR==0,])/nrow(df10)
rm<-nrow(df10[df10$DR<0,])/nrow(df10)

dbar<-data.frame(name=c("r>0", "r=0", "r<0"), value=c(rp, r0, rm))

barplot(height=dbar$value, names=dbar$name, space=1, ylim=c(0, 0.8))

#анализ дистанций 10%

dft<-df[df$DF>0.1,]

rp<-nrow(dft[dft$DR>0,])/nrow(dft)
r0<-nrow(dft[dft$DR==0,])/nrow(dft)
rm<-nrow(dft[dft$DR<0,])/nrow(dft)

dbar<-data.frame(name=c("r>0", "r=0", "r<0"), value=c(rp, r0, rm))

barplot(height=dbar$value, names=dbar$name, space=1, ylim=c(0, 0.8))






























#анализ дистанций 3%

df3<-df[df$DF<=0.03,]
nrow(df3)

rp<-nrow(df3[df3$DR>0,])/nrow(df3)
r0<-nrow(df3[df3$DR==0,])/nrow(df3)
rm<-nrow(df3[df3$DR<0,])/nrow(df3)



dbar<-data.frame(name=c("r>0", "r=0", "r<0"), value=c(rp, r0, rm))

barplot(height=dbar$value, names=dbar$name,  space=1, ylim=c(0, 0.8))

#анализ дистанций 5%

df5<-df[df$DF>0.03 & df$DF<=0.05,]
nrow(df5)

rp<-nrow(df5[df5$DR>0,])/nrow(df5)
r0<-nrow(df5[df5$DR==0,])/nrow(df5)
rm<-nrow(df5[df5$DR<0,])/nrow(df5)

dbar<-data.frame(name=c("r>0", "r=0", "r<0"), value=c(rp, r0, rm))

barplot(height=dbar$value, names=dbar$name,  space=1, ylim=c(0, 0.8))

#анализ дистанций 10%

df10<-df[df$DF>0.05 & df$DF<=0.1,]
nrow(df10)

rp<-nrow(df10[df10$DR>0,])/nrow(df10)
r0<-nrow(df10[df10$DR==0,])/nrow(df10)
rm<-nrow(df10[df10$DR<0,])/nrow(df10)

dbar<-data.frame(name=c("r>0", "r=0", "r<0"), value=c(rp, r0, rm))

barplot(height=dbar$value, names=dbar$name, space=1, ylim=c(0, 0.8))

#анализ дистанций 10%

dft<-df[df$DF>0.1,]
nrow(dft)

rp<-nrow(dft[dft$DR>0,])/nrow(dft)
r0<-nrow(dft[dft$DR==0,])/nrow(dft)
rm<-nrow(dft[dft$DR<0,])/nrow(dft)

dbar<-data.frame(name=c("r>0", "r=0", "r<0"), value=c(rp, r0, rm))

barplot(height=dbar$value, names=dbar$name, space=1, ylim=c(0, 0.8))





which(rownames(Resp_p) %in% st)

#my_palette<-colorRampPalette(c("white", "black"))(n=20)
my_palette<-colorRampPalette(c("white", "darkblue"))(n=20)

#col_row=rainbow(length(unique(Expl_row$Class)))
#names(col_row)<-unique(unique(Expl_row$Class))

Class=c() #вектордля классов

Seasonal=c(spring_under_ice="blue", spring_open_water="deepskyblue", summer="firebrick1", autumn="gold")

Basin=c(Middle="deeppink4", South="forestgreen")

annoCol=list(Class=Class, Seasonal=Seasonal, Basin=Basin)

pheatmap(Resp_p, cluster_rows=tree_row, cluster_cols=tree_col, color=my_palette, fontsize=5, 
 annotation_row=Expl_row, annotation_col=Expl, annotation_colors=annoCol, cutree_rows=1, cutree_cols=4)

#///////////////////////////////////////////////////////
#///////////////////////////////////////////////////////




fit<-hclust(as.dist(DF), method = "average")
plot(fit, hang=-1)
k<-cutree(fit, h=0.2)
km<-table(k)
km<-km[km>1]
km<-as.numeric(names(km))
k<-k[which(k %in% km)]

kn<-unique(k)

DDR<-c()
DDF<-c()

for(i in 1: length(kn))
 {
  kv<-k[k==kn[i]]
  nr<-which(rownames(DD) %in% names(kv))
  DDV<-DD[nr,nr]
  DFV<-DF[nr,nr]
  DDV<-as.dist(DDV)
  DFV<-as.dist(DFV)
  DDR<-c(DDR, as.vector(DDV))
  DDF<-c(DDF, as.vector(DFV))   
 }

df<-data.frame(DDR, DDF)
