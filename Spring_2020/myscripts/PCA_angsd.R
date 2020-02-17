##################################
###      PCA Red Spruce       ####
##################################

setwd("~/Documents/ProjetRedSpruce/Analyses/Downstream_analyses/PCA/")

library(ade4)

## PCA from covariance matrix
TAB <- read.table("Full_Sampling_covmatrix.cov")
#TAB <- read.table("redspruce_WES_fullsampling_covmatrix.cov")
pca <- eigen(TAB)

## Info about individuals
names <- read.table("names_bam.list")
info_inds<-read.table("../../Info_samples_revised.txt", header=T)
info_inds<-info_inds[match(as.vector(names[,1]), as.character(info_inds$Family)),]

## Explain variance
var <- pca$values/sum(pca$values)
barplot(pca$values[1:100])
barplot(var[1:100])

## Figure
pdf("./PCA_RedSpruce.pdf")
plot(pca$vectors[,c(1,2)], type = "n", xlim = c(-0.1,0.1), xlab = "PC1 (14.8%)", ylab = "PC2 (0.5%)", main = "Complete sampling (3.5 millions SNP)")
par(xpd=FALSE)
abline(h = 0, lty=5, col = "black")
abline(v = 0, lty=5, col = "black")
s.class(pca$vectors[which(info_inds$Region=="M"),c(1,2)], fac=as.factor(info_inds$Site[which(info_inds$Region=="M")]), cellipse = 0, col = rep("#4DAF4A", length(levels(as.factor(info_inds$Site[which(info_inds$Region=="M")])))), add.plot = T, clabel = 0.8)
s.class(pca$vectors[which(info_inds$Region=="C"),c(1,2)], fac=as.factor(info_inds$Site[which(info_inds$Region=="C")]), cellipse = 0, col = rep("#FEDF00", length(levels(as.factor(info_inds$Site[which(info_inds$Region=="C")])))), add.plot = T, clabel = 0.8)
s.class(pca$vectors[which(info_inds$Region=="E"),c(1,2)], fac=as.factor(info_inds$Site[which(info_inds$Region=="E")]), cellipse = 0, col = rep("#377EB8", length(levels(as.factor(info_inds$Site[which(info_inds$Region=="E")])))), add.plot = T, clabel = 0.8)
legend(0.07,0.26,
       legend = c("Core", "Margin", "Edge"), border = "black", pch = 21, col = "black",
       pt.bg = c("#F9D017", "#4DAF4A", "#377EB8"), title = "Regions", text.font = 4, title.adj = 0.2)
dev.off()

