setwd("~/Downloads")

list.files()

theta <- read.table("AB.thetas.idx.pestPG")

colnames(theta) = c("window","chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

theta$tWsite = theta$tW/theta$numSites
theta$tPsite = theta$tP/theta$numSites
theta$diff = theta$tPsite-theta$tWsite
theta$Ne = theta$tPsite/(4*2E-9)

sfs<-scan('AB_allsites.sfs')
sfs<-sfs[-c(1,length(sfs))]
sfs<-sfs/sum(sfs)


par(mfrow=c(2,2))
hist(theta$tWsite)
hist(theta$tPsite)
barplot(sfs,names=1:length(sfs),main='SFS')
hist(theta$tajD)
dev.off()

summary(theta)
theta[which(theta$tajD>1.5 & theta$tPsite<0.001),]

PC <- read.table("EDGE_poly.covar")
PCA <- eigen(PC)$vectors

annot <- read.table("EDGE.annot")

plot(PCA[,1],PCA[,2], xlim=c(-.1,.1),ylim=c(-.1,.1),col="lightgray")
text(PCA[,1], PCA[,2],  annot$V2,
     cex=0.65, pos=3,col="red") 

plot(PCA[,1],PCA[,2],col="lightgray")
text(PCA[,1], PCA[,2],  annot$V1,
     cex=0.65, pos=3,col="red") 

write.table(PCA[,c(1)],"EDGE_PC1.txt",quote = F,row.names = F,col.names = F)
write.table(PCA[,c(2)],"EDGE_PC2.txt",quote = F,row.names = F,col.names = F)
