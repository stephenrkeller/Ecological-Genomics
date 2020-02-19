setwd("~/Downloads")

list.files()

dat <- read.table("AB.thetas.idx.pestPG")

colnames(dat) = c("window","chrname", "wincenter", "tW", "tP", "tF", "tH", "tL", "tajD", "fulif", "fuliD", "fayH", "zengsE", "numSites")

dat[which(dat$chrname=="chloroplast"),]

dat$tWsite = dat$tW/dat$numSites
dat$tPsite = dat$tP/dat$numSites
dat$diff = dat$tPsite-dat$tWsite
dat$Ne = dat$tPsite/(4*2E-9)

par(mfrow=c(2,2))
hist(dat$tWsite)
hist(dat$tPsite)
plot(dat$tPsite,dat$tWsite)
#hist(dat$diff)
#hist(dat$Ne)
hist(dat$tajD)

summary(dat)

