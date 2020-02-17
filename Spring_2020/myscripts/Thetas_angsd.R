#####################################################################
###     Thetas - pi - Tajima's D distribution WES Red Spruce     ####
#####################################################################

library(ggplot2)
library(ggpubr)

setwd("~/Documents/ProjetRedSpruce/Analyses/Downstream_analyses/Thetas-Pi-TajimaD/")

##############
#### DATA ####

## Thetas
files <- list.files(path = "/Volumes/kellrlab/datashare/Spruce/exome_capture/WES_mapping/ANGSD/ref_Pabies/REGIONS/", pattern = "_intersect.thetas.idx.pestPG$", full.names = T, recursive = T)
thetas <- lapply(files, function(x) read.table(x, header = F)[,-1])
colnames(thetas[[1]]) <- c("Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
colnames(thetas[[2]]) <- c("Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")
colnames(thetas[[3]]) <- c("Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")

###########################
#### Distribution plot ####

table <- data.frame(do.call(rbind, thetas), variable = c(rep("CORE", nrow(thetas[[1]])), rep("EDGE", nrow(thetas[[2]])), rep("MARGIN", nrow(thetas[[3]])))) #, rep("FullSamlping", nrow(thetas[[3]]))
  
p1 <- ggplot() + 
  geom_density(data=table, aes(x=Tajima, fill=variable), lwd=0.5) +
  scale_fill_manual(values = c("#F9D017", "#377EB8", "#4DAF4A")) +
  facet_wrap(~variable, ncol = 1) +
  xlab("Tajima's D") +
  xlim(-3,3) +
  theme_bw(base_line_size = 0)

p2 <- ggplot() + 
  geom_density(data=table, aes(x=tW/nSites, fill=variable), lwd=0.5) +
  scale_fill_manual(values = c("#F9D017", "#377EB8", "#4DAF4A")) +
  xlab("ThetasW") +
  facet_wrap(~variable, ncol = 1) +
  xlim(0,0.03) +
  theme_bw(base_line_size = 0)

p3 <- ggplot() + 
  geom_density(data=table, aes(x=tP/nSites, fill=variable), lwd=0.5) +
  scale_fill_manual(values = c("#F9D017", "#377EB8", "#4DAF4A")) +
  xlab("ThetasP") +
  facet_wrap(~variable, ncol = 1) +
  xlim(0,0.05) +
  theme_bw(base_line_size = 0)

pdf("Tajima_Thetas_Pglauca.pdf", width = 7, height = 6)
ggarrange(p1,p2,p3, ncol=3, legend = "none", common.legend = T, widths = c(1,1))
dev.off()

###############
#### Menas ####

# Tajima
mean(thetas[[1]]$Tajima) # -0.06296943
mean(thetas[[2]]$Tajima) # 0.1479578
mean(thetas[[3]]$Tajima) # 0.1890214
# Waterson
mean(thetas[[1]]$tW/thetas[[1]]$nSites, na.rm = T) # 0.005236109
mean(thetas[[2]]$tW/thetas[[2]]$nSites, na.rm = T) # 0.004498055
mean(thetas[[3]]$tW/thetas[[3]]$nSites, na.rm = T) # 0.004646657
# Tajima
mean(thetas[[1]]$tP/thetas[[1]]$nSites, na.rm = T) # 0.00554475
mean(thetas[[2]]$tP/thetas[[2]]$nSites, na.rm = T) # 0.00522434
mean(thetas[[3]]$tP/thetas[[3]]$nSites, na.rm = T) # 0.00541076

