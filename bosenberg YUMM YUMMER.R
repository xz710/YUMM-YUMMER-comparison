setwd("C:/Users/Xinyi/Desktop/omics/YUMM YUMMER comparison")


library(tidyverse)
library(rhdf5)
library(dplyr)
library(plyr)
library(stringr)
library(edgeR)
library(plotly)

cell_dt <- read_tsv("htseq-count_normalized.txt")
colnames(cell_dt) <- c("GeneID","YUMM1","YUMMER1","YUMMER2", "YUMM2","YUMM3","YUMMER3","YUMM3",
                      "YUMMER4","YUMMER5","RAG_YUMMER1","RAG_YUMMER2","YUMMER6",
                       "R_R","Y_Y")
cell_dt_norm <- data.frame(cell_dt)
rownames(cell_dt_norm)<- cell_dt$GeneID
cell_dt_norm<-cell_dt_norm[,-1]
cell_dt_norm<-select(cell_dt_norm,"YUMM1","YUMM2","YUMM3","YUMM4","YUMMER1","YUMMER2","YUMMER3",
                     "YUMMER4","YUMMER5","RAG_YUMMER1","RAG_YUMMER2")


keepers <- rowSums(cell_dt_norm>1)>=4
cell_norm_filter <- cell_dt_norm[keepers,]

cell.filtered.norm.log2.cpm <- cpm(cell_norm_filter, log=TRUE)


StudyDesign <- tibble(Sample_title = colnames(cell_dt_norm), 
                      tumorType = c( "YUMM1.7", "YUMM1.7","YUMM1.7","YUMM1.7", "YUMMER1.7", "YUMMER1.7","YUMMER1.7",
                                     "YUMMER1.7", "YUMMER1.7","YUMMER1.7","YUMMER1.7"))

TumorType <- factor (StudyDesign$tumorType)

pca.res <- prcomp(t(cell.filtered.norm.log2.cpm), scale.=F, retx=T)
pca.res.df <- as_tibble(pca.res$x)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, colour=TumorType) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

png(file="C:/Users/Xinyi/Desktop/omics/YUMM YUMMER comparison/DEG pca plot.png",
    width=300, height=300)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, colour=TumorType) +
  geom_point(size=4) +
  # geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
dev.off()
