library(tidyverse)
library(rhdf5)
library(stringr)
library(edgeR)
library(plotly)
library(TCGAbiolinks)
library(TCGAbiolinksGUI.data)

#import CYT data 
melanoma_CYT <- read_tsv("C:/Users/Xinyi/Desktop/omics/melanoma patient tcga/melanoma cyt.txt")
#tidy up patient cyt data
melanoma_dt <- select(melanoma_CYT, "PatientID", "Cytolytic Activity")

write_tsv(melanoma_dt, "melanoma_CYT_tidy.txt")

#tcga data
tcga_dt <- read_tsv("C:/Users/Xinyi/Desktop/omics/HiSeqV2")
head(tcga_dt) 

colnames(tcga_dt)<-str_sub(colnames(tcga_dt),1,12)

c_id <- melanoma_dt$PatientID


tcga_patient_dt <- tcga_dt[, (which(colnames(tcga_dt) %in% c_id))] %>%
  mutate("Gene"=tcga_dt$sample) %>%
  select("Gene", everything())

##anaplerosis gene expression in patients
anaplerosis_gene <- read_tsv("C:/Users/Xinyi/Desktop/omics/gene sets/anaplerosis.txt")
anaplerosis_gene <- anaplerosis_gene[-1,]
anaplerosis_gene <- anaplerosis_gene$KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM

patient_anaplerosis_expression <- tcga_patient_dt[(which(tcga_patient_dt$Gene %in% anaplerosis_gene)),]


##aatransport gene expression in patients
aatransport_gene <- read_tsv("C:/Users/Xinyi/Desktop/omics/gene sets/amino acid transport.txt")
aatransport_gene <- aatransport_gene[-1,]
aatransport_gene <- aatransport_gene$REACTOME_AMINO_ACID_TRANSPORT_ACROSS_THE_PLASMA_MEMBRANE

patient_aatransport_expression <- tcga_patient_dt[(which(tcga_patient_dt$Gene %in% aatransport_gene)),]

patient_aagene <- rbind(patient_aatransport_expression,patient_anaplerosis_expression)
patient_aagene_expression<-t(patient_aagene)
colnames(patient_aagene_expression)<-patient_aagene_expression[1,]
patient_aagene_expression<-patient_aagene_expression[-1,]

patient_aagene_expression_df<-as.data.frame(patient_aagene_expression)
patient_aagene_expression_df <-  mutate(patient_aagene_expression_df, "PatientID"=rownames(patient_aagene_expression_df))
dim(patient_aagene_expression_df)

data <- join(melanoma_dt,patient_aagene_expression_df, by="PatientID")


##process for heatmap
#patient_aagene_hm <- pivot_longer(patient_aagene, !Gene, names_to = "patient")
#write_tsv(patient_aagene_hm,"patient_aagene_expression.txt")


##combine CYT and expression into one df
results <- data
results <- data.frame(results)
rownames(results) <- results[,1]
results <- results[,-1]


## creat correlation matirx

library(Hmisc)
library(tidyverse)
library(rhdf5)
library(edgeR)

results1 <- results
#filter <- results1

results1[results1 == 0] <- NA
results.filter <- results1 %>% select_if(~!any(is.na(.)))


res <- cor(results.filter)
round(res, 2)
cor(results1, use = "complete.obs")

res2 <- rcorr(as.matrix(results.filter))
res2
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
flattenCorrMatrix(res2$r, res2$P)

library(corrplot)

# Insignificant correlations are leaved blank
corrplot(res2$r, type="upper", col = COL2('RdYlBu', 200),
         method = "circle",
         tl.col = 'black',tl.cex = 0.55,tl.srt = 45,
         p.mat = res2$P, sig.level = 0.05, insig = "blank")





#export correlation matrix value
P_value <- data.frame(res2$P)

r_value<- data.frame(res2$r)

export_correlation_value <- t(rbind(P_value[1,],r_value[1,]))
colnames(export_correlation_value) <- c("P_value","r_value")

correlation <- data.frame(export_correlation_value)
correlation_filter <- correlation[(which(correlation$P_value <0.05)),]
correlation_filter <- mutate(correlation_filter,"name"=rownames(correlation_filter))




write_tsv(data.frame(export_correlation_value), "cyt res value.txt")

ggplot(correlation_filter, aes(x=reorder(name, desc(r_value)),y=r_value,fill = -P_value))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, size=7.5,vjust = 1, hjust=1,face="bold"),
        axis.text.y = element_text(face="bold"))

#save image
png(file="C:/Users/Xinyi/Desktop/omics/melanoma patient tcga/melanoma patient aagene cyt barplot.png",
    width=600, height=400)
ggplot(correlation_filter, aes(x=reorder(name, desc(r_value)),y=r_value,fill = -P_value))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, size=7.5,vjust = 1, hjust=1,face="bold"),
        axis.text.y = element_text(face="bold"))
dev.off()
# save image --------------------------------------------------------------


png(file="C:/Users/Xinyi/Desktop/omics/melanoma patient tcga/melanoma patient aagene cyt correlation.png",
    width=1200, height=1200)
corrplot(res2$r, type="upper", col = COL2('RdYlBu', 200),
         method = "circle",
         tl.col = 'black',tl.cex = 0.65,tl.srt = 45,
         p.mat = res2$P, sig.level = 0.05, insig = "blank")
dev.off()