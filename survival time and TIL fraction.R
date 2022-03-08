setwd("C:/Users/Xinyi/Desktop/omics/melanoma patient tcga")

library(tidyverse)
library(rhdf5)
library(edgeR)
library(ggplot2)

patient_survival_dt <- read_tsv("survival time and lymphocyte fraction.txt")
colnames(patient_survival_dt)<- c("id","day","leukocyte","lym_fraction","lymphocyte")


leukocyte_plot<- ggplot(patient_survival_dt,aes(leukocyte,day),shape = circle)+
  geom_point()+
  theme_bw()

lymphocyte_plot<- ggplot(patient_survival_dt,aes(lymphocyte,day),shape = circle)+
  geom_point()+
  theme_bw()

lmOS_time <- lm(day~leukocyte+lymphocyte, data = patient_survival_dt) #Create the linear regression
summary(lmOS_time) #Review the results

#CYT and survival time
melanoma_CYT <- read_tsv("C:/Users/Xinyi/Desktop/omics/melanoma patient tcga/melanoma cyt.txt")
#tidy up patient cyt data
melanoma_dt <- select(melanoma_CYT, "PatientID", "Cytolytic Activity")
colnames(melanoma_dt) <-c("id","cyt")

patient_cyt_survival <- join(melanoma_dt, patient_survival_dt, by="id")
patient_cyt_survival.filter <- na.omit(patient_cyt_survival) 

ggplot(patient_cyt_survival.filter,aes(cyt,day),shape = circle)+
  geom_point()+
  theme_bw()
lmcyt_time <- lm(day~cyt, data = patient_cyt_survival) #Create the linear regression
summary(lmcyt_time)


#saving images
png(file="C:/Users/Xinyi/Desktop/omics/melanoma patient tcga/melanoma patient leukocyte survival time.png",
    width=500, height=400)
ggplot(patient_survival_dt,aes(leukocyte,day),shape = circle)+
  geom_point()+
  theme_bw()
dev.off()

png(file="C:/Users/Xinyi/Desktop/omics/melanoma patient tcga/melanoma patient lymphocyte survival time.png",
    width=500, height=400)
ggplot(patient_survival_dt,aes(lymphocyte,day),shape = circle)+
  geom_point()+
  theme_bw()
dev.off()
  