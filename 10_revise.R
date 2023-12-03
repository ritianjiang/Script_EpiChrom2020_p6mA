setwd("/home/owht/KIZ/data/ML_6ma_3/Result/Revised")
jk_res<-read.table("Jackknife_species",stringsAsFactors = F)
unique(jk_res$V3)
Fly_res<-jk_res[jk_res$V3 == "Fly",]
Human_res<-jk_res[jk_res$V3 == "Human",]
Rice_res<-jk_res[jk_res$V3 == "Rice",]
Worm_res<-jk_res[jk_res$V3 == "Worm",]
library(pROC)

roc<-roc(labs,resultTun,levels = c(1,0))
pred<-read.table("../Compare/MM-6mAPred")