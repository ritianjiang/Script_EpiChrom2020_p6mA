setwd("/home/owht/KIZ/data/ML_6ma_3/Result/")
dir.create("features")
#library(RTFE)
library(stringr)

##total
rm(list = ls())
load("../Dataset/AllFeatures/total_4_feature_41.RData")
rownames(All_Feature_41)<-All_Feature_41$seq
All_Feature_41<-All_Feature_41[,-1]
grep(rownames(All_Feature_41),pattern = "Pos")->Pos

#Feature41F<-Fscore(features = All_Feature_41,Pos)
#save(All_Feature_41,file = "features/total_Feature41.RData")

#MRMD arff format
sink("total_feature_4MRMD.arff.header")
cat("@relation Total");cat("\n");cat("\n")
for(i in 1:172){
  line=paste0("@attribute ",colnames(All_Feature_41)[i],
              " real")
  cat(line);cat("\n")
}
cat("\n");cat("@attribute class {1,-1}");cat("\n")
cat("@data")
sink()

All_Feature_41$class<--1
All_Feature_41$class[1:3040]<-1
write.csv(All_Feature_41,quote=F,row.names = F,
          file="total_feature_4MRMD.arff.data")

#### MRMD #####################3
# java -jar ~/KIZ/software/mrmd/mrmd.jar 
#  -i total_feature_4MRMD.arff -o total_feature_4MRMD.mrmd 
#######################3

mrmd<-read.table("total_feature_4MRMD.mrmd",header = F,
                 stringsAsFactors = F)
mrmd$feature<-colnames(All_Feature_41)[1+as.numeric(str_split(mrmd$V2,pattern = "Fea",simplify = T)[,2])]
  
save(All_Feature_41,
     mrmd,file = "features/total_Feature41.RData")
