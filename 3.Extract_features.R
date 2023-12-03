#This script will extract the features for the 4 datasets 
setwd("/home/owht/KIZ/data/ML_6ma_3/Dataset/")
library(Biostrings)
library(RTFE)
physio<-read.csv("/home/owht/myGIT/RTFE/data/feature_6_diprogb.csv",
                 row.names = 1)


##################################For total
rm(list = ls())
#setwd("41bp/")
###read fasta files
Pos<-readDNAStringSet(filepath = "total_Positive.fasta",
                      format = "fasta",use.names = F)
Pos<-Pos %>% as.character()
grep(x=Pos,pattern = "N")->filt
if(length(filt)==0){filt <- NULL}
Pos<-Pos[-filt]

Neg<-readDNAStringSet(filepath = "total_Negative.fasta",
                      format = "fasta",use.names = F)
Neg<-Neg %>% as.character()
grep(x=Neg,pattern = "N")->filt
if(length(filt)==0){filt <- NULL}
Neg<-Neg[-filt] #if error dont worry, it means nothing needs to be filtered
### make the Benneat obj
Bean<-makeBenneat(Seqs = c(Pos,Neg),labs = c(rep("Pos",length(Pos)),
                                             rep("Neg",length(Neg))))
All_Feature_41<-data.frame(seq=c(paste0("Pos",1:length(Pos)),paste0("Neg",1:length(Neg))))

##PSTNPss
PSTNPss<-PSTNPss(object = Bean,lableA = "Pos",lableB = "Neg")
colnames(PSTNPss)<-paste0("PSTNPss_",1:ncol(PSTNPss))
All_Feature_41<-cbind(All_Feature_41,PSTNPss)
rm(PSTNPss)
##PSTNPds
PSTNPds<-PSTNPds(object = Bean,lableA = "Pos",lableB = "Neg")
colnames(PSTNPds)<-paste0("PSTNPds_",1:ncol(PSTNPds))
All_Feature_41<-cbind(All_Feature_41,PSTNPds)
rm(PSTNPds)
##Type2-PseKNC lambda=5
physio<-read.csv("/home/owht/myGIT/RTFE/data/feature_6_diprogb.csv",
                 row.names = 1,stringsAsFactors = F)
PseKNC<-getT2PseKNC(object = Bean,phychem_file = physio,
                    normalization = T,lambda = 5)
n<-ncol(PseKNC)/nrow(physio)
paste0(paste0("T2P",lapply(1:5,rep,times=6) %>% unlist),
       rep(rownames(physio),5)) -> colnames(PseKNC)
All_Feature_41<-cbind(All_Feature_41,PseKNC)
rm(PseKNC)
##EIIP
EIIP<-getEIIP(object = Bean)
colnames(EIIP)<-paste0("EIIP_",colnames(EIIP))
All_Feature_41<-cbind(All_Feature_41,EIIP)
rm(EIIP)
###Save
colnames(All_Feature_41)
save(All_Feature_41,file="total_4_feature_41.RData")
########################


