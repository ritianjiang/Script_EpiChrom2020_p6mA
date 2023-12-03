##
## USAGE: Rscript p6mA.R [MODE] [INPUT FASTA FILE] [OUTPUT FILE NAME] [CUTOFF]
## 
##        INPUT FASTA FILE: The input .fasta file. The sequences should be 41-bp and the A in 21st base
##        OUTPUT FILE NAME: The output file name. Default: Output.result
##
##        Attention: If you want to modify CUTOFF, please set the OUTPUT FILE NAME
####################################################################

#Obtain the parameters.
Args<-commandArgs(trailingOnly=T)
if(length(Args) == 0){
  cat("USAGE: Rscript p6mA.R [INPUT FASTA FILE] [OUTPUT FILE NAME]");cat("\n");cat("\n")
  #cat("MODE            : The model to test the input sequences; one of Rice, Fly, Worm, Human, Compre.n");cat("\n")
  cat("INPUT FASTA FILE: The input .fasta file. The sequences should be 41-bp and the A in 21st base");cat("\n")
  cat("OUTPUT FILE NAME: The output file name. Default: output.result");cat("\n")
  #cat("CUTOFF          : The cutoff of classifier. Sites whose prediction socre below it will be identified as non-6mA. Default: 0.5");cat("\n")
  cat("\n");
  #cat("Attention: If you want to modify CUTOFF, please set the OUTPUT FILE NAME")
  stop("Please read the usage.",call.=F)
}


#Load the packages
library(RTFE,verbose = F)
library(magrittr,verbose = F)
library(stringr,verbose = F)
library(xgboost,verbose = F)
library(Biostrings,verbose = F)

#Utilities defination
PssZ<-function(object,Z){
  P<-matrix(-1,nrow=object@seq_num,ncol=ncol(Z)) #get the Features
  for(i in 1:ncol(Z)){
    Z[lapply(object@Seqs, str_sub,start=i,end=i+2) %>% do.call(what="rbind"),i]->P[,i]
  }
  return(P %>% as.data.frame())
}

PdsZ<-function(object,Z){
  newSeq<-str_replace_all(string = object@Seqs,pattern = "T",replacement = "A")
  newSeq<-str_replace_all(string = newSeq,pattern = "G",replacement = "C")
  P<-matrix(-1,nrow=object@seq_num,ncol=ncol(Z)) #get the Features
  for(i in 1:ncol(Z)){
    Z[lapply(newSeq, str_sub,start=i,end=i+2) %>% do.call(what="rbind"),i]->P[,i]
  }
  return(P %>% as.data.frame())
}


cat("Runing");cat("\n")

output<-"output.result"
#mode<-Args[1] %>% as.character; input<-Args[2] %>% as.character
try(output<-Args[1] %>% as.character(),silent = T); 
if(!is.na(output)){output<-Args[2]}
#try(cutoff<-Args[4],silent = T)
if(is.na(cutoff)){cutoff<-0.5}

cat(paste0("The model is ",mode,". The input file is ",input));cat("\n")
cat(paste0("The cutoff is ",cutoff));cat("\n")


xgb.load("Models/total_models/total_train_model")->modelA
load("Models/total_models/Zds_41.rds")
load("Models/total_models/Zss_41.rds")
load("Models/total_models/Features.RData")


#extract features
Seq<-readDNAStringSet(input,format = "fasta",use.names = T)
result<-data.frame(Seqnames=names(Seq),
                   Score=0,type = "N",stringsAsFactors = F)
Seq<-as.character(Seq)
Bean<-makeBenneat(Seqs = Seq,labs = rep("test",length(Seq)))
All_Feature<-data.frame(seq=names(Seq))

physio<-read.csv("Models/feature_6_diprogb.csv",
                 row.names = 1,stringsAsFactors = F)
PseKNC<-getT2PseKNC(object = Bean,phychem_file = physio,
                    normalization = T,lambda = 5)
n<-ncol(PseKNC)/nrow(physio)
paste0(paste0("T2P",lapply(1:5,rep,times=6) %>% unlist),
       rep(rownames(physio),5)) -> colnames(PseKNC)
All_Feature<-cbind(All_Feature,PseKNC)
rm(PseKNC)
##EIIP
EIIP<-getEIIP(object = Bean)
colnames(EIIP)<-paste0("EIIP_",colnames(EIIP))
All_Feature<-cbind(All_Feature,EIIP)
rm(EIIP)
##PSTNPss
PSTNPss<-PssZ(object = Bean,Z = Zss)
colnames(PSTNPss)<-paste0("PSTNPss_",1:ncol(PSTNPss))
All_Feature<-cbind(All_Feature,PSTNPss)
rm(PSTNPss)
##PSTNPds
PSTNPds<-PdsZ(object = Bean,Z = Zds)
colnames(PSTNPds)<-paste0("PSTNPds_",1:ncol(PSTNPds))
All_Feature<-cbind(All_Feature,PSTNPds)
rm(PSTNPds)

All_Feature<-All_Feature[,feaN]
dtest<-xgb.DMatrix(All_Feature%>% as.matrix())

result$Score<-predict(modelA,dtest)
result[result$Score>0.5,]$type<-"6mA";
try(result[result$Score<0.5,]$type<-"non-6mA",silent = T)
cat("\n")
cat("The procedure has ended.");cat("\n")
cat(paste0("There are ",length(Seq)," sequences. And ",sum(result$Score>0.5),
           " Argines are identified as 6mA."));cat("\n")
cat(paste0("The Positive Ratio is ",sum(result$Score>0.5)/length(Seq)))
cat("\n")

write.table(result,file=output)

