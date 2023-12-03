#We will train the models by the optimalized parameter sets
setwd("/home/owht/KIZ/data/ML_6ma_3/Result/Models")
library(magrittr)
library(RTFE)
library(xgboost)
library(Biostrings)

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



#For total, param = [eta=0.3 depth=4 gamma=0.16] Acc=0.8388 top 124
#rm(list = ls())
dir.create("total_models")

Pos41<-readDNAStringSet(filepath = "/home/owht/KIZ/data/ML_6ma_3/Dataset/total_Positive.fasta",
                        format = "fasta",use.names = F)
Pos41<-Pos41 %>% as.character()
Neg41<-readDNAStringSet(filepath = "/home/owht/KIZ/data/ML_6ma_3/Dataset/total_Negative.fasta",
                        format = "fasta",use.names = F)
Neg41<-Neg41 %>% as.character()

Bean41<-makeBenneat(Seqs = c(Pos41,Neg41),labs = c(rep("Pos",length(Pos41)),
                                                   rep("Neg",length(Neg41))))

seq_A<-Bean41@Seqs[Bean41@labs == "Pos"] #Get thesequence which is labeled as lableA
seq_B<-Bean41@Seqs[Bean41@labs == "Neg"] #Get thesequence which is labeled as lableA
Dic<-expand.grid(c("A","T","C","G"),c("A","T","C","G"),c("A","T","C","G"),
                 stringsAsFactors = F) %>%
  apply(1,paste,collapse="") %>% sort #Create the name vector for 3-mer
F_A<-matrix(-1,nrow=64,ncol=Bean41@seq_length-2)
F_B<-matrix(-1,nrow=64,ncol=Bean41@seq_length-2) #Calculate the frequency of each 3-mer in each position in both label-classes
for(i in 1:Bean41@seq_length-2){
  F_A[,i]<-lapply(seq_A, str_sub,start=i,end=i+2) %>%
    do.call(what = "rbind") %>%factor(levels = Dic) %>% table / length(seq_A)
}
for(i in 1:Bean41@seq_length-2){
  F_B[,i]<-lapply(seq_B, str_sub,start=i,end=i+2) %>%
    do.call(what = "rbind") %>%factor(levels = Dic) %>% table / length(seq_A)
}
Zss<- F_A - F_B;rownames(Zss)<-Dic
rm(seq_A,seq_B,F_A,F_B,Dic)
save(Zss,file="total_models/Zss_41.rds")

newSeq<-str_replace_all(string = Bean41@Seqs,pattern = "T",replacement = "A")
newSeq<-str_replace_all(string = newSeq,pattern = "G",replacement = "C")
seq_A<-newSeq[Bean41@labs == "Pos"] #Get thesequence which is labeled as lableA
seq_B<-newSeq[Bean41@labs == "Neg"] #Get thesequence which is labeled as lableA
Dic<-expand.grid(c("A","C"),c("A","C"),c("A","C"),
                 stringsAsFactors = F) %>%
  apply(1,paste,collapse="") %>% sort #Create the name vector for 3-mer
F_A<-matrix(-1,nrow=8,ncol=Bean41@seq_length-2)
F_B<-matrix(-1,nrow=8,ncol=Bean41@seq_length-2) #Calculate the frequency of each 3-mer in each position in both label-classes
for(i in 1:Bean41@seq_length-2){
  F_A[,i]<-lapply(seq_A, str_sub,start=i,end=i+2) %>%
    do.call(what = "rbind") %>%factor(levels = Dic) %>% table / length(seq_A)
}
for(i in 1:Bean41@seq_length-2){
  F_B[,i]<-lapply(seq_B, str_sub,start=i,end=i+2) %>%
    do.call(what = "rbind") %>%factor(levels = Dic) %>% table / length(seq_A)
}
Zds<- F_A - F_B;rownames(Zds)<-Dic
rm(seq_A,seq_B,F_A,F_B,Dic,newSeq,Pos41,Neg41,Bean41)
save(Zds,file="total_models/Zds_41.rds")

load("/home/owht/KIZ/data/ML_6ma_3/Result/features/total_Feature41.RData")

top20<-All_Feature_41[,mrmd[1:124,]$feature]
feaN<-colnames(top20)
save(feaN,file="total_models/Features.RData")
labs<-rep(0,nrow(top20));
labs[grep(rownames(top20),pattern = "Pos")]<-1

dtrain<-xgb.DMatrix(top20%>% as.matrix(), label=labs)
dtrain
param<-list(eta=0.3,depth=4,
            gamma=0.16)
xgb.train(param,dtrain,nrounds = 100,
          verbose=F)->total_train_model
save(total_train_model,file="total_models/total_train_model.rds")
xgb.save(total_train_model,fname = "total_models/total_train_model")

#####################################################
#jackknife 
library(pROC)
resultTun<-rep(0,6080)
labs<-c(rep(1,3040),rep(0,3040))
for(i in 1:6080){
  dtrain<-xgb.DMatrix(top20[-i,] %>% as.matrix(), label=labs[-i])
  model<-xgb.train(data = dtrain,nrounds = 100,
                   verbose=F)
  dtest<-xgb.DMatrix(top20[i,] %>% as.matrix())
  resultTun[i]<-predict(model,dtest)
}

roc<-roc(labs,resultTun,levels = c(1,0))
pred<-read.table("../Compare/MM-6mAPred")

