#The Script should be run on Server

setwd("/data/user_data/wanght/ML_6ma_3")
library(magrittr)
library(xgboost)
###generate the parameter table
eta<-seq(0.1,0.5,0.1);depth<-seq(2,10,2);gamma<-seq(0,0.2,0.01)
partable<-expand.grid(eta,depth,gamma) %>% as.data.frame()
colnames(partable)<-c("eta","depth","gamma");partable$Acc<-0  



#For total 124 features
load("total_Feature_124.RData") #load the 41-bp features

top20<-Total_feature_124
labs<-rep(0,nrow(top20));
labs[grep(rownames(top20),pattern = "Pos")]<-1

dtrain<-xgb.DMatrix(top20%>% as.matrix(), label=labs)
for(i in 1:nrow(partable)){
  param<-list(eta=partable[i,]$eta,max_depth=partable[i,]$depth,
              gamma=partable[i,]$gamma)
  xgb.cv(param,dtrain,nrounds = 100,nfold = 10,
         metrics = {'error'},verbose=F)->res
  partable[i,]$Acc<-(1-res$evaluation_log$test_error_mean %>% tail(1))
  cat(i);cat("\n")
}
write.table(partable,file="total_124_feature.csv")

