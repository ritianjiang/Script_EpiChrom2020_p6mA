#Add features one by one according to the rank and plot the Acc lines
setwd("/home/owht/KIZ/data/ML_6ma_3/Result")
library(xgboost)
library(magrittr)

#total
rm(list = ls())
load("features/total_Feature41.RData")
All_Feature_41<-All_Feature_41[,-173]
labs<-rep(0,nrow(All_Feature_41));
labs[grep(rownames(All_Feature_41),pattern = "Pos")]<-1
All_Feature_41<-All_Feature_41[,mrmd$feature]
AccVec<-rep(0,ncol(All_Feature_41))
for(i in 1:ncol(All_Feature_41)){
  dtrain<-xgb.DMatrix(All_Feature_41[,1:i] %>% as.matrix(), label=labs)
  xgb.cv(data = dtrain,nrounds = 100,nfold = 10,
         metrics = {'error'},verbose = T,early_stopping_rounds = 60)->res
  AccVec[i]<-(1-res$evaluation_log$test_error_mean %>% tail(1))
}
order(AccVec,decreasing = T); AccVec %>% max
pdf("MRMD/total_MRMD.pdf")
plot(1:172,AccVec,ylim=c(0.45,0.85),"l",main="total",xaxt="n")
axis(side = 1,at = seq(0,180,20))
text(x =order(AccVec,decreasing = T)[1],y=0.8,
     paste0(order(AccVec,decreasing = T)[1],",",AccVec %>% max))
dev.off()
save(AccVec,file="MRMD/total_acc_vec.RData")

############### AUC ##############
library(pROC)

resultAll<-rep(0,6080)
result124<-resultAll

#All
for(i in 1:6080){
  dtrain<-xgb.DMatrix(All_Feature_41[-i,] %>% as.matrix(), label=labs[-i])
  model<-xgb.train(data = dtrain,nrounds = 100,
                   verbose=F)
  dtest<-xgb.DMatrix(All_Feature_41[i,] %>% as.matrix())
  resultAll[i]<-predict(model,dtest)
}
#124
for(i in 1:6080){
  dtrain<-xgb.DMatrix(All_Feature_41[-i,1:124] %>% as.matrix(), label=labs[-i])
  model<-xgb.train(data = dtrain,nrounds = 100,
                   verbose=F)
  dtest<-xgb.DMatrix(All_Feature_41[i,1:124] %>% as.matrix())
  result124[i]<-predict(model,dtest)
}

library(pROC)
#rocAll<-roc(labs,resultAll,levels=c(1,0)) #0.8896
#roc124<-roc(labs,result124,levels=c(1,0)) #0.8918
plot.roc(labs,resultAll)
lines.roc(labs,result124,col="#008600")
legend("bottomright",legend=c("All features","Selected features"),
       col=c("black","#008600"),lwd=2)
save(rocAll,roc124,file="./MRMD/jackknife_roc.RData")
Total_feature_124<-All_Feature_41[,1:124]
save(Total_feature_124,file="./Tuning/total_Feature_124.RData")
