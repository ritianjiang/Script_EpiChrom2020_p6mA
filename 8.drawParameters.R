#This section draw the 4 difference Accuracy
setwd("/home/owht/KIZ/data/ML_6ma_2/Result/Train_tuning")
#library(scatterplot3d)
library(rgl)

#boxplot(bp41$Acc,bp151$Acc,bp301$Acc,bp601$Acc,ylim=c(0.8,1))

#Then we select the bp301 can see the result of parameter tuning
###The top 60% performs best , sowe plot the 3d plot of top 60%

load("Train_tuning.RData");total<-partable;rm(partable)
prc <- 0.7 * total$Acc / diff(range(total$Acc))
prc<-prc-min(prc) + 0.3
hsv(prc-0.01)->prc
plot3d(total[,1:3], mar = c(5, 3, 4, 3),
       col = prc,size = 12,type="p")
par(mar=c(5, 3, 0, 3))
plot(seq(min(total$Acc), max(total$Acc), length = 100), rep(0, 100), pch = 15,
     axes = FALSE, xlab = "color code of variable \"Precision\"",
     ylab = "", col = hsv(seq(0.3, 1, length = 100)))
rgl.postscript("plot3d_tuning_total.pdf", "pdf", drawText = TRUE) 

##ETA=0.3 gamma=0.16; max_depth=4
