#This script draw the Entropy figure

setwd("/home/owht/KIZ/data/ML_6ma_2/Result/Entropy")
library(Biostrings)
library(magrittr)
library(stringr)

#Entropy Function
entropy<-function(x){
  n<-length(x)
  se<-(table(x) %>% as.numeric)/n
  return(-sum(se*log2(se)))
}

par(mfrow=c(2,3))

##Rice
Pos41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/Rice_Positive.fasta",
                      format = "fasta",use.names = F)
Pos41<-as.character(Pos41)
lapply(Pos41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Pos41),nrow=41) %>%
  t -> Pos41M
#plot(1:41,apply(Pos41M, 2, entropy),"l",col="red",ylim=c(0,3))
Neg41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/Rice_Negative.fasta",
                      format = "fasta",use.names = F)
Neg41<-as.character(Neg41)
lapply(Neg41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Neg41),nrow=41) %>%
  t -> Neg41M
#lines(1:41,apply(Neg41M, 2, entropy),"l",col="blue")

Pos41Mer<-matrix("N",nrow=length(Pos41),ncol=39)
for(i in 1:39){
  Pos41Mer[,i]<-paste0(Pos41M[,i],Pos41M[,i+1],Pos41M[,i+2])
}
plot(1:39,apply(Pos41Mer, 2, entropy),"o",col="red",pch=16,ylim=c(3,6),
     ylab="Entropy")
Neg41Mer<-matrix("N",nrow=length(Neg41),ncol=39)
for(i in 1:39){
  Neg41Mer[,i]<-paste0(Neg41M[,i],Neg41M[,i+1],Neg41M[,i+2])
}
lines(1:39,apply(Neg41Mer, 2, entropy),"o",pch=16,col="blue")
text(x=5,y=3.5,"Rice")

##Fly
Pos41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/Fly_Positive_dm3_nr0.6.fasta",
                      format = "fasta",use.names = F)
Pos41<-as.character(Pos41)
lapply(Pos41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Pos41),nrow=41) %>%
  t -> Pos41M
#plot(1:41,apply(Pos41M, 2, entropy),"l",col="red",ylim=c(0,3))
Neg41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/Fly_Negative_dm3_nr0.6.fasta",
                      format = "fasta",use.names = F)
Neg41<-as.character(Neg41)
lapply(Neg41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Neg41),nrow=41) %>%
  t -> Neg41M
#lines(1:41,apply(Neg41M, 2, entropy),"l",col="blue")

Pos41Mer<-matrix("N",nrow=length(Pos41),ncol=39)
for(i in 1:39){
  Pos41Mer[,i]<-paste0(Pos41M[,i],Pos41M[,i+1],Pos41M[,i+2])
}
plot(1:39,apply(Pos41Mer, 2, entropy),"o",col="red",pch=16,ylim=c(3,6),
     ylab="Entropy")
Neg41Mer<-matrix("N",nrow=length(Neg41),ncol=39)
for(i in 1:39){
  Neg41Mer[,i]<-paste0(Neg41M[,i],Neg41M[,i+1],Neg41M[,i+2])
}
lines(1:39,apply(Neg41Mer, 2, entropy),"o",pch=16,col="blue")
text(x=5,y=3.5,"Fly")

##Worm
Pos41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/Worm_Positive_ce10_nr0.6.fasta",
                      format = "fasta",use.names = F)
Pos41<-as.character(Pos41)
lapply(Pos41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Pos41),nrow=41) %>%
  t -> Pos41M
#plot(1:41,apply(Pos41M, 2, entropy),"l",col="red",ylim=c(0,3))
Neg41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/Worm_Negative_ce10_nr0.6.fasta",
                      format = "fasta",use.names = F)
Neg41<-as.character(Neg41)
lapply(Neg41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Neg41),nrow=41) %>%
  t -> Neg41M
#lines(1:41,apply(Neg41M, 2, entropy),"l",col="blue")

Pos41Mer<-matrix("N",nrow=length(Pos41),ncol=39)
for(i in 1:39){
  Pos41Mer[,i]<-paste0(Pos41M[,i],Pos41M[,i+1],Pos41M[,i+2])
}
plot(1:39,apply(Pos41Mer, 2, entropy),"o",col="red",pch=16,ylim=c(3,6),
     ylab="Entropy")
Neg41Mer<-matrix("N",nrow=length(Neg41),ncol=39)
for(i in 1:39){
  Neg41Mer[,i]<-paste0(Neg41M[,i],Neg41M[,i+1],Neg41M[,i+2])
}
lines(1:39,apply(Neg41Mer, 2, entropy),"o",pch=16,col="blue")
text(x=5,y=3.5,"Worm")

##Human
Pos41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/Human_Positive_nr0.6.fasta",
                      format = "fasta",use.names = F)
Pos41<-as.character(Pos41)
lapply(Pos41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Pos41),nrow=41) %>%
  t -> Pos41M
#plot(1:41,apply(Pos41M, 2, entropy),"l",col="red",ylim=c(0,3))
Neg41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/Human_Negative_nr0.6.fasta",
                      format = "fasta",use.names = F)
Neg41<-as.character(Neg41)
lapply(Neg41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Neg41),nrow=41) %>%
  t -> Neg41M
#lines(1:41,apply(Neg41M, 2, entropy),"l",col="blue")

Pos41Mer<-matrix("N",nrow=length(Pos41),ncol=39)
for(i in 1:39){
  Pos41Mer[,i]<-paste0(Pos41M[,i],Pos41M[,i+1],Pos41M[,i+2])
}
plot(1:39,apply(Pos41Mer, 2, entropy),"o",col="red",pch=16,ylim=c(3,6),
     ylab="Entropy")
Neg41Mer<-matrix("N",nrow=length(Neg41),ncol=39)
for(i in 1:39){
  Neg41Mer[,i]<-paste0(Neg41M[,i],Neg41M[,i+1],Neg41M[,i+2])
}
lines(1:39,apply(Neg41Mer, 2, entropy),"o",pch=16,col="blue")
text(x=5,y=3.5,"Human")

##Compre1
Pos41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/total_Positive.fasta",
                      format = "fasta",use.names = F)
Pos41<-as.character(Pos41)
lapply(Pos41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Pos41),nrow=41) %>%
  t -> Pos41M
#plot(1:41,apply(Pos41M, 2, entropy),"l",col="red",ylim=c(0,3))
Neg41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/total_Negative.fasta",
                      format = "fasta",use.names = F)
Neg41<-as.character(Neg41)
lapply(Neg41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Neg41),nrow=41) %>%
  t -> Neg41M
#lines(1:41,apply(Neg41M, 2, entropy),"l",col="blue")

Pos41Mer<-matrix("N",nrow=length(Pos41),ncol=39)
for(i in 1:39){
  Pos41Mer[,i]<-paste0(Pos41M[,i],Pos41M[,i+1],Pos41M[,i+2])
}
plot(1:39,apply(Pos41Mer, 2, entropy),"o",col="red",pch=16,ylim=c(3,6),
     ylab="Entropy")
Neg41Mer<-matrix("N",nrow=length(Neg41),ncol=39)
for(i in 1:39){
  Neg41Mer[,i]<-paste0(Neg41M[,i],Neg41M[,i+1],Neg41M[,i+2])
}
lines(1:39,apply(Neg41Mer, 2, entropy),"o",pch=16,col="blue")
text(x=5,y=3.5,"Comprehensive")

##Compre2
Pos41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/total_Positive_nr0.6.fasta",
                      format = "fasta",use.names = F)
Pos41<-as.character(Pos41)
lapply(Pos41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Pos41),nrow=41) %>%
  t -> Pos41M
#plot(1:41,apply(Pos41M, 2, entropy),"l",col="red",ylim=c(0,3))
Neg41<-readBStringSet("/home/owht/KIZ/data/ML_6ma_2/Dataset/total_Negative_nr0.6.fasta",
                      format = "fasta",use.names = F)
Neg41<-as.character(Neg41)
lapply(Neg41,str_split,pattern="")%>%unlist %>% 
  matrix(ncol=length(Neg41),nrow=41) %>%
  t -> Neg41M
#lines(1:41,apply(Neg41M, 2, entropy),"l",col="blue")

Pos41Mer<-matrix("N",nrow=length(Pos41),ncol=39)
for(i in 1:39){
  Pos41Mer[,i]<-paste0(Pos41M[,i],Pos41M[,i+1],Pos41M[,i+2])
}
plot(1:39,apply(Pos41Mer, 2, entropy),"o",col="red",pch=16,ylim=c(3,6),
     ylab="Entropy")
Neg41Mer<-matrix("N",nrow=length(Neg41),ncol=39)
for(i in 1:39){
  Neg41Mer[,i]<-paste0(Neg41M[,i],Neg41M[,i+1],Neg41M[,i+2])
}
lines(1:39,apply(Neg41Mer, 2, entropy),"o",pch=16,col="blue")
text(x=5,y=3.5,"Compre NR")
