#2_1.Comprehensive samples
#Now we use 'cat' to obtain the comprehensive dataset, named "total_Negative/Positive.fasta"

setwd("/home/owht/KIZ/data/ML_6ma_2/Dataset")
Pos<-readDNAStringSet("total_Positive.fasta",format = "fasta")
#Pos<-as.character(Pos)
grep(Pos,pattern = "N")

newOrder<-sample(x = 1:3040,size = 3040,replace = F)
Pos<-Pos[newOrder]
writeXStringSet(Pos,"total_Positive.fasta",format = "fasta")


Neg<-readDNAStringSet("total_Negative.fasta",format = "fasta")
#Pos<-as.character(Pos)
grep(Neg,pattern = "N")

newOrder<-sample(x = 1:3040,size = 3040,replace = F)
Neg<-Neg[newOrder]
writeXStringSet(Neg,"total_Negative.fasta",format = "fasta")
