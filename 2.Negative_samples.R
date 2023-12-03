setwd("/home/owht/KIZ/data/ML_6ma_2/Dataset")
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(BSgenome.Celegans.UCSC.ce10)
library(magrittr)
library(stringr)

####################C.elegans ce10####################################
###Generate the A pool
APool<-data.frame(seqnames = rep("chr",100000),start=0,end=0,strand="*",
                  context = "N",stringsAsFactors = F)
BSgenome.Celegans.UCSC.ce10@seqinfo %>% as.data.frame(stringsAsFactors=F)->SeqInfo

#Generate the random site. It's time consuming and I will change the step. It takes for about 3 days
for(i in 63323:100000){
  seqn<-sample(rownames(SeqInfo)[1:6],1)
  loci<-sample(1:SeqInfo[seqn,]$seqlengths,1)
  strr<-sample(c("+","-"),1)
  temp<-getSeq(BSgenome.Celegans.UCSC.ce10,
               names=seqn,
               start=loci,
               end=loci,strand=strr,as.character = T)
  APool[i,]$seqnames<-seqn;APool[i,]$start = loci;APool[i,]$end = loci
  APool[i,]$context<-temp;APool[i,]$strand<-strr
}
rm(a,i,loci,seqn,strr,temp)

#Extract the A random sites
APool<-APool[APool$context == "A",]

#get 1000 surrounding regions of SMRT
SMRT<-read.table("GFF/celegans_modifications_ce10_new_simple.gff",
                 stringsAsFactors = F,header = T);
SMRT<-SMRT[,c(1,4,5)]
#SMRT$V1<-str_replace(SMRT$V1,pattern = "Chr0",replacement = "chr")
#SMRT$V1<-str_replace(SMRT$V1,pattern = "Chr",replacement = "chr")
colnames(SMRT)<-c("seqnames","start","end")
SMRT$strand<-"*"
SMRT$start<-SMRT$start - 500; SMRT$end<-SMRT$end + 500
colnames(SMRT)[1]<-"seqnames"

#filter Neg ranges
toFilt<-findOverlaps(as(APool,"GRanges"),as(SMRT,"GRanges")) %>% as.data.frame()
BPool<-APool[-toFilt$queryHits,]

temp<-BPool
temp$start<-temp$start - 20; temp$end<-temp$end + 20
loc_41<-as(temp,"GRanges")
loc_41<-loc_41[loc_41@strand != "*" & loc_41$context == "A",]
context<-getSeq(BSgenome.Celegans.UCSC.ce10,loc_41,
                as.character=T)
grep(context,pattern = "N") #There is no sequences which contains "N"

sink(file="/home/owht/KIZ/data/ML_6ma_2/Dataset/Worm_Negative_ce10.fasta")
for(i in 1:length(loc_41)){
  titl<-paste0(">Worm_Neg_",i)
  cat(titl)
  cat("\n")
  cat(context[i])
  cat("\n")
}
sink()

####################Fly dm3####################################
###Generate the A pool
rm(list = ls())
APool<-data.frame(seqnames = rep("chr",140000),start=0,end=0,strand="*",
                  context = "N",stringsAsFactors = F)
BSgenome.Dmelanogaster.UCSC.dm3@seqinfo %>% as.data.frame(stringsAsFactors=F)->SeqInfo
SMRT<-read.table("GFF/dro_modifications_dm3_6mA_simple.gff",
                 stringsAsFactors = F,header = F);
SeqInfo1<-intersect(rownames(SeqInfo),unique(SMRT$V1) %>% unique)

#Generate the random site. It's time consuming and I will change the step. It takes for about 3 days
for(i in 105685:140000){
  seqn<-sample(SeqInfo1,1)
  loci<-sample(1:SeqInfo[seqn,]$seqlengths,1)
  strr<-sample(c("+","-"),1)
  temp<-getSeq(BSgenome.Dmelanogaster.UCSC.dm3,
               names=seqn,
               start=loci,
               end=loci,strand=strr,as.character = T)
  APool[i,]$seqnames<-seqn;APool[i,]$start = loci;APool[i,]$end = loci
  APool[i,]$context<-temp;APool[i,]$strand<-strr
}
rm(a,i,loci,seqn,strr,temp)

#Extract the A random sites
APool<-APool[APool$context == "A",]

#get 1000 surrounding regions of SMRT
SMRT<-SMRT[,c(1,4,5)]
#SMRT$V1<-str_replace(SMRT$V1,pattern = "Chr0",replacement = "chr")
#SMRT$V1<-str_replace(SMRT$V1,pattern = "Chr",replacement = "chr")
colnames(SMRT)<-c("seqnames","start","end")
SMRT$strand<-"*"
SMRT$start<-SMRT$start - 500; SMRT$end<-SMRT$end + 500
colnames(SMRT)[1]<-"seqnames"

#filter Neg ranges
toFilt<-findOverlaps(as(APool,"GRanges"),as(SMRT,"GRanges")) %>% as.data.frame()
BPool<-APool[-toFilt$queryHits,]

temp<-BPool
temp$start<-temp$start - 20; temp$end<-temp$end + 20
loc_41<-as(temp,"GRanges")
loc_41<-loc_41[loc_41@strand != "*" & loc_41$context == "A",]
context<-getSeq(BSgenome.Dmelanogaster.UCSC.dm3,loc_41,
                as.character=T)
grep(context,pattern = "N") #There is no sequences which contains "N"
context<-context[-grep(context,pattern = "N")] %>% na.omit()

sink(file="/home/owht/KIZ/data/ML_6ma_2/Dataset/Fly_Negative_dm3.fasta")
for(i in 1:length(loc_41)){
  titl<-paste0(">Fly_Neg_",i)
  cat(titl)
  cat("\n")
  cat(context[i])
  cat("\n")
}
sink()

####################human hg38####################################
###Generate the A pool
rm(list = ls())
APool<-data.frame(seqnames = rep("chr",180000),start=0,end=0,strand="*",
                  context = "N",stringsAsFactors = F)
BSgenome.Hsapiens.UCSC.hg38@seqinfo %>% as.data.frame(stringsAsFactors=F)->SeqInfo
SMRT<-read.csv("/home/owht/KIZ/data/ML_6ma/Dataset/MolCell2018/All_6mA_bySMRT.csv",
                 stringsAsFactors = F,header = T);
SeqInfo1<-intersect(rownames(SeqInfo),unique(SMRT$Chromsome) %>% unique)

#Generate the random site. It's time consuming and I will change the step. It takes for about 3 days
for(i in 111109:180000){
  seqn<-sample(SeqInfo1,1)
  loci<-sample(1:SeqInfo[seqn,]$seqlengths,1)
  strr<-sample(c("+","-"),1)
  temp<-getSeq(BSgenome.Hsapiens.UCSC.hg38,
               names=seqn,
               start=loci,
               end=loci,strand=strr,as.character = T)
  APool[i,]$seqnames<-seqn;APool[i,]$start = loci;APool[i,]$end = loci
  APool[i,]$context<-temp;APool[i,]$strand<-strr
}
rm(a,i,loci,seqn,strr,temp)

#Extract the A random sites
APool<-APool[APool$context == "A",]

#get 1000 surrounding regions of SMRT
SMRT<-SMRT[,c(1,2,2)]
#SMRT$V1<-str_replace(SMRT$V1,pattern = "Chr0",replacement = "chr")
#SMRT$V1<-str_replace(SMRT$V1,pattern = "Chr",replacement = "chr")
colnames(SMRT)<-c("seqnames","start","end")
SMRT$strand<-"*"
SMRT$start<-SMRT$start - 500; SMRT$end<-SMRT$end + 500
colnames(SMRT)[1]<-"seqnames"

#filter Neg ranges
toFilt<-findOverlaps(as(APool,"GRanges"),as(SMRT,"GRanges")) %>% as.data.frame()
BPool<-APool[-toFilt$queryHits,]

temp<-BPool
temp$start<-temp$start - 20; temp$end<-temp$end + 20
loc_41<-as(temp,"GRanges")
loc_41<-loc_41[loc_41@strand != "*" & loc_41$context == "A",]
context<-getSeq(BSgenome.Hsapiens.UCSC.hg38,loc_41,
                as.character=T)
grep(context,pattern = "N") #There is no sequences which contains "N"
context<-context[-grep(context,pattern = "N")] %>% na.omit()

sink(file="/home/owht/KIZ/data/ML_6ma_2/Dataset/Human_Negative.fasta")
for(i in 1:length(loc_41)){
  titl<-paste0(">Human_Neg_",i)
  cat(titl)
  cat("\n")
  cat(context[i])
  cat("\n")
}
sink()

####################Ara TAIR10####################################
###Generate the A pool
library(BSgenome.Athaliana.TAIR.TAIR10)
rm(list = ls())
APool<-data.frame(seqnames = rep("chr",240000),start=0,end=0,strand="*",
                  context = "N",stringsAsFactors = F)
BSgenome.Athaliana.TAIR.TAIR10@seqinfo %>% as.data.frame(stringsAsFactors=F)->SeqInfo
SMRT<-read.table("GFF/ara_TAIR10_6mA_simple.gff",
                 stringsAsFactors = F,header = F);
SMRT$V1<-paste0("araTh",SMRT$V1)
SeqInfo1<-intersect(rownames(SeqInfo),unique(SMRT$V1) %>% unique)

#Generate the random site. It's time consuming and I will change the step. It takes for about 3 days
for(i in 326000:440000){
  seqn<-sample(SeqInfo1,1)
  loci<-sample(1:SeqInfo[seqn,]$seqlengths,1)
  strr<-sample(c("+","-"),1)
  temp<-getSeq(BSgenome.Athaliana.TAIR.TAIR10,
               names=seqn,
               start=loci,
               end=loci,strand=strr,as.character = T)
  APool[i,]$seqnames<-seqn;APool[i,]$start = loci;APool[i,]$end = loci
  APool[i,]$context<-temp;APool[i,]$strand<-strr
}
rm(a,i,loci,seqn,strr,temp)

#Extract the A random sites
APool1<-APool[APool$context == "A",]

#get 1000 surrounding regions of SMRT
SMRT<-SMRT[,c(1,4,5)]
#SMRT$V1<-str_replace(SMRT$V1,pattern = "Chr0",replacement = "chr")
#SMRT$V1<-str_replace(SMRT$V1,pattern = "Chr",replacement = "chr")
colnames(SMRT)<-c("seqnames","start","end")
SMRT$strand<-"*"
SMRT$start<-SMRT$start - 500; SMRT$end<-SMRT$end + 500
colnames(SMRT)[1]<-"seqnames"

#filter Neg ranges
toFilt<-findOverlaps(as(APool1 %>% na.omit,"GRanges"),as(SMRT,"GRanges")) %>% as.data.frame()
BPool<-APool1[-toFilt$queryHits,] %>% na.omit()

temp<-BPool
temp$start<-temp$start - 20; temp$end<-temp$end + 20
loc_41<-as(temp,"GRanges")
loc_41<-loc_41[loc_41@strand != "*" & loc_41$context == "A",]
context<-getSeq(BSgenome.Athaliana.TAIR.TAIR10,loc_41,
                as.character=T)
grep(context,pattern = "N") #There is no sequences which contains "N"
context<-context[-grep(context,pattern = "N")] %>% na.omit()

sink(file="/home/owht/KIZ/data/ML_6ma_2/Dataset/Ara_Negative_TAIR10.fasta")
for(i in 1:length(loc_41)){
  titl<-paste0(">Ara_Neg_",i)
  cat(titl)
  cat("\n")
  cat(context[i])
  cat("\n")
}
sink()
