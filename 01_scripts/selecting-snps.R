#### SNP Selection ####

## DESCRIPTION: Choosing SNPs to include on the array for inferring ancestry
## AUTHORS: Lucas Eckert

## packages
library(tidyverse)

#importing data
data<-read.csv("PopSpecificSNPs.csv", header = T)
data$LG<-as.factor(data$LG)
data<-subset(data, LG!="chrUn")
data<-subset(data, LG!="chrM")
data<-subset(data, LG!="chrXIX")

#calculate the frequency of all the SNPs
#take the two highest frequency SNPs per chromosome per population

#FG
FG<-subset(data, FG_N_A!=0)
FG<-subset(FG, FG_Nreads>56)
FG$pop<-as.factor("FG")
FG$freq<-FG$FG_N_A/FG$FG_Nreads
FGbyChr<-FG %>% group_by(LG)
topFG<-FGbyChr %>% slice_max(freq, n=3)

#LG
LG<-subset(data, LG_N_A!=0)
LG<-subset(LG, LG_Nreads>56)
LG$pop<-as.factor("LG")
LG$freq<-LG$LG_N_A/LG$LG_Nreads
LGbyChr<-LG %>% group_by(LG)
topLG<-LGbyChr %>% slice_max(freq, n=3) 

#SR
SR<-subset(data, SR_N_A!=0)
SR<-subset(SR, SR_Nreads>56)
SR$pop<-as.factor("SR")
SR$freq<-SR$SR_N_A/SR$SR_Nreads
SRbyChr<-SR %>% group_by(LG)
topSR<-SRbyChr %>% slice_max(freq, n=3) 

#SL
SL<-subset(data, SL_N_A!=0)
SL<-subset(SL, SL_Nreads>56)
SL$pop<-as.factor("SL")
SL$freq<-SL$SL_N_A/SL$SL_Nreads
SLbyChr<-SL %>% group_by(LG)
topSL<-SLbyChr %>% slice_max(freq, n=3)

#TL
TL<-subset(data, TL_N_A!=0)
TL<-subset(TL, TL_Nreads>56)
TL$pop<-as.factor("TL")
TL$freq<-TL$TL_N_A/TL$TL_Nreads
TLbyChr<-TL %>% group_by(LG)
topTL<-TLbyChr %>% slice_max(freq, n=3) 

#WB
WB<-subset(data, WB_N_A!=0)
WB<-subset(WB, WB_Nreads>56)
WB$pop<-as.factor("WB")
WB$freq<-WB$WB_N_A/WB$WB_Nreads
WBbyChr<-WB %>% group_by(LG)
topWB<-WBbyChr %>% slice_max(freq, n=3) 

#WT
WT<-subset(data, WT_N_A!=0)
WT<-subset(WT, WT_Nreads>56)
WT$pop<-as.factor("WT")
WT$freq<-WT$WT_N_A/WT$WT_Nreads
WTbyChr<-WT %>% group_by(LG)
topWT<-WTbyChr %>% slice_max(freq, n=3) 

#WK
WK<-subset(data, WK_N_A!=0)
WK<-subset(WK, WK_Nreads>56)
WK$pop<-as.factor("WK")
WK$freq<-WK$WK_N_A/WK$WK_Nreads
WKbyChr<-WK %>% group_by(LG)
topWK<-WKbyChr %>% slice_max(freq, n=3) 

#put em all together
topSNPs<-rbind(topFG,topLG,topSR,topSL,topTL,topWB,topWT,topWK)
with(topSNPs, table(pop,LG))
write.csv(topSNPs, file ="topPopSNPsFiltered.csv")

#calculate the mean and sd of total reads, used for determing appropriate minimum reads
tot.reads<-c(data$FG_Nreads,data$LG_Nreads,data$SR_Nreads,data$SL_Nreads,data$TL_Nreads,data$WB_Nreads,data$WT_Nreads,data$WK_Nreads)
hist(tot.reads, breaks = 100)
mean(tot.reads)-2*sd(tot.reads)
