#### Alexis Loon Ancestry ####

## load packages
library(tidyverse)

#### BENTHIC ####

## read in pop snp data
popsnp <- read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-final-snps.csv")

## get rid of bad snps
FG_bad<-c("12")
TL_bad<-c("28","33","41","36")
WB_bad<-c("50","62","61","70","52","71")
WT_bad<-c("86","74","92","83","78","96","91","90","80","93","94")
bad_snps<-c(FG_bad,TL_bad,WB_bad,WT_bad)
popsnp<-filter(popsnp, !SNP.ID %in% bad_snps)
## 74 snps with average freq = 80

## average freq by pop
popfreq<-popsnp %>% group_by(pop) %>% 
  summarise(pop_freq=mean(freq))
popsnp<-mutate(popsnp, pop_freq = case_when(pop == "FG" ~ popfreq[[1,2]],
                                            pop == "TL" ~ popfreq[[2,2]],
                                            pop == "WB" ~ popfreq[[3,2]],
                                            pop == "WT" ~ popfreq[[4,2]]))
FG_freq<-popfreq[[1,2]]
TL_freq<-popfreq[[2,2]]
WB_freq<-popfreq[[3,2]]
WT_freq<-popfreq[[4,2]]

## read in fish snp data
fish1<-read_csv("./00_raw-data/02_snp-array-data/loon-data/loon-benthic-1.csv")
fish2<-read_csv("./00_raw-data/02_snp-array-data/loon-data/loon-benthic-2.csv")
fish<-bind_rows(fish1,fish2)
fish <- fish %>%
  unique()
colnames(fish) <- str_replace_all(colnames(fish), "X", "")
colnames(fish)[1] = "SAMPLE.ID"

## filter SNPs
fish<-select(fish, -all_of(bad_snps))

## replace NCs
fish[fish == "NC"] <- "NN"

## count NN's
FG_snps<-filter(popsnp, pop=="FG")$SNP.ID
TL_snps<-filter(popsnp, pop=="TL")$SNP.ID
WB_snps<-filter(popsnp, pop=="WB")$SNP.ID
WT_snps<-filter(popsnp, pop=="WT")$SNP.ID

fishN<-fish
fishN$FG_N<-length(FG_snps)-rowSums(fishN[,FG_snps]=="NN")
fishN$TL_N<-length(TL_snps)-rowSums(fishN[,TL_snps]=="NN")
fishN$WB_N<-length(WB_snps)-rowSums(fishN[,WB_snps]=="NN")
fishN$WT_N<-length(WT_snps)-rowSums(fishN[,WT_snps]=="NN")

ben_N<-fishN %>% select(SAMPLE.ID,FG_N,TL_N,WB_N,WT_N) %>% as.data.frame()

## scores
fish0<-t(fish)
colnames(fish0)<-fish0[1,]
fish1<-as.data.frame(fish0[-1,])
fish1<-fish1[1:74,]
fish2<-fish1[order(match(rownames(fish1), popsnp$SNP.ID)),]
fish3<-data.frame(t(fish2))
colnames(fish3)<-rownames(fish2)

## list of SNP bases
snps<-popsnp$Base_A

## scores for bases
for (i in 1:ncol(fish3)) {
  fish3[[paste0(names(fish3)[i], "_score")]] <- sapply(fish3[[i]], function(x) str_count(x, snps[i]))
}
fish5<-select(fish3,'1_score':'95_score')

## summing scores and homozygosity
ben_data_raw<-data.frame(fish=rownames(fish5),
                     FG = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="FG",]))]),
                     TL = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="TL",]))]),
                     WB = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="WB",]))]),
                     WT = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="WT",]))]))

#### LIMNETIC ####

## read in pop snp data
popsnp <- read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-final-snps.csv")

## get rid of bad SNPs
LG_bad<-c("116","99","97","104")
SL_bad<-c("126","130","143","136","144","127","142")
SR_bad<-c()
WK_bad<-c("171")
bad_snps<-c(LG_bad,SL_bad,SR_bad,WK_bad)
popsnp<-filter(popsnp, !SNP.ID %in% bad_snps)
## 84 snps with average freq = 83

## average freq by pop
popfreq<-popsnp %>% group_by(pop) %>% 
  summarise(pop_freq=mean(freq))
popsnp<-mutate(popsnp, pop_freq = case_when(pop == "LG" ~ popfreq[[1,2]],
                                            pop == "SL" ~ popfreq[[2,2]],
                                            pop == "SR" ~ popfreq[[3,2]],
                                            pop == "WK" ~ popfreq[[4,2]]))
LG_freq<-popfreq[[1,2]]
SL_freq<-popfreq[[2,2]]
SR_freq<-popfreq[[3,2]]
WK_freq<-popfreq[[4,2]]

## read in fish snp data
fish1<-read_csv("./00_raw-data/02_snp-array-data/loon-data/loon-limnetic-1.csv")
fish2<-read_csv("./00_raw-data/02_snp-array-data/loon-data/loon-limnetic-2.csv")
fish<-bind_rows(fish1,fish2)

fish <- fish %>%
  unique()
colnames(fish) <- str_replace_all(colnames(fish), "X", "")
colnames(fish)[1] = "SAMPLE.ID"

## filter SNPs
fish<-select(fish, -all_of(bad_snps))

## replace NCs
fish[fish == "NC"] <- "NN"

## count NN's
LG_snps<-filter(popsnp, pop=="LG")$SNP.ID
SL_snps<-filter(popsnp, pop=="SL")$SNP.ID
SR_snps<-filter(popsnp, pop=="SR")$SNP.ID
WK_snps<-filter(popsnp, pop=="WK")$SNP.ID

fish2<-fish
fish2$LG_N<-length(LG_snps)-rowSums(fish2[,LG_snps]=="NN")
fish2$SL_N<-length(SL_snps)-rowSums(fish2[,SL_snps]=="NN")
fish2$SR_N<-length(SR_snps)-rowSums(fish2[,SR_snps]=="NN")
fish2$WK_N<-length(WK_snps)-rowSums(fish2[,WK_snps]=="NN")

lim_N<-fish2 %>% select(SAMPLE.ID,LG_N,SL_N,SR_N,WK_N) %>% as.data.frame()

## scores
fish0<-t(fish)
colnames(fish0)<-fish0[1,]
fish1<-as.data.frame(fish0[-1,])
fish1<-fish1[1:84,]
fish2<-fish1[order(match(rownames(fish1), popsnp$SNP.ID)),]
fish3<-data.frame(t(fish2))
colnames(fish3)<-rownames(fish2)

## list of SNP bases
snps<-popsnp$Base_A

## scores for bases
for (i in 1:ncol(fish3)) {
  fish3[[paste0(names(fish3)[i], "_score")]] <- sapply(fish3[[i]], function(x) str_count(x, snps[i]))
}
fish5<-select(fish3,'98_score':'192_score')

## summing scores and homozygosity
lim_data_raw<-data.frame(fish=rownames(fish5),
                     LG = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="LG",]))]),
                     SL = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="SL",]))]),
                     SR = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="SR",]))]),
                     WK = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="WK",]))]))

#### COMBINING ####

## combine data
colnames(ben_N)[1]<-"fish"
colnames(lim_N)[1]<-"fish"
data_raw<-full_join(ben_data_raw,lim_data_raw)
N<-full_join(ben_N,lim_N)

## adjust counts by number of SNPs and frequencies
full_data<-full_join(data_raw,N) %>% 
  mutate(FG=FG*(0.125/(FG_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N)))*(1/FG_freq),
         TL=TL*(0.125/(TL_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N)))*(1/TL_freq),
         WB=WB*(0.125/(WB_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N)))*(1/WB_freq),
         WT=WT*(0.125/(WT_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N)))*(1/WT_freq),
         LG=LG*(0.125/(LG_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N)))*(1/LG_freq),
         SL=SL*(0.125/(SL_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N)))*(1/SL_freq),
         SR=SR*(0.125/(SR_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N)))*(1/SR_freq),
         WK=WK*(0.125/(WK_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N)))*(1/WK_freq)) %>% 
  select(fish,FG,TL,WB,WT,LG,SL,SR,WK)

## check outliers in sums
boxplot(rowSums(full_data[,2:9]))
boxplot.stats(rowSums(full_data[,2:9]))$out

## fix one outlier
full_data[rowSums(full_data[,2:9])<20,2:9]<-rep(NA,8)

## proportions
data<-full_data %>% mutate(FG_prop=FG/(FG+TL+WB+WT+LG+SL+SR+WK),
                               TL_prop=TL/(FG+TL+WB+WT+LG+SL+SR+WK),
                               WB_prop=WB/(FG+TL+WB+WT+LG+SL+SR+WK),
                               WT_prop=WT/(FG+TL+WB+WT+LG+SL+SR+WK),
                               LG_prop=LG/(FG+TL+WB+WT+LG+SL+SR+WK),
                               SL_prop=SL/(FG+TL+WB+WT+LG+SL+SR+WK),
                               SR_prop=SR/(FG+TL+WB+WT+LG+SL+SR+WK),
                               WK_prop=WK/(FG+TL+WB+WT+LG+SL+SR+WK)) %>% 
  select(fish,FG_prop,TL_prop,WB_prop,WT_prop,LG_prop,SL_prop,SR_prop,WK_prop)

## rounding loon
maxGen<-3
cats<-2^maxGen

## function for rounding while preserving the sum of the rows
sum.round <- function(x, digits = 0) {
  x<-as.numeric(x)
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

## rounding the proportions and checking sums
data[is.na(data)] <- 0
data_final<-as.data.frame(t(apply(data[,2:9]*cats,1,sum.round)))/cats
colnames(data_final)=c("FG","TL","WB","WT","LG","SL","SR","WK")
rowSums(data_final)
data_final[rowSums(data_final)==0,]<-c(rep(NA,8))
data_final$fish<-data$fish
data_final<-select(data_final, fish, FG:WK)

data_raw<-left_join(data_final,data)

## success
success<-100*colSums(data_final[,2:9], na.rm = T)/(nrow(data_final)-1)

## writing csv
write.csv(data_final, file = "./02_tidy-data/alexis-loon-ancestry.csv",
          quote = F, row.names = F)
write.csv(data_raw, file = "./02_tidy-data/alexis-loon-ancestry-raw.csv",
          quote = F, row.names = F)



