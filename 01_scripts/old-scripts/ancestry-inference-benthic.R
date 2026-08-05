#### Ancestry Inference ####

## DESCRIPTION: Inferring proportional ancestry from genotyping results from the small array

## load packages
library(tidyverse)

#### ALL GOOD SNPS ####

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
fish1<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-1.csv") %>% select(-CallRate)
fish2<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-2.csv") %>% select(-Call.Rate)
fish3<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-3.csv") %>% select(-Call.Rate)
fish4<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-4.csv") %>% select(-CallRate)
fish5<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-5.csv") %>% select(-CallRate)
fish6<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-6.csv") %>% select(-CallRate)
fish7<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-7.csv") %>% select(-CallRate)
fish8<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-8.csv") %>% select(-CallRate)
fish9<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-9.csv") %>% select(-CallRate)
fish10<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-10.csv") %>% select(-CallRate)
fish11<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-11.csv") %>% select(-CallRate)
fish<-bind_rows(fish1,fish2,fish3,fish4,fish5,fish6,fish7,fish8,fish9,fish10,fish11)
fish <- fish %>%
  unique()
colnames(fish) <- str_replace_all(colnames(fish), "X", "")
colnames(fish)[1] = "SAMPLE.ID"

## filter lakes
fish$lake<-factor(substr(fish$SAMPLE.ID,1,2))
lakes<-c("CC","LL","LP","LO")
fish<-filter(fish, lake %in% lakes)

## dupes
fish<-fish %>% group_by(SAMPLE.ID) %>%
  mutate(SAMPLE.ID = paste0(SAMPLE.ID, "_", row_number())) %>%
  ungroup()

## filter SNPs
fish<-select(fish, -all_of(bad_snps))

## replace NCs
fish[fish == "NC"] <- "NN"
fish[fish == "No Call"] <- "NN"
fish[fish == "Invalid"] <- "NN"
fish<-fish %>% mutate(across(everything(), ~ str_replace_all(., ":", "")))

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

weighting<-fishN %>% mutate(FG_W=0.25/(FG_N/(FG_N+TL_N+WB_N+WT_N)),
                            TL_W=0.25/(TL_N/(FG_N+TL_N+WB_N+WT_N)),
                            WB_W=0.25/(WB_N/(FG_N+TL_N+WB_N+WT_N)),
                            WT_W=0.25/(WT_N/(FG_N+TL_N+WB_N+WT_N))) %>% 
  select(FG_W,TL_W,WB_W,WT_W) %>% 
  as.data.frame()
rownames(weighting)<-fishN$SAMPLE.ID

weighting$FG_W<-weighting$FG_W*(1/FG_freq)
weighting$TL_W<-weighting$TL_W*(1/TL_freq)
weighting$WB_W<-weighting$WB_W*(1/WB_freq)
weighting$WT_W<-weighting$WT_W*(1/WT_freq)

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
data_raw<-data.frame(fish=rownames(fish5),
                     FG = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="FG",]))]),
                     TL = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="TL",]))]),
                     WB = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="WB",]))]),
                     WT = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="WT",]))]),
                     FG_homo = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="FG",]))]>1),
                     TL_homo = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="TL",]))]>1),
                     WB_homo = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="WB",]))]>1),
                     WT_homo = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="WT",]))]>1))

## weighting
rownames(data_raw)==rownames(weighting)
data<-data_raw[,2:5]*weighting
data_raw_w<-data_raw[,2:5]*weighting

## proportions
data<-data %>% mutate(FG_prop=FG/(FG+TL+WB+WT),
                      TL_prop=TL/(FG+TL+WB+WT),
                      WB_prop=WB/(FG+TL+WB+WT),
                      WT_prop=WT/(FG+TL+WB+WT))

## fractions based on max generations
maxGen<-2
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
data_final<-as.data.frame(t(apply(data[,5:8]*cats,1,sum.round)))/cats
colnames(data_final)=colnames(data[,1:4])
rowSums(data_final)
data_final[rowSums(data_final)==0,]

## adding back homozygosity and raw counts
data_final<-cbind(data_final, data_raw[,6:9], data_raw[,2:5])
data_final$fish<-rownames(data_final)
colnames(data_final)[9:12]<-c("FG_raw","TL_raw","WB_raw","WT_raw")

write.csv(data_final, file = "./02_tidy-data/benthic-ancestry.csv",quote = F, row.names = F)











## df for checking homozygosity
no_loon<-data_final %>% mutate(lake=factor(substr(fish,1,2))) %>% 
  filter(lake != "LO")

## homozygosity of pure types
pures<-no_loon[rowSums(no_loon[1:4]==1)>0,]
pures[rowSums(pures[,5:8]>1)==0,] ## LL51 might be an error
rownames(pures[rowSums(pures[,5:8]>1)==0,])
##data_final[rownames(pures[rowSums(pures[,5:8]>1)==0,]),1:4]<-c(NA,NA,NA,NA)

## homozygosity of 50/50s
fifties<-no_loon[rowSums(no_loon[1:4]==0.5)==2,]
fifties2<-fifties[rowSums(fifties[,5:8]>2)==1,]
fifties2[rowSums(fifties2[,5:8]>0)<2,] ## none

## homozygosity of 75/25
backs<-no_loon[rowSums(no_loon[,1:4]==0.75)>0,]
backs[rowSums(backs[,5:8]>1)==0,] ## LP02, CC36, CC26, LP62 might be 50/50
backs[rowSums(backs[,5:8]>1)==2,] ## LL74, LL97 might be 50/50

## homozygosity of 50/25/25
f2<-no_loon[rowSums(no_loon[,1:4]==0.25)>1,]
outcross<-f2[rowSums(f2[,1:4]==0.50)>0,]
outcross[rowSums(outcross[,5:8]>1)>1,] ## none

## homozygosity of 25/25/25/25
truef2<-f2[rowSums(f2[,1:4]==0.25)==4,]
truef2[rowSums(truef2[,5:8]>1)>0,] ## none

## checking sums
no_loon2<-filter(data_raw_w,!grepl("LO", rownames(data_raw_w)))
boxplot(rowSums(no_loon2[1:4]))
data_final[names(boxplot.stats(rowSums(no_loon2[1:4]))$out),]
## data_final[names(boxplot.stats(rowSums(no_loon1[1:4]))$out),1:4]<-c(NA,NA,NA,NA)

## checking
data_old<-read.csv("./02_tidy-data/final-ancestry-allsnps.csv", stringsAsFactors = T)
data_old<-filter(data_old, lake %in% c("CC","LL","LP"))
data_old[order(data_old$fish),3:6]
data_f<-data_final %>% mutate(lake = substr(fish,1,2)) %>% 
  filter(lake != "LO")
data_f[order(data_f$fish),1:4]
sum(data_old[order(data_old$fish),3:6] != data_f[order(data_f$fish),1:4],na.rm=T)

#### OLD ####

## make sure all fish snps are in population snp data 
fish_snps <- str_replace_all(colnames(fish), "X", "")
unique(popsnp$SNP.ID)
which(!fish_snps %in% unique(popsnp$SNP.ID))

## subset to columns we care about 
popsnp <- select(popsnp, SNP.ID, pop, Base_A, pop_freq)

## transpose the fish data so that each row = 1 snp, each column = 1 fish
fish = t(fish)
colnames(fish) = fish[1,]
fish = as.data.frame(fish)
fish$SAMPLE.ID = row.names(fish)
row.names(fish) = NULL
colnames(fish) = fish[1, ]
fish <- fish[-1, ]
fish = select(fish, SAMPLE.ID, everything())
colnames(fish)[1] = "SNP.ID"

## check duplicated sample IDs:
colnames(fish)[which(duplicated(colnames(fish)))]

## split each column into 2:
length(fish[1,c(2:ncol(fish))])
fishlist = as.list(fish[,c(2:ncol(fish))])

first_char = sapply(fishlist, FUN = function(x) {
  return(substring(x, 1, 1))
})

sec_char = sapply(fishlist, FUN = function(x) {
  return(substring(x, 2, 2))
})

colnames(first_char) <- paste(colnames(first_char), "_1char", sep = "")
colnames(sec_char) <- paste(colnames(sec_char), "_2char", sep = "")

first_char = as.data.frame(first_char)
sec_char = as.data.frame(sec_char)

first_char$SNP.ID = fish$SNP.ID
sec_char$SNP.ID = fish$SNP.ID

## join: 
fishsplit = full_join(first_char, sec_char)%>%
  select(SNP.ID, everything())

## join with snps:
popsnp$SNP.ID<-as.character(popsnp$SNP.ID)
mega <- left_join(popsnp, fishsplit)

## replace with 1 or 0 based on if they match, weighted by frequency
for (i in 1:(ncol(mega)-4)){
  
  data<-mega[,(4+i)]
  
  for (n in 1:length(data)){
    
    if (mega$Base_A[n]==data[n]){mega[n,(4+i)]<-1/mega$pop_freq[n]}
    if (identical(mega$Base_A[n],data[n])==FALSE) {mega[n,(4+i)]<-0}
  }
}

## sum columns for each sample
combined.mega<-mega[,1:4]
for (i in 1:length(unique(gsub("_.char","",names(mega)[5:ncol(mega)])))){

  col=unique(gsub("_.char","",names(mega)[5:ncol(mega)]))[i]
  combined.mega[,i+4]<-as.numeric(mega[,c(paste0(col,"_1char"))])+as.numeric(mega[,c(paste0(col,"_2char"))])
  names(combined.mega)[i+4]<-col
}

## exploring SNPs
# combined.mega$sum<-rowSums(combined.mega[,5:ncol(combined.mega)])
# select(combined.mega, freq, pop, sum)
# ggplot(combined.mega, aes(x=freq, y=sum, label=SNP.ID))+
#   geom_point()+
#   geom_text(hjust=-0.25)+
#   theme_bw()+
#   facet_wrap(~pop)
# filter(combined.mega, sum==0)$SNP.ID

## just TL fish
# ggplot(filter(combined.mega, pop=="WT"), aes(x=freq, y=sum, label=SNP.ID))+
#   geom_point()+
#   geom_text(hjust=-0.25)+
#   theme_bw()

## sum each pop for each sample
data=matrix(data=NA,ncol=(ncol(combined.mega)-3),nrow=length(unique(combined.mega$pop)))
colnames(data)=c("pop",names(combined.mega[,5:ncol(combined.mega)]))
for (i in 1:length(unique(combined.mega$pop))){

  pop=unique(combined.mega$pop)[i]
  data[i,]= cbind(pop,t(colSums(combined.mega[which(combined.mega$pop==pop),c(5:ncol(combined.mega))])))
}

## cleaning up dataframe
data<-as.data.frame(data)
data<-t(data)
colnames(data)<-data[1,]
data<-as.data.frame(data[-1,])
data <- data %>% mutate_if(is.character, as.numeric)
rowSums(data)
data[rowSums(data)==0,]

## weight by number of SNPs per pop
data_raw<-data
colnames(data_raw)<-c("FG_raw","TL_raw","WB_raw","WT_raw")
sum(rownames(data)!=rownames(weighting))
dim(data)==dim(weighting)
data<-data*weighting

## fixing errors
boxplot(rowSums(filter(data,!grepl("LO", rownames(data)))))
boxplot.stats(rowSums(filter(data,!grepl("LO", rownames(data)))))$out
data[names(boxplot.stats(rowSums(filter(data,!grepl("LO", rownames(data)))))$out),]
data[!grepl("LO", rownames(data)) & rowSums(data)<16,]<-c(NA,NA,NA,NA)
summary(data)

## proportions
data$FG_prop<-data$FG/(data$FG+data$TL+data$WB+data$WT)
data$TL_prop<-data$TL/(data$FG+data$TL+data$WB+data$WT)
data$WB_prop<-data$WB/(data$FG+data$TL+data$WB+data$WT)
data$WT_prop<-data$WT/(data$FG+data$TL+data$WB+data$WT)

## fractions based on max generations
maxGen<-2
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
data_final<-as.data.frame(t(apply(data[,5:8]*cats,1,sum.round)))/cats
colnames(data_final)=colnames(data[,1:4])
rowSums(data_final)
data_final[rowSums(data_final)==0,]

## adding back the raw scores
data_final<-cbind(data_final, data_raw[,1:4], fish2[,c("FG_N","TL_N","WB_N","WT_N")])
data_final$fish<-rownames(data_final)

## refixing errors
data_final[!grepl("LO", rownames(data_final)) & rowSums(data[1:4])<25,]
data_final[!grepl("LO", rownames(data_final)) & rowSums(data[1:4])<25,1:4]<-c(NA,NA,NA,NA)

## save results
write.csv(data_final, file = "02_tidy-data/benthic-ancestry-allsnps.csv",
          quote = F, row.names = F)


#### 12 BEST SNPS ####

## read in pop snp data
popsnp <- read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-final-snps.csv")

## choose best snps
FG_snps<-c("16","22","20","17","21","2","7","23","3","6","4","19")
TL_snps<-c("25","26","27","29","31","33","34","39","40","43","45_2","47")
WB_snps<-c("66","53","54","68","69","60_2","72","57","49","51","67","56_2")
WT_snps<-c("84","82","83","85","95","76","73","77","89","88","79","81")
best_snps<-c(FG_snps,TL_snps,WB_snps,WT_snps)
popsnp<-filter(popsnp, SNP.ID %in% best_snps)
## 48 snps with average freq = 83.5

## average freq by pop
popfreq<-popsnp %>% group_by(pop) %>% 
  summarise(pop_freq=mean(freq))
popsnp<-mutate(popsnp, pop_freq = case_when(pop == "FG" ~ popfreq[[1,2]],
                                            pop == "TL" ~ popfreq[[2,2]],
                                            pop == "WB" ~ popfreq[[3,2]],
                                            pop == "WT" ~ popfreq[[4,2]]))

## read in fish snp data
fish1<-read_csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-1.csv")
fish2<-read_csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-2.csv")
fish3<-read_csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-3.csv")

fish<-bind_rows(fish1,fish2,fish3)

fish <- fish %>%
  unique()
colnames(fish) <- str_replace_all(colnames(fish), "X", "")
colnames(fish)[1] = "SAMPLE.ID"

## filter lakes
fish$lake<-factor(substr(fish$SAMPLE.ID,1,2))
lakes<-c("CC","LL","LP","LO")
fish<-filter(fish, lake %in% lakes)

## replace NCs
fish[fish == "NC"] <- "NN"

## make sure all fish snps are in population snp data 
fish_snps <- str_replace_all(colnames(fish), "X", "")
unique(popsnp$SNP.ID)
which(!fish_snps %in% unique(popsnp$SNP.ID))

## subset to columns we care about 
popsnp <- select(popsnp, SNP.ID, pop, Base_A, pop_freq)

## transpose the fish data so that each row = 1 snp, each column = 1 fish
fish = t(fish)
colnames(fish) = fish[1,]
fish = as.data.frame(fish)
fish$SAMPLE.ID = row.names(fish)
row.names(fish) = NULL
colnames(fish) = fish[1, ]
fish <- fish[-1, ]
fish = select(fish, SAMPLE.ID, everything())
colnames(fish)[1] = "SNP.ID"

## check duplicated sample IDs:
colnames(fish)[which(duplicated(colnames(fish)))]

## split each column into 2:
length(fish[1,c(2:ncol(fish))])
fishlist = as.list(fish[,c(2:ncol(fish))])

first_char = sapply(fishlist, FUN = function(x) {
  return(substring(x, 1, 1))
})

sec_char = sapply(fishlist, FUN = function(x) {
  return(substring(x, 2, 2))
})

colnames(first_char) <- paste(colnames(first_char), "_1char", sep = "")
colnames(sec_char) <- paste(colnames(sec_char), "_2char", sep = "")

first_char = as.data.frame(first_char)
sec_char = as.data.frame(sec_char)

first_char$SNP.ID = fish$SNP.ID
sec_char$SNP.ID = fish$SNP.ID

## join: 
fishsplit = full_join(first_char, sec_char)%>%
  select(SNP.ID, everything())

## join with snps:
popsnp$SNP.ID<-as.character(popsnp$SNP.ID)
mega <- left_join(popsnp, fishsplit)

## replace with 1 or 0 based on if they match, weighted by frequency
for (i in 1:(ncol(mega)-4)){
  
  data<-mega[,(4+i)]
  
  for (n in 1:length(data)){
    
    if (mega$Base_A[n]==data[n]){mega[n,(4+i)]<-1/mega$pop_freq[n]}
    if (identical(mega$Base_A[n],data[n])==FALSE) {mega[n,(4+i)]<-0}
  }
}

## sum columns for each sample
combined.mega<-mega[,1:4]
for (i in 1:length(unique(gsub("_.char","",names(mega)[5:ncol(mega)])))){
  
  col=unique(gsub("_.char","",names(mega)[5:ncol(mega)]))[i]
  combined.mega[,i+4]<-as.numeric(mega[,c(paste0(col,"_1char"))])+as.numeric(mega[,c(paste0(col,"_2char"))])
  names(combined.mega)[i+4]<-col
}

## exploring SNPs
# combined.mega$sum<-rowSums(combined.mega[,5:ncol(combined.mega)])
# select(combined.mega, freq, pop, sum)
# ggplot(combined.mega, aes(x=freq, y=sum, label=SNP.ID))+
#   geom_point()+
#   geom_text(hjust=-0.25)+
#   theme_bw()+
#   facet_wrap(~pop)

## sum each pop for each sample
data=matrix(data=NA,ncol=(ncol(combined.mega)-3),nrow=length(unique(combined.mega$pop)))
colnames(data)=c("pop",names(combined.mega[,5:ncol(combined.mega)]))
for (i in 1:length(unique(combined.mega$pop))){
  
  pop=unique(combined.mega$pop)[i]
  data[i,]= cbind(pop,t(colSums(combined.mega[which(combined.mega$pop==pop),c(5:ncol(combined.mega))])))
}

## cleaning up dataframe
data<-as.data.frame(data)
data<-t(data)
colnames(data)<-data[1,]
data<-as.data.frame(data[-1,])
data <- data %>% mutate_if(is.character, as.numeric)

## fixing errors
boxplot(rowSums(filter(data,!grepl("LO", rownames(data)))))
boxplot.stats(rowSums(filter(data,!grepl("LO", rownames(data)))))$out
data[names(boxplot.stats(rowSums(filter(data,!grepl("LO", rownames(data)))))$out),]

data[!grepl("LO", rownames(data)) & rowSums(data)<10,]<-c(NA,NA,NA,NA)
summary(data)

## proportions
data$FG_prop<-data$FG/(data$FG+data$TL+data$WB+data$WT)
data$TL_prop<-data$TL/(data$FG+data$TL+data$WB+data$WT)
data$WB_prop<-data$WB/(data$FG+data$TL+data$WB+data$WT)
data$WT_prop<-data$WT/(data$FG+data$TL+data$WB+data$WT)

## fractions based on max generations
maxGen<-2
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
data_final<-as.data.frame(t(apply(data[,5:8]*cats,1,sum.round)))/cats
colnames(data_final)=colnames(data[,1:4])
rowSums(data_final)
data_final[rowSums(data_final)==0,]

## adding back the raw scores
data_final<-cbind(data_final, data[,1:4])
data_final$fish<-rownames(data_final)
data_final <- data_final[, c(9, 1:4, 5:8)]
colnames(data_final)[6:9]<-c("FG_raw","TL_raw","WB_raw","WT_raw")

## save results
write.csv(data_final, file = "02_tidy-data/benthic-ancestry-bestsnps.csv",
          quote = F, row.names = F)
