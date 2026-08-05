#### Ancestry Inference ####

## DESCRIPTION: Inferring proportional ancestry from genotyping results from the small array

## load packages
library(tidyverse)

#### ALL GOOD SNPS ####

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
fish1<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-1.csv") %>% select(-CallRate)
fish2<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-2.csv") %>% select(-CallRate)
fish3<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-3.csv") %>% select(-CallRate)
fish4<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-4.csv") %>% select(-CallRate)
fish5<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-5.csv") %>% select(-CallRate)
fish6<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-6.csv") %>% select(-CallRate)
fish7<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-7.csv") %>% select(-CallRate)
fish8<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-8.csv") %>% select(-CallRate)
fish9<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-9.csv") %>% select(-CallRate)
fish10<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-10.csv") %>% select(-CallRate)
fish11<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-11.csv") %>% select(-CallRate)
fish12<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-12.csv") %>% select(-CallRate)
fish13<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-13.csv") %>% select(-CallRate)
fish14<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-14.csv") %>% select(-CallRate)
fish15<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-15.csv") %>% select(-CallRate)
fish16<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-16.csv") %>% select(-CallRate)

fish<-bind_rows(fish1,fish2,fish3,fish4,fish5,fish6,fish7,fish8,fish9,fish10,fish11,fish12,fish13,fish14,fish15,fish16)

fish <- fish %>%
  unique()
colnames(fish) <- str_replace_all(colnames(fish), "X", "")
colnames(fish)[1] = "SAMPLE.ID"

## filter lakes
fish$lake<-factor(substr(fish$SAMPLE.ID,1,2))
lakes<-c("CL","FL","HL","RL","LO")
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
LG_snps<-filter(popsnp, pop=="LG")$SNP.ID
SL_snps<-filter(popsnp, pop=="SL")$SNP.ID
SR_snps<-filter(popsnp, pop=="SR")$SNP.ID
WK_snps<-filter(popsnp, pop=="WK")$SNP.ID

fish2<-fish
fish2$LG_N<-length(LG_snps)-rowSums(fish2[,LG_snps]=="NN")
fish2$SL_N<-length(SL_snps)-rowSums(fish2[,SL_snps]=="NN")
fish2$SR_N<-length(SR_snps)-rowSums(fish2[,SR_snps]=="NN")
fish2$WK_N<-length(WK_snps)-rowSums(fish2[,WK_snps]=="NN")

weighting<-fish2 %>% mutate(LG_W=0.25/(LG_N/(LG_N+SL_N+SR_N+WK_N)),
                            SL_W=0.25/(SL_N/(LG_N+SL_N+SR_N+WK_N)),
                            SR_W=0.25/(SR_N/(LG_N+SL_N+SR_N+WK_N)),
                            WK_W=0.25/(WK_N/(LG_N+SL_N+SR_N+WK_N))) %>% 
  select(LG_W,SL_W,SR_W,WK_W) %>% 
  as.data.frame()
rownames(weighting)<-fish2$SAMPLE.ID

weighting$LG_W<-weighting$LG_W*(1/LG_freq)
weighting$SL_W<-weighting$SL_W*(1/SL_freq)
weighting$SR_W<-weighting$SR_W*(1/SR_freq)
weighting$WK_W<-weighting$WK_W*(1/WK_freq)

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
data_raw<-data.frame(fish=rownames(fish5),
                     LG = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="LG",]))]),
                     SL = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="SL",]))]),
                     SR = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="SR",]))]),
                     WK = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="WK",]))]),
                     LG_homo = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="LG",]))]>1),
                     SL_homo = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="SL",]))]>1),
                     SR_homo = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="SR",]))]>1),
                     WK_homo = rowSums(fish5[,as.numeric(rownames(popsnp[popsnp$pop=="WK",]))]>1))

## weighting
rownames(data_raw)==rownames(weighting)
data<-data_raw[,2:5]*weighting
data_raw_w<-data_raw[,2:5]*weighting

## proportions
data<-data %>% mutate(LG_prop=LG/(LG+SL+SR+WK),
                      SL_prop=SL/(LG+SL+SR+WK),
                      SR_prop=SR/(LG+SL+SR+WK),
                      WK_prop=WK/(LG+SL+SR+WK))

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
colnames(data_final)[9:12]<-c("LG_raw","SL_raw","SR_raw","WK_raw")

write.csv(data_final, file = "./02_tidy-data/limnetic-ancestry.csv",quote = F,row.names = F)


## df for checking homozygosity
no_loon<-data_final %>% mutate(lake=factor(substr(fish,1,2))) %>% 
  filter(lake != "LO")

## homozygosity of pure types
pures<-no_loon[rowSums(no_loon[1:4]==1)>0,]
pures[rowSums(pures[,5:8]>1)==0,] ## HL08, HL27, HL33 are errors
rownames(pures[rowSums(pures[,5:8]>1)==0,])
data_final[rownames(pures[rowSums(pures[,5:8]>1)==0,]),1:4]<-c(NA,NA,NA,NA)

## homozygosity of 50/50s
fifties<-no_loon[rowSums(no_loon[1:4]==0.5)==2,]
fifties2<-fifties[rowSums(fifties[,5:8]>2)==1,]
fifties2[rowSums(fifties2[,5:8]>0)<2,] ## CL20, HL85, RL12, RL16 might be 75/25
data_raw_w[rownames(fifties2[rowSums(fifties2[,5:8]>0)<2,]),]
##data_final[c("CL21L020","HL21L085"),1:4] <- c(0.25,0,0.75,0)
##data_final[c("RL21L012","RL21L016"),1:4] <- c(0.25,0,0,0.75)

## homozygosity of 75/25
backs<-no_loon[rowSums(no_loon[,1:4]==0.75)>0,]
backs[rowSums(backs[,5:8]>1)==0,] ## CL84, HL44, HL87, RL01, RL73 might be 50/50
backs[rowSums(backs[,5:8]>1)==2,] ## HL84

## homozygosity of 50/25/25
f2<-no_loon[rowSums(no_loon[,1:4]==0.25)>1,]
outcross<-f2[rowSums(f2[,1:4]==0.50)>0,]
outcross[rowSums(outcross[,5:8]>1)>1,] ## HL77, RL52, RL90, RL100

## homozygosity of 25/25/25/25
truef2<-f2[rowSums(f2[,1:4]==0.25)==4,]
truef2[rowSums(truef2[,5:8]>1)>0,]

## checking sums
boxplot(rowSums(no_loon[9:12]))
names(boxplot.stats(rowSums(no_loon[9:12]))$out)
data_final[names(boxplot.stats(rowSums(no_loon[9:12]))$out),1:4]<-c(NA,NA,NA,NA)

## checking
data_old<-read.csv("./02_tidy-data/final-ancestry-allsnps.csv", stringsAsFactors = T)
data_old<-filter(data_old, lake %in% c("CL","FL","HL","RL"))
data_old[order(data_old$fish),7:10]

data_f<-data_final %>% mutate(lake = substr(fish,1,2)) %>% 
  filter(lake != "LO")
data_f[order(data_f$fish),1:4]

sum(data_old[order(data_old$fish),7:10] != data_f[order(data_f$fish),1:4],na.rm=T)












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
# select(combined.mega, SNP.ID,freq, pop, sum)[order(combined.mega$sum),]
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
rowSums(data)
data[rowSums(data)==0,]

## weight by number of SNPs per pop
data_raw<-data
colnames(data_raw)<-c("LG_raw","SL_raw","SR_raw","WK_raw")
sum(rownames(data)!=rownames(weighting))
dim(data)==dim(weighting)
data<-data*weighting

## fixing errors
boxplot(rowSums(filter(data,!grepl("LO", rownames(data)))))
boxplot.stats(rowSums(filter(data,!grepl("LO", rownames(data)))))$out
data[names(boxplot.stats(rowSums(filter(data,!grepl("LO", rownames(data)))))$out),]
data[!grepl("LO", rownames(data)) & rowSums(data)<20,]<-c(NA,NA,NA,NA)
summary(data)

## proportions
data$LG_prop<-data$LG/(data$LG+data$SL+data$SR+data$WK)
data$SL_prop<-data$SL/(data$LG+data$SL+data$SR+data$WK)
data$SR_prop<-data$SR/(data$LG+data$SL+data$SR+data$WK)
data$WK_prop<-data$WK/(data$LG+data$SL+data$SR+data$WK)

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
data_final<-cbind(data_final, data_raw[,1:4], fish2[,c("LG_N","SL_N","SR_N","WK_N")])
data_final$fish<-rownames(data_final)

## refixing errors
data_final[!grepl("LO", rownames(data_final)) & rowSums(data[1:4])<20,]
data_final[!grepl("LO", rownames(data_final)) & rowSums(data[1:4])<20,1:4]<-c(NA,NA,NA,NA)

## save results
write.csv(data_final, file = "02_tidy-data/limnetic-ancestry-allsnps.csv",
          quote = F, row.names = F)

#### 12 BEST SNPS ####

## read in pop snp data
popsnp <- read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-final-snps.csv")

## choose best SNPs
LG_snps<-c("110","102","111","100","98","103","108","113","107","101","105","106")
SL_snps<-c("131","135","122","123","125","128","129","124","140","134_2","137","132_2")
SR_snps<-c("145","146","154","167","166_2","149","157","159","164","151","152","150")
WK_snps<-c("190","172","174","187","178","175","169","188","176","177","184","173")
best_snps<-c(LG_snps,SL_snps,SR_snps,WK_snps)
popsnp<-filter(popsnp, SNP.ID%in%best_snps)
## 48 snps with average freq = 88

## average freq by pop
popfreq<-popsnp %>% group_by(pop) %>% 
  summarise(pop_freq=mean(freq))
popsnp<-mutate(popsnp, pop_freq = case_when(pop == "LG" ~ popfreq[[1,2]],
                                            pop == "SL" ~ popfreq[[2,2]],
                                            pop == "SR" ~ popfreq[[3,2]],
                                            pop == "WK" ~ popfreq[[4,2]]))

## read in fish snp data
fish1<-read_csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-1.csv")
fish2<-read_csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-2.csv")
fish3<-read_csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-3.csv")
fish4<-read_csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-4.csv")
fish5<-read_csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-5.csv")
fish6<-read_csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-6.csv")

fish<-bind_rows(fish1,fish2,fish3,fish4,fish5,fish6)

fish <- fish %>%
  unique()
colnames(fish) <- str_replace_all(colnames(fish), "X", "")
colnames(fish)[1] = "SAMPLE.ID"

## filter lakes
fish$lake<-factor(substr(fish$SAMPLE.ID,1,2))
lakes<-c("CL","FL","HL","RL","LO")
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
# select(combined.mega, SNP.ID,freq, pop, sum)[order(combined.mega$sum),]
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
rowSums(data)
data[rowSums(data)==0,]

## fixing errors
boxplot(rowSums(filter(data,!grepl("LO", rownames(data)))))
boxplot.stats(rowSums(filter(data,!grepl("LO", rownames(data)))))$out
data[names(boxplot.stats(rowSums(filter(data,!grepl("LO", rownames(data)))))$out),]

data[!grepl("LO", rownames(data)) & rowSums(data)<10,]<-c(NA,NA,NA,NA)
summary(data)

## proportions
data$LG_prop<-data$LG/(data$LG+data$SL+data$SR+data$WK)
data$SL_prop<-data$SL/(data$LG+data$SL+data$SR+data$WK)
data$SR_prop<-data$SR/(data$LG+data$SL+data$SR+data$WK)
data$WK_prop<-data$WK/(data$LG+data$SL+data$SR+data$WK)

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
data_final<-cbind(data_final, data_raw[,1:4], fish2[,"LG_N":"WK_N"])
data_final$fish<-rownames(data_final)
data_final <- data_final[, c(9, 1:4, 5:8)]
colnames(data_final)[6:9]<-c("LG_raw","SL_raw","SR_raw","WK_raw")

## save results
write.csv(data_final, file = "02_tidy-data/limnetic-ancestry-bestsnps.csv",
          quote = F, row.names = F)
