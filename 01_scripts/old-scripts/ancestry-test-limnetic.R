#### Ancestry Inference Tests ####

## DESCRIPTION: Inferring proportional ancestry from genotyping results from the small array
## AUTHORS: Lucas Eckert w/ help from Nikki + Isaac 

## load packages
library(tidyverse)

## read in pop snp data
popsnp <- read.csv("./00_raw-data/01_snp-array-tests/03_limnetic-test-snps.csv") %>% 
  filter(pop!="SR")

## read in fish snp data
fish<-read_csv("./00_raw-data/01_snp-array-tests/02_limnetic-test-data.csv") %>% 
  filter(sample_pop!="CO")

fish <- fish %>%
  unique()
colnames(fish) <- str_replace_all(colnames(fish), "X", "")
colnames(fish)[1] = "SAMPLE.ID"

## replace NCs
fish[fish == "NC"] <- "NN"

## make sure all fish snps are in population snp data 
fish_snps <- str_replace_all(colnames(fish), "X", "")
unique(popsnp$SNP.ID)
which(!fish_snps %in% unique(popsnp$SNP.ID))

## subset to columns we care about 
popsnp <- select(popsnp, SNP.ID, pop, Base_A, freq)

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
    
    if (mega$Base_A[n]==data[n]){mega[n,(4+i)]<-1}
    if (identical(mega$Base_A[n],data[n])==FALSE) {mega[n,(4+i)]<-0}
  }
}

## to weight by frequency /mega$freq[n]

## sum columns for each sample
combined.mega<-mega[,1:4]
for (i in 1:length(unique(gsub("_.char","",names(mega)[5:ncol(mega)])))){
  
  col=unique(gsub("_.char","",names(mega)[5:ncol(mega)]))[i]
  combined.mega[,i+4]<-as.numeric(mega[,c(paste0(col,"_1char"))])+as.numeric(mega[,c(paste0(col,"_2char"))])
  names(combined.mega)[i+4]<-col
}

## exploring SNPs
combined.mega$sum<-rowSums(combined.mega[,5:ncol(combined.mega)])
select(combined.mega, freq, pop, sum)
ggplot(combined.mega, aes(x=freq, y=sum, label=SNP.ID))+
  geom_point()+
  geom_text(hjust=-0.25, vjust=-0.25)+
  geom_abline(slope = 36, intercept = 0)+
  theme_bw()+
  facet_wrap(~pop)

filter(combined.mega, sum==0)$SNP.ID

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
colnames(data_final)=colnames(data[,5:8])
rowSums(data_final)