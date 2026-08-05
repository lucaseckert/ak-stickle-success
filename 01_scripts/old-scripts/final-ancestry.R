## COMBINING ANCESTRY DATA ##

## packages
library(tidyverse)

## loading benthic data
benthic_all<-read.csv("./02_tidy-data/benthic-ancestry.csv") %>% 
  mutate(lake=factor(substr(fish,1,2)))

## loading limnetic data
limnetic_all<-read.csv("./02_tidy-data/limnetic-ancestry.csv") %>% 
  mutate(lake=factor(substr(fish,1,2)))

## combining data
lakes<-c("CC","CL","LL","LP","RL","HL","FL")
data_all<-full_join(benthic_all,limnetic_all, by=c("fish","lake")) %>% 
  filter(lake %in% lakes) %>% 
  select(fish, lake, FG,TL,WB,WT,LG,SL,SR,WK)
summary(data_all)

## fixing loon
loonL_all<- limnetic_all %>% filter(lake == "LO")
loonB_all<- benthic_all %>% filter(lake == "LO")

loon_all<- full_join(loonB_all,loonL_all, by=c("fish","lake"))

loon_all<-loon_all %>% 
  mutate(FG=FG_raw*(0.125/(FG_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N))),
         TL=TL_raw*(0.125/(TL_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N))),
         WB=WB_raw*(0.125/(WB_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N))),
         WT=WT_raw*(0.125/(WT_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N))),
         LG=LG_raw*(0.125/(LG_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N))),
         SL=SL_raw*(0.125/(SL_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N))),
         SR=SR_raw*(0.125/(SR_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N))),
         WK=WK_raw*(0.125/(WK_N/(FG_N+TL_N+WB_N+WT_N+LG_N+SL_N+SR_N+WK_N)))) %>% 
  select(fish,lake,FG,TL,WB,WT,LG,SL,SR,WK)

boxplot(rowSums(loon_all[,3:10]))
boxplot.stats(rowSums(loon_all[,3:10]))$out
loon_all[rowSums(loon_all[,3:10])<16,3:10]<-c(NA,NA,NA,NA)
 
loon_all<-loon_all %>%  mutate(FG_prop=FG/(FG+TL+WB+WT+LG+SL+SR+WK),
         TL_prop=TL/(FG+TL+WB+WT+LG+SL+SR+WK),
         WB_prop=WB/(FG+TL+WB+WT+LG+SL+SR+WK),
         WT_prop=WT/(FG+TL+WB+WT+LG+SL+SR+WK),
         LG_prop=LG/(FG+TL+WB+WT+LG+SL+SR+WK),
         SL_prop=SL/(FG+TL+WB+WT+LG+SL+SR+WK),
         SR_prop=SR/(FG+TL+WB+WT+LG+SL+SR+WK),
         WK_prop=WK/(FG+TL+WB+WT+LG+SL+SR+WK)) %>% 
  select(fish,lake,FG_prop,TL_prop,WB_prop,WT_prop,LG_prop,SL_prop,SR_prop,WK_prop)

## rounding loon
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
loon_all[is.na(loon_all)] <- 0
loon_all_final<-as.data.frame(t(apply(loon_all[,3:10]*cats,1,sum.round)))/cats
colnames(loon_all_final)=c("FG","TL","WB","WT","LG","SL","SR","WK")
rowSums(loon_all_final)
loon_all_final[rowSums(loon_all_final)==0,]<-c(rep(NA,8))

## fixing loon
loon_all_final$fish<-loon_all$fish
loon_all_final$lake<-factor("LO")
loon_all_final<-select(loon_all_final,fish,lake,FG,TL,WB,WT,LG,SL,SR,WK)

## adding loon back
data_all <- bind_rows(data_all,loon_all_final)
summary(data_all)

## checking zero sums
data_all[rowSums(data_all[,3:10], na.rm = T)!=1,]

## writing csv
write.csv(data_all, file = "./02_tidy-data/final-ancestry.csv",
          quote = F, row.names = F)

