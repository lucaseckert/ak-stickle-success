##### STRUCTURE Comparison ####

## packages
library(tidyverse)

## data
mine<-read.csv("./02_tidy-data/final-ancestry-allsnps.csv", stringsAsFactors = T)
struc<-read.csv(file.choose(), stringsAsFactors = T)

## filter mine
mine<-mine %>% filter(lake %in% c("CL","FL","HL","RL")) %>% 
  select(fish,LG,SL,SR,WK)

## round mine
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
struc_final<-as.data.frame(t(apply(struc[,2:5]*cats,1,sum.round)))/cats
colnames(struc_final)=colnames(struc[,2:5])
struc_final$fish<-struc$fish

struc<-struc_final %>% filter(fish %in% mine$fish) %>% 
  select(fish,LG,SL,SR,WK)
comp<-struc[order(struc$fish),2:5] == mine[order(mine$fish),2:5]
view(comp)

sum(rowSums(comp)!=4, na.rm = T)
nrow(comp)


