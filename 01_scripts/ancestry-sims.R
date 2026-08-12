#### ANCESTRY SIMULATIONS #### 

## Description: Simulations to estimate the accuracy of ancestry inference from the chosen SNPs
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## packages
library(tidyverse)

## general parameters
N_loci <- 100
Popsize <- 1000
Nfish_per_source <- 250
N_pops_mixed <- 4

## function
anc_inf_sim <- function(N, freq, gen){
  
  ## genomes
  popA_fun<-function(freq){
    matrix(c(sample(c(rep("A",2*N_loci*freq[1]),rep("a",(2*N_loci)-(2*N_loci*freq[1])))),rep("X",600)),nrow = 2)
  }
  popA<-c(replicate(Nfish_per_source, popA_fun(freq)))
  
  popB_fun<-function(freq){
    matrix(c(rep("X",200),sample(c(rep("B",2*N_loci*freq[2]),rep("b",(2*N_loci)-(2*N_loci*freq[2])))),rep("X",400)),nrow = 2)
  }
  popB<-c(replicate(Nfish_per_source, popB_fun(freq)))
  
  popC_fun<-function(freq){
    matrix(c(rep("X",400),sample(c(rep("C",2*N_loci*freq[3]),rep("c",(2*N_loci)-(2*N_loci*freq[3])))),rep("X",200)),nrow = 2)
  }
  popC<-c(replicate(Nfish_per_source, popC_fun(freq)))
  
  popD_fun<-function(freq){
    matrix(c(rep("X",600),sample(c(rep("D",2*N_loci*freq[4]),rep("d",(2*N_loci)-(2*N_loci*freq[4]))))),nrow = 2)
  }
  popD<-c(replicate(Nfish_per_source, popD_fun(freq)))
  
  popAll<-c(popA,popB,popC,popD)
  
  NewPopulation<-lapply(split(popAll, (seq_along(popAll) - 1) %/% 800), matrix, nrow = 2, ncol = 400)
  
  ## new generation
  reproduce <- function(ParentPopulation){
    FX_population <- list()
    for(i in 1:Popsize){
      dam <- ParentPopulation[sample(c(1:length(ParentPopulation)),1)][[1]]
      sire <- ParentPopulation[sample(c(1:length(ParentPopulation)),1)][[1]]
      egghaploid <- c()
      spermhaploid <- c()
      for(l in 1:(4*N_loci)){
        egghaploid <- c(egghaploid, dam[rbinom(1,1,0.5)+1,l])
        spermhaploid <- c(spermhaploid, sire[rbinom(1,1,0.5)+1,l])
      }
      FX_population[[i]] <- as.matrix(rbind(egghaploid,spermhaploid))
    }
    return(FX_population)
  }
  ParentPopulation <- NewPopulation
  
  ## offspring become population
  ParentPopulation <- NewPopulation
  for(g in 1:gen){
    Newgeneration <- reproduce(ParentPopulation)
    ParentPopulation <- Newgeneration
  }
  
  ## actual ancestry
  Actual_ancestry <- data.frame(matrix(ncol = N_pops_mixed, nrow = Popsize))
  for(i in 1:Popsize){
    individualgeno <- Newgeneration[[i]]
    Actual_ancestry[i,1] <- 4*sum(individualgeno == "a"|individualgeno == "A")/(length(individualgeno))
    Actual_ancestry[i,2] <- 4*sum(individualgeno == "b"|individualgeno == "B")/(length(individualgeno))
    Actual_ancestry[i,3] <- 4*sum(individualgeno == "c"|individualgeno == "C")/(length(individualgeno))
    Actual_ancestry[i,4] <- 4*sum(individualgeno == "d"|individualgeno == "D")/(length(individualgeno))
  }
  colnames(Actual_ancestry) <- c("Prop_A","Prop_B", "Prop_C","Prop_D")
  
  ##  inferred ancestry
  Inferred_ancestry <- data.frame(matrix(ncol = N_pops_mixed, nrow = Popsize))
  
  focal_loci_A <- sample(c(1:100), size = N[1])
  focal_loci_B <- sample(c(101:200), size = N[2])
  focal_loci_C <- sample(c(201:300), size = N[3])
  focal_loci_D <- sample(c(301:400), size = N[4])
  N_focal_loci<-length(c(focal_loci_A,focal_loci_B,focal_loci_C,focal_loci_D))
  
  for(i in 1:Popsize){
    individualgeno <- Newgeneration[[i]]
    Inferred_ancestry[i,1] <- sum(individualgeno[,focal_loci_A] == "A")*(1/freq[1])*((N_focal_loci/4)/length(focal_loci_A))
    Inferred_ancestry[i,2] <- sum(individualgeno[,focal_loci_B] == "B")*(1/freq[2])*((N_focal_loci/4)/length(focal_loci_B))
    Inferred_ancestry[i,3] <- sum(individualgeno[,focal_loci_C] == "C")*(1/freq[3])*((N_focal_loci/4)/length(focal_loci_C))
    Inferred_ancestry[i,4] <- sum(individualgeno[,focal_loci_D] == "D")*(1/freq[4])*((N_focal_loci/4)/length(focal_loci_D))
  }
  colnames(Inferred_ancestry) <- c("Score_A","Score_B", "Score_C","Score_D")
  
  Inferred_ancestry<-mutate(Inferred_ancestry, Prop_A = Score_A/(Score_A+Score_B+Score_C+Score_D),
                            Prop_B = Score_B/(Score_A+Score_B+Score_C+Score_D),
                            Prop_C = Score_C/(Score_A+Score_B+Score_C+Score_D),
                            Prop_D = Score_D/(Score_A+Score_B+Score_C+Score_D)) %>% 
    select(Prop_A,Prop_B,Prop_C,Prop_D)
  
  
  ## r squared
  actual<-Actual_ancestry %>% pivot_longer(everything()) %>% select(name, actual=value)
  inferred<-Inferred_ancestry %>% pivot_longer(everything()) %>% select(name, inferred=value)
  comp<-data.frame(actual=actual$actual, inferred=inferred$inferred)
  lm<-lm(actual~inferred, data = comp)
  sum<-summary(lm)
  r_squared<-sum$adj.r.squared
  
  ## rounding
  maxGen<- gen
  cats<-2^maxGen

  sum.round <- function(x, digits = 0) {
    x<-as.numeric(x)
    up <- 10 ^ digits
    x <- x * up
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y / up
  }

  actual_ancestry_round<-as.data.frame(t(apply(Actual_ancestry*cats,1,sum.round)))/cats
  colnames(actual_ancestry_round)<-colnames(Actual_ancestry)
  
  inferred_ancestry_round<-as.data.frame(t(apply(Inferred_ancestry*cats,1,sum.round)))/cats
  colnames(inferred_ancestry_round)<-colnames(Inferred_ancestry)
  
  ## how many are wrong
  Nwrong_bins <- rowSums(1*as.matrix(actual_ancestry_round != inferred_ancestry_round))
  wrong<-sum(Nwrong_bins > 0)/nrow(actual_ancestry_round)*100
  
  Freq_dropout <-  rowSums(1*as.matrix((actual_ancestry_round > 0) & (inferred_ancestry_round == 0)))
  dropped<-(sum(Freq_dropout>0) / length(Freq_dropout))*100
  
  ## output
  print(c(wrong,dropped,r_squared))
}

## limnetic all snp parameters
N_loci_A_l<-20
freq_A_l<-0.70
N_loci_B_l<-17
freq_B_l<-0.90
N_loci_C_l<-24
freq_C_l<-0.90
N_loci_D_l<-23
freq_D_l<-0.81

## benthic all snp parameters
N_loci_A_b<-23
freq_A_b<-0.70
N_loci_B_b<-20
freq_B_b<-1
N_loci_C_b<-18
freq_C_b<-0.75
N_loci_D_b<-13
freq_D_b<-0.73

## inputs
N_ben<-c(N_loci_A_b, N_loci_B_b, N_loci_C_b, N_loci_D_b)
freq_ben<-c(freq_A_b, freq_B_b, freq_C_b, freq_D_b)
N_lim<-c(N_loci_A_l, N_loci_B_l, N_loci_C_l, N_loci_D_l)
freq_lim<-c(freq_A_l, freq_B_l, freq_C_l, freq_D_l)

## ci function
ci<-function(sample){
  mean(sample) + qt(0.975,df=length(sample)-1)*sd(sample)/sqrt(length(sample))
}

## limnetic F1
## lim_F1<-replicate(100,anc_inf_sim(N_lim,freq_lim,1))
## saveRDS(lim_F1, file = "02_tidy-data/ancestry-sims/lim_F1.RDS")
lim_F1<-readRDS("02_tidy-data/ancestry-sims/lim_F1.RDS")

lim_F1_wrong_mean<-mean(lim_F1[1,])
lim_F1_wrong_ci <-ci(lim_F1[1,])
lim_F1_miss_mean<-mean(lim_F1[2,])
lim_F1_miss_ci<-ci(lim_F1[2,])

## limnetic F2
lim_F2<-replicate(100,anc_inf_sim(N_lim,freq_lim,2))
saveRDS(lim_F2, file = "02_tidy-data/ancestry-sims/lim_F2.RDS")
lim_F2<-readRDS("02_tidy-data/ancestry-sims/lim_F2.RDS")

lim_F2_wrong_mean<-mean(lim_F2[1,])
lim_F2_wrong_ci <-ci(lim_F2[1,])
lim_F2_miss_mean<-mean(lim_F2[2,])
lim_F2_miss_ci<-ci(lim_F2[2,])

hist(lim_F2[1,])
hist(lim_F2[2,])

## limnetic F3
lim_F3<-replicate(100,anc_inf_sim(N_lim,freq_lim,3))
saveRDS(lim_F3, file = "02_tidy-data/ancestry-sims/lim_F3.RDS")
lim_F3<-readRDS("02_tidy-data/ancestry-sims/lim_F3.RDS")

lim_F3_wrong_mean<-mean(lim_F3[1,])
lim_F3_wrong_ci <-ci(lim_F3[1,])
lim_F3_miss_mean<-mean(lim_F3[2,])
lim_F3_miss_ci<-ci(lim_F3[2,])

## benthic F1
## ben_F1<-replicate(100,anc_inf_sim(N_ben,freq_ben,1))
## saveRDS(ben_F1, file = "02_tidy-data/ancestry-sims/ben_F1.RDS")
ben_F1<-readRDS("02_tidy-data/ancestry-sims/ben_F1.RDS")

ben_F1_wrong_mean<-mean(ben_F1[1,])
ben_F1_wrong_ci <-ci(ben_F1[1,])
ben_F1_miss_mean<-mean(ben_F1[2,])
ben_F1_miss_ci<-ci(ben_F1[2,])

## benthic F2
## ben_F2<-replicate(100,anc_inf_sim(N_ben,freq_ben,2))
## saveRDS(ben_F2, file = "02_tidy-data/ancestry-sims/ben_F2.RDS")
ben_F2<-readRDS("02_tidy-data/ancestry-sims/ben_F2.RDS")

ben_F2_wrong_mean<-mean(ben_F2[1,])
ben_F2_wrong_ci <-ci(ben_F2[1,])
ben_F2_miss_mean<-mean(ben_F2[2,])
ben_F2_miss_ci<-ci(ben_F2[2,])

## benthic F3
## ben_F3<-replicate(100,anc_inf_sim(N_ben,freq_ben,3))
## saveRDS(ben_F3, file = "02_tidy-data/ancestry-sims/ben_F3.RDS")
ben_F3<-readRDS("02_tidy-data/ancestry-sims/ben_F3.RDS")

ben_F3_wrong_mean<-mean(ben_F3[1,])
ben_F3_wrong_ci <-ci(ben_F3[1,])
ben_F3_miss_mean<-mean(ben_F3[2,])
ben_F3_miss_ci<-ci(ben_F3[2,])

## F1 summary
F1_wrong_mean<-mean(c(lim_F1_wrong_mean,ben_F1_wrong_mean))
F1_wrong_ci<-mean(c(lim_F1_wrong_ci,ben_F1_wrong_ci))
F1_miss_mean<-mean(c(lim_F1_miss_mean, ben_F1_miss_mean))
F1_miss_ci<-mean(c(lim_F1_miss_ci, ben_F1_miss_ci))

## F2 summary
F2_wrong_mean<-mean(c(lim_F2_wrong_mean,ben_F2_wrong_mean))
F2_wrong_ci<-mean(c(lim_F2_wrong_ci,ben_F2_wrong_ci))
F2_miss_mean<-mean(c(lim_F2_miss_mean, ben_F2_miss_mean))
F2_miss_ci<-mean(c(lim_F2_miss_ci, ben_F2_miss_ci))

## F3 summary
F3_wrong_mean<-mean(c(lim_F3_wrong_mean,ben_F3_wrong_mean))
F3_wrong_ci<-mean(c(lim_F3_wrong_ci,ben_F3_wrong_ci))
F3_miss_mean<-mean(c(lim_F3_miss_mean, ben_F3_miss_mean))
F3_miss_ci<-mean(c(lim_F3_miss_ci, ben_F3_miss_ci))

## comparison with using just the best snps
N_loci_A_l_best<-12
freq_A_l_best<-0.74
N_loci_B_l_best<-12
freq_B_l_best<-0.94
N_loci_C_l_best<-12
freq_C_l_best<-0.97
N_loci_D_l_best<-12
freq_D_l_best<-0.87

N_lim_best<-c(N_loci_A_l_best, N_loci_B_l_best, N_loci_C_l_best, N_loci_D_l_best)
freq_lim_best<-c(freq_A_l_best, freq_B_l_best, freq_C_l_best, freq_D_l_best)

## lim_F2_best<-replicate(100,anc_inf_sim(N_lim_best,freq_lim_best,2))
## saveRDS(lim_F2_best, file = "02_tidy-data/ancestry-sims/lim_F2_best.RDS")
lim_F2_best<-readRDS("02_tidy-data/ancestry-sims/lim_F2_best.RDS")

lim_F2_best_wrong_mean<-mean(lim_F2_best[1,])
lim_F2_best_wrong_ci <-ci(lim_F2_best[1,])
lim_F2_best_miss_mean<-mean(lim_F2_best[2,])
lim_F2_best_miss_ci<-ci(lim_F2_best[2,])

#### LOON ####

## loon parameters
N_loci <- 100
Popsize <- 1000
Nfish_per_source <- 125
N_pops_mixed <- 8

## loon input
N_all<-c(N_ben,N_lim)
freq_all<-c(freq_ben,freq_lim)

## loon function
anc_inf_sim_loon <- function(N, freq, gen){
  
  ## genomes
  popA_fun<-function(freq){
    matrix(c(sample(c(rep("A",2*N_loci*freq[1]),rep("a",(2*N_loci)-(2*N_loci*freq[1])))),rep("X",1400)),nrow = 2)
  }
  popA<-c(replicate(Nfish_per_source, popA_fun(freq)))
  
  popB_fun<-function(freq){
    matrix(c(rep("X",200),sample(c(rep("B",2*N_loci*freq[2]),rep("b",(2*N_loci)-(2*N_loci*freq[2])))),rep("X",1200)),nrow = 2)
  }
  popB<-c(replicate(Nfish_per_source, popB_fun(freq)))
  
  popC_fun<-function(freq){
    matrix(c(rep("X",400),sample(c(rep("C",2*N_loci*freq[3]),rep("c",(2*N_loci)-(2*N_loci*freq[3])))),rep("X",1000)),nrow = 2)
  }
  popC<-c(replicate(Nfish_per_source, popC_fun(freq)))
  
  popD_fun<-function(freq){
    matrix(c(rep("X",600),sample(c(rep("D",2*N_loci*freq[4]),rep("d",(2*N_loci)-(2*N_loci*freq[4])))),rep("X",800)),nrow = 2)
  }
  popD<-c(replicate(Nfish_per_source, popD_fun(freq)))
  
  popE_fun<-function(freq){
    matrix(c(rep("X",800),sample(c(rep("E",2*N_loci*freq[5]),rep("e",(2*N_loci)-(2*N_loci*freq[5])))),rep("X",600)),nrow = 2)
  }
  popE<-c(replicate(Nfish_per_source, popE_fun(freq)))
  
  popF_fun<-function(freq){
    matrix(c(rep("X",1000),sample(c(rep("F",2*N_loci*freq[6]),rep("f",(2*N_loci)-(2*N_loci*freq[6])))),rep("X",400)),nrow = 2)
  }
  popF<-c(replicate(Nfish_per_source, popF_fun(freq)))
  
  popG_fun<-function(freq){
    matrix(c(rep("X",1200),sample(c(rep("G",2*N_loci*freq[7]),rep("g",(2*N_loci)-(2*N_loci*freq[7])))),rep("X",200)),nrow = 2)
  }
  popG<-c(replicate(Nfish_per_source, popG_fun(freq)))
  
  popH_fun<-function(freq){
    matrix(c(rep("X",1400),sample(c(rep("H",2*N_loci*freq[8]),rep("h",(2*N_loci)-(2*N_loci*freq[8]))))),nrow = 2)
  }
  popH<-c(replicate(Nfish_per_source, popH_fun(freq)))
  
  popAll<-c(popA,popB,popC,popD,popE,popF,popG,popH)
  
  NewPopulation<-lapply(split(popAll, (seq_along(popAll) - 1) %/% 1600), matrix, nrow = 2, ncol = 800)
  
  ## new generation
  reproduce <- function(ParentPopulation){
    FX_population <- list()
    for(i in 1:Popsize){
      dam <- ParentPopulation[sample(c(1:length(ParentPopulation)),1)][[1]]
      sire <- ParentPopulation[sample(c(1:length(ParentPopulation)),1)][[1]]
      egghaploid <- c()
      spermhaploid <- c()
      for(l in 1:(8*N_loci)){
        egghaploid <- c(egghaploid, dam[rbinom(1,1,0.5)+1,l])
        spermhaploid <- c(spermhaploid, sire[rbinom(1,1,0.5)+1,l])
      }
      FX_population[[i]] <- as.matrix(rbind(egghaploid,spermhaploid))
    }
    return(FX_population)
  }
  ParentPopulation <- NewPopulation
  
  ## offspring become population
  ParentPopulation <- NewPopulation
  for(g in 1:gen){
    Newgeneration <- reproduce(ParentPopulation)
    ParentPopulation <- Newgeneration
  }
  
  ## actual ancestry
  Actual_ancestry <- data.frame(matrix(ncol = N_pops_mixed, nrow = Popsize))
  for(i in 1:Popsize){
    individualgeno <- Newgeneration[[i]]
    Actual_ancestry[i,1] <- 8*sum(individualgeno == "a"|individualgeno == "A")/(length(individualgeno))
    Actual_ancestry[i,2] <- 8*sum(individualgeno == "b"|individualgeno == "B")/(length(individualgeno))
    Actual_ancestry[i,3] <- 8*sum(individualgeno == "c"|individualgeno == "C")/(length(individualgeno))
    Actual_ancestry[i,4] <- 8*sum(individualgeno == "d"|individualgeno == "D")/(length(individualgeno))
    Actual_ancestry[i,5] <- 8*sum(individualgeno == "e"|individualgeno == "E")/(length(individualgeno))
    Actual_ancestry[i,6] <- 8*sum(individualgeno == "f"|individualgeno == "F")/(length(individualgeno))
    Actual_ancestry[i,7] <- 8*sum(individualgeno == "g"|individualgeno == "G")/(length(individualgeno))
    Actual_ancestry[i,8] <- 8*sum(individualgeno == "h"|individualgeno == "H")/(length(individualgeno))
  }
  colnames(Actual_ancestry) <- c("Prop_A","Prop_B", "Prop_C","Prop_D","Prop_E","Prop_F","Prop_G","Prop_H")
  
  ##  inferred ancestry
  Inferred_ancestry <- data.frame(matrix(ncol = N_pops_mixed, nrow = Popsize))
  
  focal_loci_A <- sample(c(1:100), size = N[1])
  focal_loci_B <- sample(c(101:200), size = N[2])
  focal_loci_C <- sample(c(201:300), size = N[3])
  focal_loci_D <- sample(c(301:400), size = N[4])
  focal_loci_E <- sample(c(401:500), size = N[1])
  focal_loci_F <- sample(c(501:600), size = N[2])
  focal_loci_G <- sample(c(601:700), size = N[3])
  focal_loci_H <- sample(c(701:800), size = N[4])
  N_focal_loci<-length(c(focal_loci_A,focal_loci_B,focal_loci_C,focal_loci_D,focal_loci_E,focal_loci_F,focal_loci_G,focal_loci_H))
  
  for(i in 1:Popsize){
    individualgeno <- Newgeneration[[i]]
    Inferred_ancestry[i,1] <- sum(individualgeno[,focal_loci_A] == "A")*(1/freq[1])*((N_focal_loci/8)/length(focal_loci_A))
    Inferred_ancestry[i,2] <- sum(individualgeno[,focal_loci_B] == "B")*(1/freq[2])*((N_focal_loci/8)/length(focal_loci_B))
    Inferred_ancestry[i,3] <- sum(individualgeno[,focal_loci_C] == "C")*(1/freq[3])*((N_focal_loci/8)/length(focal_loci_C))
    Inferred_ancestry[i,4] <- sum(individualgeno[,focal_loci_D] == "D")*(1/freq[4])*((N_focal_loci/8)/length(focal_loci_D))
    Inferred_ancestry[i,5] <- sum(individualgeno[,focal_loci_E] == "E")*(1/freq[5])*((N_focal_loci/8)/length(focal_loci_E))
    Inferred_ancestry[i,6] <- sum(individualgeno[,focal_loci_F] == "F")*(1/freq[6])*((N_focal_loci/8)/length(focal_loci_F))
    Inferred_ancestry[i,7] <- sum(individualgeno[,focal_loci_G] == "G")*(1/freq[7])*((N_focal_loci/8)/length(focal_loci_G))
    Inferred_ancestry[i,8] <- sum(individualgeno[,focal_loci_H] == "H")*(1/freq[8])*((N_focal_loci/8)/length(focal_loci_H))
  }
  colnames(Inferred_ancestry) <- c("Score_A","Score_B", "Score_C","Score_D","Score_E","Score_F", "Score_G","Score_H")
  
  Inferred_ancestry<-mutate(Inferred_ancestry, Prop_A = Score_A/(Score_A+Score_B+Score_C+Score_D+Score_E+Score_F+Score_G+Score_H),
                            Prop_B = Score_B/(Score_A+Score_B+Score_C+Score_D+Score_E+Score_F+Score_G+Score_H),
                            Prop_C = Score_C/(Score_A+Score_B+Score_C+Score_D+Score_E+Score_F+Score_G+Score_H),
                            Prop_D = Score_D/(Score_A+Score_B+Score_C+Score_D+Score_E+Score_F+Score_G+Score_H),
                            Prop_E = Score_E/(Score_A+Score_B+Score_C+Score_D+Score_E+Score_F+Score_G+Score_H),
                            Prop_F = Score_F/(Score_A+Score_B+Score_C+Score_D+Score_E+Score_F+Score_G+Score_H),
                            Prop_G = Score_G/(Score_A+Score_B+Score_C+Score_D+Score_E+Score_F+Score_G+Score_H),
                            Prop_H = Score_H/(Score_A+Score_B+Score_C+Score_D+Score_E+Score_F+Score_G+Score_H)) %>% 
    select(Prop_A,Prop_B,Prop_C,Prop_D,Prop_E,Prop_F,Prop_G,Prop_H)
  
  
  ## r squared
  actual<-Actual_ancestry %>% pivot_longer(everything()) %>% select(name, actual=value)
  inferred<-Inferred_ancestry %>% pivot_longer(everything()) %>% select(name, inferred=value)
  comp<-data.frame(actual=actual$actual, inferred=inferred$inferred)
  lm<-lm(inferred~actual, data = comp)
  sum<-summary(lm)
  r_squared<-sum$adj.r.squared
  
  ## missed
  comp<-comp %>% mutate(missed = inferred==0 & actual>0)
  missed_df<-as.data.frame(matrix(comp$missed, ncol=4, byrow = T))
  missed<-sum(rowSums(missed_df))/nrow(missed_df) 
  
  ## plot
  # ggplot(comp, aes(x=actual, y=inferred))+
  #   geom_point(color="grey25", alpha=0.5)+
  #   geom_smooth(method="lm", color="black")+
  #   xlim(0,1)+
  #   ylim(0,1)+
  #   theme_bw()
  
  ## shared
  shared_actual<-data.frame(mapply(pmin,Actual_ancestry[1:500,],Actual_ancestry[501:1000,])) %>% rowSums()
  hist(shared_actual)
  mean(shared_actual)
  shared_inferred<-data.frame(mapply(pmin,Inferred_ancestry[1:500,],Inferred_ancestry[501:1000,])) %>% rowSums()
 # hist(shared_inferred)
  mean(shared_inferred)
  
  shared_df<-data.frame(actual=shared_actual, inferred=shared_inferred) %>% pivot_longer(cols = everything())
  #ggplot(shared_df, aes(x=name, y=value))+
  #  geom_violin()
  
  ## rounding
  maxGen<- gen
  cats<-2^maxGen
  
  sum.round <- function(x, digits = 0) {
    x<-as.numeric(x)
    up <- 10 ^ digits
    x <- x * up
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y / up
  }
  
  actual_ancestry_round<-as.data.frame(t(apply(Actual_ancestry*cats,1,sum.round)))/cats
  colnames(actual_ancestry_round)<-colnames(Actual_ancestry)
  
  inferred_ancestry_round<-as.data.frame(t(apply(Inferred_ancestry*cats,1,sum.round)))/cats
  colnames(inferred_ancestry_round)<-colnames(Inferred_ancestry)
  
  ## how many are wrong
  Nwrong_bins <- rowSums(1*as.matrix(actual_ancestry_round != inferred_ancestry_round))
  wrong<-sum(Nwrong_bins > 0)/nrow(actual_ancestry_round)*100
  
  Freq_dropout <-  rowSums(1*as.matrix((actual_ancestry_round > 0) & (inferred_ancestry_round == 0)))
  dropped<-(sum(Freq_dropout>0) / length(Freq_dropout))*100
  
  ## output
  print(c(wrong,dropped,r_squared,missed))
}

## loon F1
loon_F1<-replicate(100,anc_inf_sim_loon(N_all,freq_all,1))
saveRDS(loon_F1, file = "02_tidy-data/ancestry-sims/loon_F1.RDS")
loon_F1<-readRDS("02_tidy-data/ancestry-sims/loon_F1.RDS")

loon_F1_wrong_mean<-mean(loon_F1[1,])
loon_F1_wrong_ci <-ci(loon_F1[1,])
loon_F1_miss_mean<-mean(loon_F1[2,])
loon_F1_miss_ci<-ci(loon_F1[2,])

hist(loon_F1[1,])
hist(loon_F1[2,])

## loon F2
loon_F2<-replicate(100,anc_inf_sim_loon(N_all,freq_all,2))
saveRDS(loon_F2, file = "02_tidy-data/ancestry-sims/loon_F2.RDS")
loon_F2<-readRDS("02_tidy-data/ancestry-sims/loon_F2.RDS")

loon_F2_r_mean<-mean(loon_F2[3,])
loon_F2_miss_mean<-mean(loon_F2[4,])
hist(loon_F2[1,])


## loon F3
loon_F3<-replicate(100,anc_inf_sim_loon(N_all,freq_all,3))
saveRDS(loon_F3, file = "02_tidy-data/ancestry-sims/loon_F3.RDS")
loon_F3<-readRDS("02_tidy-data/ancestry-sims/loon_F3.RDS")

loon_F3_r_mean<-mean(loon_F3[3,])
hist(loon_F3[3,])
loon_F3_miss_mean<-mean(loon_F3[4,])


