#### ANCESTRY INFERENCE SIMULATIONS #### 

## packages
library(tidyverse)

## general parameters
N_loci <- 100
generations <- 3
Popsize <- 1000
Nfish_per_source <- 250
N_pops_mixed <- 4

## limnetic all snp parameters
# N_loci_A<-20
# freq_A<-0.70
# N_loci_B<-17
# freq_B<-0.90
# N_loci_C<-24
# freq_C<-0.90
# N_loci_D<-23
# freq_D<-0.81

## benthic all snp parameters
N_loci_A<-23
freq_A<-0.70
N_loci_B<-20
freq_B<-1
N_loci_C<-18
freq_C<-0.75
N_loci_D<-13
freq_D<-0.73

## genomes for each population NEED TO FIX THIS
popA <- matrix(c(sample(c(rep("A",2*N_loci*freq_A),rep("a",(2*N_loci)-(2*N_loci*freq_A)))),rep("X",600)),nrow = 2)
popB <- matrix(c(rep("X",200),sample(c(rep("B",2*N_loci*freq_B),rep("b",(2*N_loci)-(2*N_loci*freq_B)))),rep("X",400)),nrow = 2)
popC <- matrix(c(rep("X",400),sample(c(rep("C",2*N_loci*freq_C),rep("c",(2*N_loci)-(2*N_loci*freq_C)))),rep("X",200)),nrow = 2)
popD <- matrix(c(rep("X",600),sample(c(rep("D",2*N_loci*freq_D),rep("d",(2*N_loci)-(2*N_loci*freq_D))))),nrow = 2)
Genotypes <- list(popA, popB, popC, popD)

## creating the mixed population
NewPopulation <- list()
counter <- 1
for(p in 1:N_pops_mixed)
{
  for(i in 1:Nfish_per_source){
    NewPopulation[counter] <- Genotypes[p]
    counter <- counter + 1
  }
}

## creating a new generation based on random mating of the parent generation
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

## the new generation becomes the next parent generation
ParentPopulation <- NewPopulation
for(g in 1:generations){
  Newgeneration <- reproduce(ParentPopulation)
  ParentPopulation <- Newgeneration
}

## determining the actual ancestry of the most recent generation
Actual_ancestry <- data.frame(matrix(ncol = N_pops_mixed, nrow = Popsize))
for(i in 1:Popsize){
  individualgeno <- Newgeneration[[i]]
  Actual_ancestry[i,1] <- 4*sum(individualgeno == "a"|individualgeno == "A")/(length(individualgeno))
  Actual_ancestry[i,2] <- 4*sum(individualgeno == "b"|individualgeno == "B")/(length(individualgeno))
  Actual_ancestry[i,3] <- 4*sum(individualgeno == "c"|individualgeno == "C")/(length(individualgeno))
  Actual_ancestry[i,4] <- 4*sum(individualgeno == "d"|individualgeno == "D")/(length(individualgeno))
}
colnames(Actual_ancestry) <- c("Prop_A","Prop_B", "Prop_C","Prop_D")

## estimating the inferred ancestry according to our focal loci 
Inferred_ancestry <- data.frame(matrix(ncol = N_pops_mixed, nrow = Popsize))

focal_loci_A <- sample(c(1:100), size = N_loci_A)
focal_loci_B <- sample(c(101:200), size = N_loci_B)
focal_loci_C <- sample(c(201:300), size = N_loci_C)
focal_loci_D <- sample(c(301:400), size = N_loci_D)
N_focal_loci<-length(c(focal_loci_A,focal_loci_B,focal_loci_C,focal_loci_D))

for(i in 1:Popsize){
  individualgeno <- Newgeneration[[i]]
  Inferred_ancestry[i,1] <- sum(individualgeno[,focal_loci_A] == "A")*(1/freq_A)*((N_focal_loci/4)/length(focal_loci_A))
  Inferred_ancestry[i,2] <- sum(individualgeno[,focal_loci_B] == "B")*(1/freq_B)*((N_focal_loci/4)/length(focal_loci_B))
  Inferred_ancestry[i,3] <- sum(individualgeno[,focal_loci_C] == "C")*(1/freq_C)*((N_focal_loci/4)/length(focal_loci_C))
  Inferred_ancestry[i,4] <- sum(individualgeno[,focal_loci_D] == "D")*(1/freq_D)*((N_focal_loci/4)/length(focal_loci_D))
}
colnames(Inferred_ancestry) <- c("Score_A","Score_B", "Score_C","Score_D")

Inferred_ancestry<-mutate(Inferred_ancestry, Prop_A = Score_A/(Score_A+Score_B+Score_C+Score_D),
                          Prop_B = Score_B/(Score_A+Score_B+Score_C+Score_D),
                          Prop_C = Score_C/(Score_A+Score_B+Score_C+Score_D),
                          Prop_D = Score_D/(Score_A+Score_B+Score_C+Score_D)) %>% 
  select(Prop_A,Prop_B,Prop_C,Prop_D)

# correlation of actual and inferred ancestry
corresult <- cor.test(Inferred_ancestry$Prop_A , Actual_ancestry$Prop_A)
plot(Inferred_ancestry$Prop_A ~ Actual_ancestry$Prop_A, main = paste("Correlation r = " , round(corresult$estimate, 5), " in Gen ", generations))

## fractions based on max generations
maxGen<- generations
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

## rounding both ancestry
actual_ancestry_round<-as.data.frame(t(apply(Actual_ancestry*cats,1,sum.round)))/cats
colnames(actual_ancestry_round)<-colnames(Actual_ancestry)

inferred_ancestry_round<-as.data.frame(t(apply(Inferred_ancestry*cats,1,sum.round)))/cats
colnames(inferred_ancestry_round)<-colnames(Inferred_ancestry)

## determining the number that are wrong
table(inferred_ancestry_round$Prop_A, actual_ancestry_round$Prop_A)
Nwrong_bins <- rowSums(1*as.matrix(actual_ancestry_round != inferred_ancestry_round))
sum(Nwrong_bins > 0)/nrow(actual_ancestry_round)*100

## number where we missed a pop
Freq_dropout <-  rowSums(1*as.matrix((actual_ancestry_round > 0) & (inferred_ancestry_round == 0)))
(sum(Freq_dropout>0) / length(Freq_dropout))*100

lim_F1_wrong<-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
lim_F1_drop<-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

ben_F1_wrong<-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
ben_F1_drop<-c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

lim_F2_wrong<-c(4.9, 3.8, 4.0, 3.8, 4.3, 5.2, 6.4, 4.5, 5.3, 5.1)
lim_F2_drop<-c(2.8, 2.0, 2.2, 2.0, 2.2, 2.4, 3.2, 2.2, 3.1, 2.4)
mean(lim_F2_wrong)
mean(lim_F2_drop)

ben_F2_wrong<-c(8.8, 8.8, 6.4, 9.3, 8.2, 9.0, 12.8, 11.8, 6.3, 7.9)
ben_F2_drop<-c(4.8, 5.5, 3.8, 3.9, 3.4, 4.7, 5.7, 6.2, 3.7, 3.5)
mean(ben_F2_wrong)
mean(ben_F2_drop)


