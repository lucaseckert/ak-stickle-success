#### STRCUTURE ANALYSIS ####

## packages
library(tidyverse)

#### DATA PREP ####

#### BENTHIC ####

## benthic trial data
ben_trial<-read.csv("./00_raw-data/01_snp-array-tests/00_benthic-test-data.csv")
colnames(ben_trial) <- str_replace_all(colnames(ben_trial), "X", "")
colnames(ben_trial)[1] = "sample"
ben_trial$sample<-str_replace_all(ben_trial$sample,"B","")

## adding pop data
pops<-read.csv("./02_tidy-data/extractions.csv")
FG<-as.character(filter(pops, pop=="FG")$id)
TL<-as.character(filter(pops, pop=="TL")$id)
WB<-as.character(filter(pops, pop=="WB")$id)
WT<-as.character(filter(pops, pop=="WT")$id)
ben_trial<-mutate(ben_trial, pop = factor(case_when(sample %in% FG ~ "FG",
                                             sample %in% TL ~ "TL",
                                             sample %in% WB ~ "WB",
                                             sample %in% WT ~ "WT"))) %>% 
  filter(is.na(pop)==F)

## getting rid of snps that didn't work and duplicates
badsnps<-c("10","13","30_2","56","75")
ben_trial<-select(ben_trial, -all_of(badsnps)) %>% 
  distinct(sample, .keep_all = T)

## benthic real data
ben_real1<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-1.csv")
ben_real2<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-2.csv")
ben_real3<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-3.csv")
ben_real<-bind_rows(ben_real1,ben_real2,ben_real3)
colnames(ben_real) <- str_replace_all(colnames(ben_real), "X", "")
ben_real$pop<-substr(ben_real$sample,1,2)

ben_real<-filter(ben_real, pop %in% c("CC","LL","LP"))

## benthic data
ben<-bind_rows(ben_trial,ben_real) %>% 
  select(-c("CallRate","Call.Rate"))
ben[ben=="NC"]<-NA
ben[ben=="NN"]<-NA

## split columns
split_columns <- function(df) {
  new_df <- df
  for (col_name in colnames(df)) {
    new_df <- new_df %>%
      mutate(
        !!paste0(col_name, "_a") := substr(.data[[col_name]], 1, 1),
        !!paste0(col_name, "_b") := substr(.data[[col_name]], 2, 2)
      ) %>%
      select(-{{col_name}})
  }
  return(new_df)
}

ben_split<-split_columns(select(ben, -c(sample,pop)))
ben_split$sample<-ben$sample
ben_split$pop<-ben$pop

## replace each base with number and NA with -9
ben_split[ben_split=="A"]<-1
ben_split[ben_split=="C"]<-2
ben_split[ben_split=="G"]<-3
ben_split[ben_split=="T"]<-4
ben_split[is.na(ben_split)]<--9

## pop flag
ben_split$PopFlag<-ifelse(ben_split$pop %in% c("FG","TL","WB","WT"),1,0)

## final
ben_final<-ben_split %>% mutate(pop = case_when(pop == "FG" ~ 1,
                                                pop == "TL" ~ 2,
                                                pop == "WB" ~ 3,
                                                pop == "WT" ~ 4,
                                                pop == "CC" ~ 5,
                                                pop == "LL" ~ 6,
                                                pop == "LP" ~ 7)) %>% 
  select(Label=sample,PopData=pop, PopFlag, "1_1":"75_2_2")

## data file
write.table(ben_final, file = "./01_scripts/structure/infile.txt", quote = F,
            row.names = F, col.names = F)

## snp names
snps<-t(data.frame(colnames(select(ben, -c(sample,pop)))))
write.table(snps, file = "./01_scripts/structure/snps.txt", quote = F,
            row.names = F, col.names = F)


#### LIMNETIC ####

## lim trial data
lim_trial<-read.csv("./00_raw-data/01_snp-array-tests/02_limnetic-test-data.csv")
colnames(lim_trial) <- str_replace_all(colnames(lim_trial), "X", "")
colnames(lim_trial)[1] = "sample"
lim_trial$sample<-str_replace_all(lim_trial$sample,"L","")

## getting pops
lim_trial<-lim_trial %>%  mutate(pop=sample_pop) %>% 
  filter(pop %in% c("LG","SL","SR","WK")) %>% 
  select(-c(sample_pop,CallRate))

## NAs
lim_trial[lim_trial=="NC"]<-NA
lim_trial[lim_trial=="NN"]<-NA

colSums(is.na(lim_trial))
badsnps<-c("109","118","158","165","166","182")
lim_trial<-select(lim_trial, -all_of(badsnps))

## limnetic real data
lim_real1<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-1.csv", stringsAsFactors = T)
lim_real2<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-2.csv", stringsAsFactors = T)
lim_real3<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-3.csv", stringsAsFactors = T)
lim_real4<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-4.csv", stringsAsFactors = T)
lim_real5<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-5.csv", stringsAsFactors = T)
lim_real6<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-6.csv", stringsAsFactors = T)
lim_real<-bind_rows(lim_real1,lim_real2,lim_real3,lim_real4,lim_real5,lim_real6)
colnames(lim_real) <- str_replace_all(colnames(lim_real), "X", "")
lim_real$pop<-substr(lim_real$sample,1,2)
lim_real$pop[lim_real$pop=="S1"]<-"SR"
lim_real<-filter(lim_real, pop %in% c("CL","FL","HL","RL","SR"))

## combine
lim<-bind_rows(lim_trial,lim_real) %>% 
  select(-CallRate)
lim[lim=="NC"]<-NA
lim[lim=="NN"]<-NA

## split
lim_split<-split_columns(select(lim, -c(sample,pop)))
lim_split$sample<-lim$sample
lim_split$pop<-lim$pop

## replace each base with number and NA with -9
lim_split[lim_split=="A"]<-1
lim_split[lim_split=="C"]<-2
lim_split[lim_split=="G"]<-3
lim_split[lim_split=="T"]<-4
lim_split[is.na(lim_split)]<--9

## pop flag
lim_split$PopFlag<-ifelse(lim_split$pop %in% c("LG","SL","SR","WK"),1,0)

## final
lim_final<-lim_split %>% mutate(pop = case_when(pop == "LG" ~ 1,
                                                pop == "SL" ~ 2,
                                                pop == "SR" ~ 3,
                                                pop == "WK" ~ 4,
                                                pop == "CL" ~ 5,
                                                pop == "FL" ~ 6,
                                                pop == "HL" ~ 7,
                                                pop == "RL" ~ 8)) %>% 
  select(Label=sample,PopData=pop, PopFlag, "97_a":"182_2_b")

## data file
write.table(lim_final, file = "./02_tidy-data/structure/lim-infile.txt", quote = F,
            row.names = F, col.names = F)

## snp names
snps<-t(data.frame(colnames(select(lim, -c(sample,pop)))))
write.table(snps, file = "./02_tidy-data/structure/lim-snps.txt", quote = F,
            row.names = F, col.names = F)

#### BOTH ####

## benthic trial data
ben_trial<-read.csv("./00_raw-data/01_snp-array-tests/00_benthic-test-data.csv")
colnames(ben_trial) <- str_replace_all(colnames(ben_trial), "X", "")
colnames(ben_trial)[1] = "sample"
ben_trial$sample<-str_replace_all(ben_trial$sample,"B","")

## adding pop data
pops<-read.csv("./02_tidy-data/extractions.csv")
FG<-as.character(filter(pops, pop=="FG")$id)
TL<-as.character(filter(pops, pop=="TL")$id)
WB<-as.character(filter(pops, pop=="WB")$id)
WT<-as.character(filter(pops, pop=="WT")$id)
LG<-as.character(filter(pops, pop=="LG")$id)
SL<-as.character(filter(pops, pop=="SL")$id)
WK<-as.character(filter(pops, pop=="WK")$id)
ben_trial<-mutate(ben_trial, pop = factor(case_when(sample %in% FG ~ "FG",
                                                    sample %in% TL ~ "TL",
                                                    sample %in% WB ~ "WB",
                                                    sample %in% WT ~ "WT",
                                                    sample %in% LG ~ "LG",
                                                    sample %in% SL ~ "SL",
                                                    sample %in% WK ~ "WK"))) %>% 
  filter(is.na(pop)==F)

## getting rid of snps that didn't work and duplicates
badsnps<-c("10","13","30_2","56","75")
ben_trial<-select(ben_trial, -all_of(badsnps)) %>% 
  distinct(sample, .keep_all = T) %>% 
  select(-CallRate)

## benthic real data
ben_real1<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-1.csv")
ben_real2<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-2.csv")
ben_real3<-read.csv("./00_raw-data/02_snp-array-data/00_benthic-array/benthic-data/benthic-data-3.csv")
ben_real<-bind_rows(ben_real1,ben_real2,ben_real3)
colnames(ben_real) <- str_replace_all(colnames(ben_real), "X", "")
ben_real$pop<-substr(ben_real$sample,1,2)
ben_real$pop[ben_real$pop=="S1"]<-"SR"
ben_real<-filter(ben_real, pop %in% c("SR","LO","GL")) %>% 
  select(-CallRate,-Call.Rate)

## lim trial data
lim_trial<-read.csv("./00_raw-data/01_snp-array-tests/02_limnetic-test-data.csv")
colnames(lim_trial) <- str_replace_all(colnames(lim_trial), "X", "")
colnames(lim_trial)[1] = "sample"
lim_trial$sample<-str_replace_all(lim_trial$sample,"L","")

## getting pops
lim_trial<-lim_trial %>%  mutate(pop = case_when(sample %in% FG ~ "FG",
                                                 sample %in% TL ~ "TL",
                                                 sample %in% WB ~ "WB",
                                                 sample %in% WT ~ "WT",
                                                 sample %in% LG ~ "LG",
                                                 sample %in% SL ~ "SL",
                                                 sample %in% WK ~ "WK")) %>% 
  filter(is.na(pop)==F) %>% 
  select(-c(sample_pop,CallRate))

## NAs
lim_trial[lim_trial=="NC"]<-NA
lim_trial[lim_trial=="NN"]<-NA

colSums(is.na(lim_trial))
badsnps<-c("109","118","158","165","166","182")
lim_trial<-select(lim_trial, -all_of(badsnps))

## limnetic real data
lim_real1<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-1.csv", stringsAsFactors = T)
lim_real2<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-2.csv", stringsAsFactors = T)
lim_real3<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-3.csv", stringsAsFactors = T)
lim_real4<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-4.csv", stringsAsFactors = T)
lim_real5<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-5.csv", stringsAsFactors = T)
lim_real6<-read.csv("./00_raw-data/02_snp-array-data/01_limnetic-array/limnetic-data/limnetic-data-6.csv", stringsAsFactors = T)
lim_real<-bind_rows(lim_real1,lim_real2,lim_real3,lim_real4,lim_real5,lim_real6)
colnames(lim_real) <- str_replace_all(colnames(lim_real), "X", "")
lim_real$pop<-substr(lim_real$sample,1,2)
lim_real$pop[lim_real$pop=="S1"]<-"SR"
lim_real<-filter(lim_real, pop %in% c("GL","LO","SR")) %>% 
  select(-CallRate)

## combine
trial<-full_join(lim_trial,ben_trial)
real<-full_join(lim_real,ben_real)
all<-bind_rows(trial,real)

## NA
all[all=="NC"]<-NA
all[all=="NN"]<-NA

## split
all_split<-split_columns(select(all, -c(sample,pop)))
all_split$sample<-all$sample
all_split$pop<-all$pop

## replace each base with number and NA with -9
all_split[all_split=="A"]<-1
all_split[all_split=="C"]<-2
all_split[all_split=="G"]<-3
all_split[all_split=="T"]<-4
all_split[is.na(all_split)]<--9

## pop flag
all_split$PopFlag<-ifelse(all_split$pop %in% c("LG","SL","SR","WK","FG","TL","WB","WT"),1,0)

## final
all_final<-all_split %>% mutate(pop = case_when(pop == "FG" ~ 1,
                                                pop == "TL" ~ 2,
                                                pop == "WB" ~ 3,
                                                pop == "WT" ~ 4,
                                                pop == "LG" ~ 5,
                                                pop == "SL" ~ 6,
                                                pop == "SR" ~ 7,
                                                pop == "WK" ~ 8,
                                                pop == "LO" ~ 9,
                                                pop == "GL" ~ 10)) %>% 
  select(Label=sample,PopData=pop, PopFlag, "97_a":"75_2_b")

## data file
write.table(all_final, file = "./02_tidy-data/structure/all-infile.txt", quote = F,
            row.names = F, col.names = F)

## snp names
snps<-t(data.frame(colnames(select(all, -c(sample,pop)))))
write.table(snps, file = "./02_tidy-data/structure/all-snps.txt", quote = F,
            row.names = F, col.names = F)

#### RUNNING STRUCTURE ####

## directory
## dir = "./Library/CloudStorage/OneDrive-McGillUniversity/PHDTHESIS/R/ak-stickle-ancestry/02_tidy-data/structure"
## cd "$dir"

## benthic commands
## /Users/lucaseckert/Library/CloudStorage/OneDrive-McGillUniversity/PHD\ THESIS/R/ak-stickle-ancestry/02_tidy-data/structure/structure  -m /Users/lucaseckert/Library/CloudStorage/OneDrive-McGillUniversity/PHD\ THESIS/R/ak-stickle-ancestry/02_tidy-data/structure/ben-mainparams

## limnetic commands
## /Users/lucaseckert/Library/CloudStorage/OneDrive-McGillUniversity/PHD\ THESIS/R/ak-stickle-ancestry/02_tidy-data/structure/structure  -m /Users/lucaseckert/Library/CloudStorage/OneDrive-McGillUniversity/PHD\ THESIS/R/ak-stickle-ancestry/02_tidy-data/structure/lim-mainparams

## all commands
## /Users/lucaseckert/Library/CloudStorage/OneDrive-McGillUniversity/PHD\ THESIS/R/ak-stickle-ancestry/02_tidy-data/structure/structure  -m /Users/lucaseckert/Library/CloudStorage/OneDrive-McGillUniversity/PHD\ THESIS/R/ak-stickle-ancestry/02_tidy-data/structure/all-mainparams

#### PARSING RESULTS ####

## benthic results
ben_text<-readLines("./02_tidy-data/structure/ben-outfile_f")[140:428]
ben_text<-substring(ben_text, first = 5)
ben_list <- lapply(ben_text, function(line) {
  unlist(strsplit(line, "\\s+"))
})
ben_results <- as.data.frame(do.call(rbind, ben_list))[,c(1,5:8)]
colnames(ben_results)<-c("fish","FG","TL","WB","WT")
ben_results[,2:5] <- apply(ben_results[,2:5], 2, as.numeric)


## limnetic results
lim_text<-readLines("./02_tidy-data/structure/lim-outfile_f")[c(123:199,218:526)]
lim_text<-substring(lim_text, first = 5)
lim_list <- lapply(lim_text, function(line) {
  unlist(strsplit(line, "\\s+"))
})
lim_results <- as.data.frame(do.call(rbind, lim_list))[,c(1,5:8)]
colnames(lim_results)<-c("fish","LG","SL","SR","WK")
lim_results[,2:5] <- apply(lim_results[,2:5], 2, as.numeric)

## all results
all_text<-readLines("./02_tidy-data/structure/all-outfile_f")[223:377]
all_text<-c(substring(all_text[c(1:95,154,155)], first = 5),substring(all_text[96:153], first = 6))
all_list <- lapply(all_text, function(line) {
  unlist(strsplit(line, "\\s+"))
})
all_results <- as.data.frame(do.call(rbind, all_list))[,c(1,5:12)]
colnames(all_results)<-c("fish","FG","TL","WB","WT","LG","SL","SR","WK")
all_results[,2:9] <- apply(all_results[,2:9], 2, as.numeric)

## max gens
maxGen<-2
cats<-2^maxGen

## function for rounding while preserving the sum
sum.round <- function(x, digits = 0) {
  x<-as.numeric(x)
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

## rounding benthic
ben_final<-as.data.frame(t(apply(ben_results[,2:5]*cats,1,sum.round)))/cats
colnames(ben_final)=colnames(ben_results[,2:5])
ben_final$fish<-ben_results$fish

## rounding limnetic
lim_final<-as.data.frame(t(apply(lim_results[,2:5]*cats,1,sum.round)))/cats
colnames(lim_final)=colnames(lim_results[,2:5])
lim_final$fish<-lim_results$fish

## rounding all
all_final<-as.data.frame(t(apply(all_results[,2:9]*cats,1,sum.round)))/cats
colnames(all_final)=colnames(all_results[,2:9])
all_final$fish<-all_results$fish

#### COMPARING ####

## data
mine<-read.csv("./02_tidy-data/final-ancestry-allsnps.csv", stringsAsFactors = T)
ben_mine<-mine %>% filter(lake %in% c("CC","LL","LP")) %>% 
  select(fish,FG,TL,WB,WT)
lim_mine<-mine %>% filter(lake %in% c("CL","FL","HL","RL")) %>% 
  select(fish,LG,SL,SR,WK)
all_mine<-mine %>% filter(lake == "LO") %>% 
  select(fish,FG,TL,WB,WT,LG,SL,SR,WK)

## benthic comp
ben_comp<-ben_final[order(ben_final$fish),1:4] == ben_mine[order(ben_mine$fish),2:5]
sum(rowSums(ben_comp)!=4, na.rm = T)/nrow(ben_final)

## limnetic comp
lim_comp<-lim_final[order(lim_final$fish),1:4] == lim_mine[order(lim_mine$fish),2:5]
sum(rowSums(lim_comp)!=4, na.rm = T)/nrow(lim_final)

## loon comp
loon<-mutate(all_final, lake=substr(fish,1,2)) %>% filter(lake=="LO")
all_comp<-loon[order(loon$fish),1:8] == all_mine[order(all_mine$fish),2:9]
sum(rowSums(all_comp)!=8, na.rm = T)/nrow(all_comp)

sum(rowSums(loon==0.25)>0)/nrow(all_comp)
sum(rowSums(all_mine==0.25)>0, na.rm = T)/nrow(all_comp)


