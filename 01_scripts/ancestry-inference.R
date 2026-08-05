#### ANCESTRY INFERENCE ####

## Description: Final script for inferring source ancestry from SNP genotyping data
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## packages 
library(tidyverse)

## benthic snps
benthic_snps<-read.csv("./00_raw-data/snp-array-data/benthic-array/benthic-final-snps.csv")

## get rid of bad snps
FG_bad<-c("12")
TL_bad<-c("28","33","41","36")
WB_bad<-c("50","62","61","70","52","71")
WT_bad<-c("86","74","92","83","78","96","91","90","80","93","94")
bad_snps_ben<-c(FG_bad,TL_bad,WB_bad,WT_bad)
benthic_snps<-filter(benthic_snps, !SNP.ID %in% bad_snps_ben) ## left with 74 snps with average freq of .80

## average freq by pop
benthic_freq<-benthic_snps %>% group_by(pop) %>% 
  summarise(pop_freq=mean(freq))
benthic_snps<-mutate(benthic_snps, pop_freq = case_when(pop == "FG" ~ benthic_freq[[1,2]],
                                                        pop == "TL" ~ benthic_freq[[2,2]],
                                                        pop == "WB" ~ benthic_freq[[3,2]],
                                                        pop == "WT" ~ benthic_freq[[4,2]]))
FG_freq<-benthic_freq[[1,2]]
TL_freq<-benthic_freq[[2,2]]
WB_freq<-benthic_freq[[3,2]]
WT_freq<-benthic_freq[[4,2]]

## limnetic snps 
limnetic_snps<-read.csv("./00_raw-data/snp-array-data/limnetic-array/limnetic-final-snps.csv")

## get rid of bad SNPs
LG_bad<-c("116","99","97","104")
SL_bad<-c("126","130","143","136","144","127","142")
SR_bad<-c()
WK_bad<-c("171")
bad_snps_lim<-c(LG_bad,SL_bad,SR_bad,WK_bad)
limnetic_snps<-filter(limnetic_snps, !SNP.ID %in% bad_snps_lim) ## 84 snps with average freq of .83

## average freq by pop
limnetic_freq<-limnetic_snps %>% group_by(pop) %>% 
  summarise(pop_freq=mean(freq))
limnetic_snps<-mutate(limnetic_snps, pop_freq = case_when(pop == "LG" ~ limnetic_freq[[1,2]],
                                            pop == "SL" ~ limnetic_freq[[2,2]],
                                            pop == "SR" ~ limnetic_freq[[3,2]],
                                            pop == "WK" ~ limnetic_freq[[4,2]]))
LG_freq<-limnetic_freq[[1,2]]
SL_freq<-limnetic_freq[[2,2]]
SR_freq<-limnetic_freq[[3,2]]
WK_freq<-limnetic_freq[[4,2]]

## benthic snp data 
benthic_files<-list.files("./00_raw-data/snp-array-data/benthic-array/benthic-data/", full.names = T)
benthic_data<-list()
for (i in 1:length(benthic_files)){
  benthic_data[[i]]<-read.csv(benthic_files[i]) %>% select(-any_of(c("Call.Rate","CallRate","Call Rate")))
}
benthic_data<-bind_rows(benthic_data)

## fixing loon 22 
LO22B<-filter(benthic_data,substr(sample,1,1)=="T") %>%
  mutate(sample=paste0("LO22",sample))

## combining
benthic_data<-bind_rows(benthic_data,LO22B) %>% mutate(lake=factor(substr(sample,1,2))) %>% 
  filter(lake %in% c("CC","CL","FL","HL","LL","LO","LP","RL","GL"))
colnames(benthic_data)<-str_replace_all(colnames(benthic_data), "X", "")

## limnetic snp data
limnetic_files<-list.files("./00_raw-data/snp-array-data/limnetic-array/limnetic-data/", full.names = T)
limnetic_data<-list()
for (i in 1:length(limnetic_files)){
  limnetic_data[[i]]<-read.csv(limnetic_files[i]) %>% select(-any_of(c("Call.Rate","CallRate")))
}
limnetic_data<-bind_rows(limnetic_data)

## fixing loon 22 
LO22L<-filter(limnetic_data,substr(sample,1,1)=="T") %>%
  mutate(sample=paste0("LO22",sample))

## combining
limnetic_data<-bind_rows(limnetic_data,LO22L) %>% mutate(lake=factor(substr(sample,1,2))) %>% 
  filter(lake %in% c("CC","CL","FL","HL","LL","LO","LP","RL","GL"))
colnames(limnetic_data)<-str_replace_all(colnames(limnetic_data), "X", "")

## dupes
benthic_data<-benthic_data %>% group_by(sample) %>%
  mutate(sample=paste0(sample, "_", row_number())) %>%
  ungroup()
limnetic_data<-limnetic_data %>% group_by(sample) %>%
  mutate(sample=paste0(sample, "_", row_number())) %>%
  ungroup()

## combine?
fish_data<-full_join(benthic_data,limnetic_data) %>% mutate(year=factor(substr(sample,3,4)))

## filter SNPs
fish_data<-select(fish_data, -all_of(c(bad_snps_ben,bad_snps_lim)))

## replace NCs
fish_data<-fish_data %>% mutate(across(where(is.character), ~ ifelse(.x %in% c(c("NC","No Call","Invalid")), NA, .x)))
fish_data<-fish_data %>% mutate(across(everything(), ~ str_replace_all(., ":", "")))

## counting snp calls - too lazy to clean this up
FG_snps<-filter(benthic_snps, pop=="FG")$SNP.ID
TL_snps<-filter(benthic_snps, pop=="TL")$SNP.ID
WB_snps<-filter(benthic_snps, pop=="WB")$SNP.ID
WT_snps<-filter(benthic_snps, pop=="WT")$SNP.ID
LG_snps<-filter(limnetic_snps, pop=="LG")$SNP.ID
SL_snps<-filter(limnetic_snps, pop=="SL")$SNP.ID
SR_snps<-filter(limnetic_snps, pop=="SR")$SNP.ID
WK_snps<-filter(limnetic_snps, pop=="WK")$SNP.ID
fish_data<-fish_data %>% mutate(FG_N = rowSums(!is.na(across(all_of(FG_snps)))),
                                TL_N = rowSums(!is.na(across(all_of(TL_snps)))),
                                WB_N = rowSums(!is.na(across(all_of(WB_snps)))),
                                WT_N = rowSums(!is.na(across(all_of(WT_snps)))),
                                LG_N = rowSums(!is.na(across(all_of(LG_snps)))),
                                SL_N = rowSums(!is.na(across(all_of(SL_snps)))),
                                SR_N = rowSums(!is.na(across(all_of(SR_snps)))),
                                WK_N = rowSums(!is.na(across(all_of(WK_snps)))))

## weighting
fish_data<-fish_data %>% mutate(FG_W = FG_N/rowSums(across(ends_with("_N"))),
                                TL_W = TL_N/rowSums(across(ends_with("_N"))),
                                WB_W = WB_N/rowSums(across(ends_with("_N"))),
                                WT_W = WT_N/rowSums(across(ends_with("_N"))),
                                LG_W = LG_N/rowSums(across(ends_with("_N"))),
                                SL_W = SL_N/rowSums(across(ends_with("_N"))),
                                SR_W = SR_N/rowSums(across(ends_with("_N"))),
                                WK_W = WK_N/rowSums(across(ends_with("_N")))) %>% 
                        mutate(FG_W =FG_W*FG_freq,
                                     TL_W = TL_W*TL_freq,
                                     WB_W = WB_W*WB_freq,
                                     WT_W = WT_W*WT_freq,
                                     LG_W = LG_W*LG_freq,
                                     SL_W = SL_W*SL_freq,
                                     SR_W = SR_W*SR_freq,
                                     WK_W = WK_W*WK_freq)

## snp list
snp_list<-c(benthic_snps$Base_A,limnetic_snps$Base_A)
names(snp_list)<-c(benthic_snps$SNP.ID,limnetic_snps$SNP.ID)
pops<-c("FG","TL","WB","WT","LG","SL","SR","WK")
snps_both<-bind_rows(benthic_snps,limnetic_snps)

## ancestry counts
for (i in seq_along(pops)) {
  pop_i <- pops[i]
  
  snps <- filter(snps_both, pop == pop_i)$Base_A
  names(snps) <- filter(snps_both, pop == pop_i)$SNP.ID
  
  count_col <- paste0(pop_i, "_count")
  homo_col <- paste0(pop_i, "_homo")
  
  fish_data[[count_col]] <- apply(fish_data, 1, function(row) {
    sum(sapply(names(snps), function(snp) {
      allele <- snps[[snp]]
      genotype <- row[[snp]]
      if (is.na(genotype)) return(0)
      str_count(genotype, allele)
    }))
  })
  
  fish_data[[homo_col]] <- apply(fish_data, 1, function(row) {
    sum(sapply(names(snps), function(snp) {
      allele <- snps[[snp]]
      genotype <- row[[snp]]
      if (is.na(genotype)) return(0)
      genotype == paste0(allele, allele)
    }))
  })
}

## weighting
for (pop in pops) {
  count_col <- paste0(pop, "_count")
  weight_col <- paste0(pop, "_W")
  score_col <- paste0(pop, "_score")
  fish_data[[score_col]] <- ifelse(fish_data[[weight_col]] == 0, 0, fish_data[[count_col]] / fish_data[[weight_col]])
}

## proportions
score_cols<-paste0(pops,"_score")
row_total<-rowSums(fish_data[score_cols], na.rm = TRUE)
for (pop in pops) {
  score_col <- paste0(pop, "_score")
  prop_col <- paste0(pop, "_prop")
  fish_data[[prop_col]] <- ifelse(row_total == 0, 0, fish_data[[score_col]] / row_total)
}

## checking dupes in raw data
dupes<-fish_data[str_detect(fish_data$sample,"_2"),]
dupes<-dupes%>% mutate(clean_id=substr(sample,1,8))
dupes_data<-fish_data[str_detect(fish_data$sample, str_c(dupes$clean_id, collapse = "|")), ] %>% 
  mutate(clean_id=substr(sample,1,8))
genotype_cols<-colnames(dupes_data)[2:160]

diffs<-dupes_data %>%
  group_by(clean_id) %>%
  summarise(
    across(all_of(genotype_cols), ~ {
      val1 <- .[1]
      val2 <- .[2]
      
      if (is.na(val1) && is.na(val2)) {
        "same"
      } else if (is.na(val1) || is.na(val2)) {
        "na_vs_value"
      } else if (val1 == val2) {
        "same"
      } else {
        # Compare character by character
        alleles1 <- strsplit(val1, "")[[1]]
        alleles2 <- strsplit(val2, "")[[1]]
        diff_count <- sum(alleles1 != alleles2)
        if (diff_count == 1) {
          "one_allele_diff"
        } else if (diff_count == 2) {
          "both_alleles_diff"
        } else {
          "other"  # Should not happen, safety fallback
        }
      }
    }), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    na_vs_value = sum(c_across(all_of(genotype_cols)) == "na_vs_value"),
    one_allele_diff = sum(c_across(all_of(genotype_cols)) == "one_allele_diff"),
    both_alleles_diff = sum(c_across(all_of(genotype_cols)) == "both_alleles_diff")) %>%
  select(clean_id, na_vs_value, one_allele_diff, both_alleles_diff)

## mean errors
mean_valid<-sum(!is.na(dupes_data[genotype_cols]))/nrow(dupes_data)
nocall_diffs<-mean(diffs$na_vs_value)
one_diff<-mean(diffs$one_allele_diff)
two_diff<-mean(diffs$both_alleles_diff)
error_rate<-(nocall_diffs+one_diff+two_diff)/mean_valid

## paring down cols
data<-fish_data %>% mutate(lake=factor(lake),year=factor(year)) %>% 
  select(fish_id=sample,lake,year,FG_count:WK_prop)

## function for rounding while preserving sum
round_sum<-function(x, step) {
  x <- as.numeric(x)
  up <- 1 / step
  x_scaled <- x * up
  y <- floor(x_scaled)
  indices <- tail(order(x_scaled - y), round(sum(x_scaled)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

## rounding everything but LO22 to nearest 0.25
prop_cols <- paste0(pops, "_prop")
rounded_matrix<-t(apply(data[prop_cols], 1, function(x) round_sum(x, step = 0.5)))
round_cols<-paste0(pops, "_round")
data[round_cols]<-rounded_matrix

## checking 2019
data19<-filter(data,year=="19")
weirdos19<-data19[rowSums(data19[round_cols]==0.25)+rowSums(data19[round_cols]==0.50)>0,]
## quite a few odd ones, either leave as is, round to 1, or omit them

## checking 2020
data20<-filter(data,year=="20")
weirdos20<-data20[rowSums(data20[round_cols]==0.25)+rowSums(data20[round_cols]==0.75)>0,]
## all should be 50-50 according to homozygosity - but HL45 is a bit weird still

## re-rounding
data_f0<-filter(data,year=="19"|(year=="22" & lake=="GL"))
rounded_matrix<-t(apply(data_f0[prop_cols], 1, function(x) round_sum(x, step = 1)))
data_f0[round_cols]<-rounded_matrix

data_f1<-filter(data,year=="20")
rounded_matrix<-t(apply(data_f1[prop_cols], 1, function(x) round_sum(x, step = 0.5)))
data_f1[round_cols]<-rounded_matrix

data_f2<-filter(data,year=="21" | (year=="24" & lake=="GL"))
rounded_matrix<-t(apply(data_f2[prop_cols], 1, function(x) round_sum(x, step = 0.25)))
data_f2[round_cols]<-rounded_matrix

data_f3<-filter(data,year=="22"&lake=="LO")
rounded_matrix<-t(apply(data_f3[prop_cols], 1, function(x) round_sum(x, step = 0.125)))
data_f3[round_cols]<-rounded_matrix

data<-bind_rows(data_f0,data_f1,data_f2,data_f3) %>% 
  mutate(fish_id_ext=fish_id) %>% 
  mutate(fish_id = word(fish_id_ext, 1, sep = "_"))

## checking overall sums
count_cols<-paste0(pops,"_count")
data<-data %>% mutate(count_sum = rowSums(across(all_of(count_cols)), na.rm=T),
                      score_sum = rowSums(across(all_of(score_cols)), na.rm=T))
low_sums<-filter(data,score_sum<50) ## a few sus ones maybe
high_sums<-filter(data,score_sum>400) ## these actually all look legit

## checking duplicates
dupes<-data[str_detect(data$fish_id_ext,"_2"),]
dupes_data<-data[str_detect(data$fish_id_ext, str_c(dupes$fish_id, collapse = "|")), ]
dupes_check<-dupes_data %>%
  group_by(fish_id) %>%
  mutate(rounds_match = all(apply(across(all_of(round_cols)), 2, function(x) length(unique(x)) == 1)[round_cols])) %>%
  ungroup() ## maximum deviance is 3 bases

## getting rid of dupes
data<-data[str_detect(data$fish_id_ext,"_1"),] %>% 
  select(-fish_id_ext)

## checking homozygosity for pure types
homo_cols<-paste0(pops,"_homo")
pures<-filter(data,year %in% c("20","21"))
pures<-pures[rowSums(pures[round_cols]==1)>0,]
pures<-pures[rowSums(pures[homo_cols]>1)==0,] ## LL21L051 and LO21L011 look odd

## 50/50s
fifties<-data[rowSums(data[round_cols]==0.5)==2,]
fifties<-fifties %>% filter(if_any(all_of(homo_cols), ~ .x > 2))
fifties<-fifties[rowSums(fifties[homo_cols]>0)<2,] ## a few clear errors here

## back crosses
backs<-data[rowSums(data[round_cols]==0.75)>0,]
backs[rowSums(backs[homo_cols]>0)==0,] 
backs[rowSums(backs[homo_cols]>1)==2,] ## a few errors here too

## 50/25/25
f2<-data[rowSums(data[round_cols]==0.25)>1,]
triples<-f2[rowSums(f2[round_cols]==0.50)>0,]
triples[rowSums(triples[homo_cols]>2)>1,] ## a few oddities but no clear fixes

## quads
quads<-f2[rowSums(f2[round_cols]==0.25)==4,]
quads[rowSums(quads[homo_cols]>1)>0,] ## all good

## saving
write.csv(data,file = "./02_tidy-data/final-ancestry.csv",row.names = F,quote = F)

## g lake
g_lake<-data %>% filter(year=="24") %>% select(-ends_with("_count"), -ends_with("_homo"), -ends_with("_score"))
write.csv(g_lake,file = "./02_tidy-data/g-ancestry.csv",row.names = F,quote = F)

