#### ASSORTMENT TESTS ####

## pacakges
library(tidyverse)
library(pegas)
library(adegenet)
library(hierfstat)

## constants
source(here::here("02_tidy-data/constants.R"))

## freq data
freq_data<-read.csv("02_tidy-data/assortment/freq-data.csv") %>% 
  mutate(type=factor(if_else(pure==T,"pure","admixed")))

## better hwe test
hwe_data<-freq_data %>%
  mutate(
    allele1 = substr(genotype, 1, 2),
    allele2 = substr(genotype, 3, 4)) %>%
  mutate(geno=paste0(allele1,"/",allele2)) %>%
  uncount(count_obs)

## ecotypes
allele_ecotype <- c(
  setNames(rep("B", length(benthic_sources)), benthic_sources),
  setNames(rep("L", length(limnetic_sources)), limnetic_sources))

## regions
allele_region <- c(
  setNames(rep("K", length(benthic_sources)), kenai_sources),
  setNames(rep("M", length(limnetic_sources)), matsu_sources))

## genind obj
gen_data<-hwe_data %>% 
  mutate(
    ecotype1 = allele_ecotype[allele1],
    ecotype2 = allele_ecotype[allele2],
    eco_geno = paste0(ecotype1, "/", ecotype2),
    region1 = allele_region[allele1],
    region2 = allele_region[allele2],
    reg_geno = paste0(region1, "/", region2))

## by pop
genind_df<-gen_data %>% 
  select(geno) %>% 
  df2genind(
    sep="/",
    ncode=2,
    pop=hwe_data$lake,
    ploidy=2)

## f stats
hf_data <- genind2hierfstat(genind_df)
basic_stats<-basic.stats(hf_data)
fis<- basic_stats$Fis[1,]

## tests
hw.test(genind_df, B = 10000)
pop_split<-seppop(genind_df)
hwe_by_lake<-lapply(pop_split, function(x) hw.test(x, B = 10000))

## summary
pop_summary<-data.frame(
  lake = names(basic_stats$Fis[1,]),
  He = basic_stats$Hs[1,],
  Ho = basic_stats$Ho[1,],
  Fis = basic_stats$Fis[1,],
  p_value = sapply(hwe_by_lake, function(x) x[, "Pr.exact"]),
  chi_sq   = sapply(hwe_by_lake, function(x) x[, "chi^2"]),
  df      = sapply(hwe_by_lake, function(x) x[, "df"])) %>% 
  mutate(p_adj=p.adjust(p_value, method = "BH"),
         group="source")

pop_summary %>% 
  select(lake,Ho,He,Fis,p_adj)

## by region
genind_region<-gen_data %>% 
  select(reg_geno) %>% 
  df2genind(
    sep="/",
    ncode=1,
    pop=gen_data$lake,
    ploidy=2)

## f stats
hf_reg <- genind2hierfstat(genind_region)
basic_stats_reg <-basic.stats(hf_reg)
fis_reg<- basic_stats_reg$Fis[1,]

## tests
pop_split_reg<-seppop(genind_region)
hwe_by_lake_reg<-lapply(pop_split_reg, function(x) hw.test(x, B = 10000))

## summary
reg_summary<-data.frame(
  lake = names(basic_stats_reg$Fis[1,]),
  He = basic_stats_reg$Hs[1,],
  Ho = basic_stats_reg$Ho[1,],
  Fis = basic_stats_reg$Fis[1,],
  p_value = sapply(hwe_by_lake_reg, function(x) x[, "Pr.exact"]),
  chi_sq   = sapply(hwe_by_lake_reg, function(x) x[, "chi^2"]),
  df      = sapply(hwe_by_lake_reg, function(x) x[, "df"])) %>% 
  mutate(p_adj=p.adjust(p_value, method = "BH"),
         group="region")

## by ecotype
eco_data<-gen_data %>% 
  filter(lake %in% c("LO","GL"))
genind_eco<-gen_data %>% 
  filter(lake %in% c("LO","GL")) %>%
  select(eco_geno) %>% 
  df2genind(
    sep="/",
    ncode=1,
    pop=eco_data$lake,
    ploidy=2)

## f stats
hf_eco <- genind2hierfstat(genind_eco)
basic_stats_eco <-basic.stats(hf_eco)
fis_eco<- basic_stats_eco$Fis[1,]

## tests
pop_split_eco<-seppop(genind_eco)
hwe_by_lake_eco<-lapply(pop_split_eco, function(x) hw.test(x, B = 10000))

## summary
eco_summary<-data.frame(
  lake = names(basic_stats_eco$Fis[1,]),
  He = basic_stats_eco$Hs[1,],
  Ho = basic_stats_eco$Ho[1,],
  Fis = basic_stats_eco$Fis[1,],
  p_value = sapply(hwe_by_lake_eco, function(x) x[, "Pr.exact"]),
  chi_sq   = sapply(hwe_by_lake_eco, function(x) x[, "chi^2"]),
  df      = sapply(hwe_by_lake_eco, function(x) x[, "df"])) %>% 
  mutate(p_adj=p.adjust(p_value, method = "BH"),
         group="ecotype")

## combined summary
assortment_summary<-bind_rows(pop_summary, reg_summary, eco_summary) %>% 
  select(group, lake, Ho, He, Fis, p_value, p_adj, chi_sq, df)

write.csv(assortment_summary, "02_tidy-data/assortment/assortment-summary.csv", row.names = F)

