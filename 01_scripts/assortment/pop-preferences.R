#### POPULUTION PREFERENCES ####

## packages
library(tidyverse)

## constants
source(here::here("02_tidy-data/constants.R"))

## freq data
freq_data<-read.csv("02_tidy-data/assortment/freq-data.csv") %>% 
  mutate(type=factor(if_else(pure==T,"pure","admixed")),
         o_e = count_obs - count_exp_raw) %>% 
  mutate(ecotype=factor(case_when(
    lake %in% limnetic_lakes ~ "limnetic",
    lake %in% benthic_lakes ~ "benthic",
    lake %in% both_lakes ~ "both"
  ), levels = c("benthic","limnetic","both"))) %>% 
  mutate(genotype_source=type,
         allele_1=substr(genotype,1,2),
         allele_2=substr(genotype,3,4),
         allele_ecotype_1=allele_ecotype[allele_1],
         allele_ecotype_2=allele_ecotype[allele_2],
         allele_region_1=allele_region[allele_1],
         allele_region_2=allele_region[allele_2],
         genotype_region=if_else(allele_region_1==allele_region_2,"pure","admixed"),
         genotype_ecotype=if_else(allele_ecotype_1==allele_ecotype_2,"pure","admixed"),
         genotype_ecotype_2=case_when(
           allele_ecotype_1=="B" & allele_ecotype_2=="B" ~ "Benthic",
           allele_ecotype_1=="L" & allele_ecotype_2=="L" ~ "Limnetic",
           TRUE ~ "Hybrid"),
         o_e_raw=count_obs-count_exp_raw)


## sig lakes
sig_lakes<-c("RL","FL","LP","GL")

pop_data<-data.frame()

for(source in sources) {
  source<-source
  with_df<-filter(freq_data, allele_1==source | allele_2==source) %>%
    summarize(mean_o_e=mean(o_e_raw))
  without_df<-filter(freq_data, allele_1!=source & allele_2!=source) %>%
    summarize(mean_o_e=mean(o_e_raw))
  urc<-data.frame(source=source,
                  with_mean_o_e=with_df$mean_o_e,
                  without_mean_o_e=without_df$mean_o_e)
  pop_data<-bind_rows(pop_data,urc)
}

pop_data_lakes<-data.frame()
for(source in sources) {
  for (lake_name in lakes){
  source<-source
  with_df<-filter(freq_data, allele_1==source | allele_2==source,
                  lake==lake_name) %>%
    summarize(mean_o_e=mean(o_e_raw))
  without_df<-filter(freq_data, allele_1!=source & allele_2!=source,
                     lake==lake_name) %>%
    summarize(mean_o_e=mean(o_e_raw))
  urc<-data.frame(source=source,
                  lake=lake_name,
                  with_mean_o_e=with_df$mean_o_e,
                  without_mean_o_e=without_df$mean_o_e)
  pop_data_lakes<-bind_rows(pop_data_lakes,urc)
}
}
pop_data_lakes<-pop_data_lakes %>% 
  filter(!is.na(with_mean_o_e), !is.na(without_mean_o_e))

pop_data_lakes %>% filter(lake %in% sig_lakes) %>% view()



