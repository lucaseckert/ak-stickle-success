#### EXPECTED FREQUENCIES ####

## packages
library(tidyverse)
library(cowplot)

## constants
source(here::here("02_tidy-data/constants.R"))

## all f1 data
data_f1<-read.csv("02_tidy-data/assortment/final-f1-data.csv")

## sums
source_sums<-data_f1 %>%
  group_by(lake) %>%
  summarise(across(ends_with("round"), \(x) sum(x, na.rm = TRUE)))

## counts
lake_counts<-data_f1 %>% 
  group_by(lake) %>% 
  summarise(n=n())

## proportions
source_props<-source_sums %>% 
  rowwise() %>%
  mutate(across(ends_with("round"),
                ~ .x / sum(c_across(ends_with("round")), na.rm = TRUE))) %>%
  ungroup()
colnames(source_props)<-str_remove(colnames(source_props),"_round")

## expected frequencies
freq_expected <- source_props %>%
  rowwise() %>%
  do({
    lake <- .$lake
    alleles <- names(.)[-1]          
    freqs <- unlist(.[-1])
    combos <- expand.grid(alleles, alleles, stringsAsFactors = FALSE)
    combos <- combos %>%
      mutate(genotype = ifelse(Var1 < Var2,
                               paste0(Var1, Var2),
                               paste0(Var2, Var1)),
             exp_freq = ifelse(Var1 == Var2,
                               freqs[match(Var1, alleles)]^2,
                               2 * freqs[match(Var1, alleles)] * freqs[match(Var2, alleles)])) %>%
      distinct(genotype, .keep_all = TRUE) %>%
      select(genotype, exp_freq)
    tibble(lake = lake) %>%
      bind_cols(pivot_wider(combos, names_from = genotype, values_from = exp_freq))
  }) %>%
  ungroup() %>% 
  pivot_longer(cols = -lake, names_to = "genotype", values_to = "freq")


## observed frequencies
colnames(data_f1)<-str_remove(colnames(data_f1),"_round")
obs_genos <- data_f1 %>% select(-c(fish_id,year,f2,cluster,std_length)) %>%
  rowwise() %>%
  mutate(
    alleles = list(names(pick(-lake))[pick(-lake) > 0]),
    genotype = if(length(alleles) == 2) {
      paste0(sort(alleles), collapse = "")
    } else if(length(alleles) == 1) {
      paste0(alleles, alleles)
    } else {
      NA_character_
    }
  ) %>%
  ungroup() %>%
  select(lake, genotype)

## frequencies
freq_observed <- obs_genos %>%
  group_by(lake, genotype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(lake) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

## joining
freq_data<-left_join(
  rename(freq_expected, freq_exp=freq),
  rename(freq_observed, freq_obs=freq))
freq_data[is.na(freq_data)]<-0

## cleaning
freq_data<-freq_data %>% 
  select(-n) %>% 
  mutate(pure=if_else(substr(genotype,1,2)==substr(genotype,3,4),T,F)) %>% 
  filter(freq_exp>0) %>% 
  mutate(dev_freq=freq_obs-freq_exp) %>% 
  left_join(lake_counts) %>% 
  mutate(
    count_exp=round(freq_exp*n),
    count_exp_raw=freq_exp*n,
    count_obs=freq_obs*n,
    o_e=count_obs-count_exp)

## saving
write.csv(freq_data, file = "02_tidy-data/assortment/freq-data.csv",
          quote = F, row.names = F)

## plotting all
all_plot<-ggplot(freq_data,aes(x=sqrt(count_exp),y=sqrt(count_obs),fill=pure))+
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey50")+
  geom_jitter(shape=21, size=3, alpha=0.6, width=0.05, height=0.05)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x="Expected Count (sqrt)", y="Observed Count (sqrt)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.225,0.875),
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=16))

## pures
pures<-freq_data %>% filter(pure==T) %>% 
  mutate(source=factor(substr(genotype,1,2), levels = sources))

## pure plot
pure_plot<-ggplot(pures, aes(x=sqrt(count_exp),y=sqrt(count_obs),fill=source))+
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey")+
  geom_jitter(shape=21, size=3, alpha=0.75, width=0.1, height=0.1)+
  scale_fill_manual(values = source_cols_alt)+
  labs(x="Expected Count (sqrt)", y="Observed Count (sqrt)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 14),
        axis.text = element_text(size=12))
