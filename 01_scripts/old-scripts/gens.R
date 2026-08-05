#### Generation Binning ####

## packages
library(tidyverse)

## ancestry data
ancestry_data<-read.csv("./02_tidy-data/final-ancestry.csv")

## for factors
sources<-c("SR","WK","LG","SL","WB","FG","WT","TL")
lakes<-c("CL","FL","HL","RL","CC","LL","LP","LO","GL")
limnetic_lakes<-c("CL","FL","HL","RL")
benthic_lakes<-c("CC","LL","LP")
both_lakes<-c("LO","GL")

## fix year
ancestry_data<-ancestry_data %>% mutate(year=case_when(
  year==19 ~ 0,
  year==20 ~ 1,
  year==21 ~ 2,
  year==22 & lake!="GL" ~ 3,
  year==22 & lake=="GL" ~ 0.5,
  year==23 & lake=="GL" ~ 1,
  year==24 & lake=="GL" ~ 2
))

## lethal data
lethal_data<-read.csv("./02_tidy-data/lethal-data.csv") %>% 
  rename(fish_id=fishID)

data<-ancestry_data %>% filter(year %in% c(1,2,3)) %>% 
  select(fish_id,lake,year, ends_with("_round")) %>% 
  left_join(select(lethal_data, -year))


data2<-data %>% filter(year==2) %>% 
  mutate(gen = factor(case_when(
    if_any(ends_with("round"), ~ .x == 0.25) ~ "F2",
    #if_any(ends_with("round"), ~ .x == 0.5) ~ "F1",
    TRUE ~ NA_character_)))

ggplot(data2, aes(x=std_length, fill=gen))+
  geom_histogram()+
  facet_wrap(~lake)

marqus<-data %>% filter(lake %in% c("LO","GL")) %>% 
  select(fish_id,lake,year,ends_with("round"))

write.csv(marqus, file = "./02_tidy-data/ancestry-data-marqus.csv", quote=F, row.names = F)

lethal19<-filter(lethal_data, year==2019) %>% 
  mutate(lake=substr(fish_id,1,2))

count(lethal19,lake,gravid)

