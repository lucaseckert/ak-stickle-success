#### AGE CLASSES ####

## packages 
library(tidyverse)
library(cowplot)

## constants
source(here::here("02_tidy-data/constants.R"))

## ancestry data
data<-read.csv("02_tidy-data/final-ancestry.csv") %>% 
  filter(year %in% c(20,21) | (lake=="GL" & year==24)) %>% 
  select(fish_id,lake,year,ends_with("round")) %>% 
  mutate(year=case_when(
    lake=="GL" ~ 2,
    lake!="GL" & year==20 ~ 1,
    lake!="GL" & year==21 ~ 2))

## std length data
length_lethal<-read.csv("02_tidy-data/other-data/lethal-data.csv") %>% 
  rename(fish_id=fishID) %>% 
  select(fish_id,std_length)

length_g<-read.csv("02_tidy-data/other-data/gill-rakers.csv") %>% 
  mutate(fish_id=paste0("GL24H",GL24L)) %>% 
  rename(std_length=SL) %>% 
  select(fish_id, std_length)

length_data<-bind_rows(length_lethal,length_g)

data<-left_join(data,length_data)

## year 1
data_1<-filter(data, year==1)

## year 2
data_2 <- data %>% filter(year == 2) %>%
  mutate(f2 = if_any(ends_with("round"), ~ .x == 0.25))

## plot by length
f2_plot<-ggplot(data_2, aes(x=std_length, fill=f2))+
  geom_histogram()+
  scale_fill_manual(values = c("grey","black"), labels=c("Unknown","Confirmed F2"))+
  labs(x="Standard Length (mm)", y="Count")+
  theme_bw()+
  facet_wrap(~lake)+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12, face="bold"),
        strip.background = element_rect(fill = "grey95"))

## cluster by std length
data_2 <- data_2 %>%
  filter(!is.na(std_length)) %>%
  group_by(lake) %>%
  group_modify(~ {
    km <- kmeans(.x$std_length, centers = 2)
    means <- tapply(.x$std_length, km$cluster, mean)
    clust <- ifelse(km$cluster == which.min(means), 1, 2)
    .x$cluster <- clust
    .x
  }) %>%
  ungroup()

## plot clusters
clust_plot<-ggplot(data_2, aes(x=std_length, fill=factor(cluster)))+
  geom_histogram()+
  scale_fill_manual(values = c("black","grey"), labels=c("Small Cluster","Large Cluster"))+
  labs(x="Standard Length (mm)", y="Count")+
  theme_bw()+
  facet_wrap(~lake)+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12, face="bold"),
        strip.background = element_rect(fill = "grey95"))

f1_plot<-plot_grid(f2_plot,clust_plot)
ggsave(f1_plot, filename = "03_figures/assortment/cluster-fig.tiff",
       device = "tiff",
       width = 1200, height=600, units = "px",
       dpi=300, scale=3)

## 2020 data
length_20<-read.csv("02_tidy-data/assortment/length-data-2020.csv")

data_1<-filter(data, year==1)

data_1<-data_1 %>%
  mutate(specimen=substr(fish_id,5,8),
         specimen_new=str_replace_all(specimen,"O","0"),
         fish_id=paste0(substr(fish_id,1,4),specimen_new)) %>%
  select(-specimen,-specimen_new)

data_1<-data_1 %>% select(-std_length) %>% 
  full_join(select(length_20,fish_id,std_length,bent), by="fish_id") %>% 
  mutate(f1=if_any(ends_with("round"), ~ .x == 0.5))

data_1 %>% filter(!is.na(std_length), !is.na(f1)) %>% 
  ggplot(aes(x=std_length, fill=f1))+
  geom_histogram()+
  scale_fill_manual(values = c("grey","black"), labels=c("Unknown","Confirmed F1"))+
  labs(x="Standard Length (mm)", y="Count")+
  theme_bw()+
  facet_wrap(~lake)+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12, face="bold"),
        strip.background = element_rect(fill = "grey95"))

ggsave(filename = "03_figures/assortment/cluster-2020-fig.tiff",
       device = "tiff",
       width = 800, height=600, units = "px",
       dpi=300, scale=3)

## outliers
data_1<-data_1 %>% 
  mutate(outlier=case_when(
    lake=="CC" & std_length>60 ~ T,
    lake=="FL" & std_length>70 ~ T,
    TRUE ~ F))

## final f1 data
data_f1<-data_2 %>% 
  filter(lake %in% c("LP","FL","CC") |
           !lake %in% c("LP","FL","CC") & cluster==2) %>% 
  filter(f2==F) %>% 
  bind_rows(filter(data_1, outlier==F, !is.na(FG_round)))

## saving
write.csv(data_f1, file = "02_tidy-data/assortment/final-f1-data.csv",
          quote = F, row.names = F)
