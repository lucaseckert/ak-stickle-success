#### ENVIRONMENTAL DATA CLEANING ####

## packages
library(tidyverse)

## factors
source_lakes<-c("SR","WK","LG","SL","WB","FG","WT","TL")
recipient_lakes<-c("CL","FL","HL","RL","CC","LL","LP","LO","GL")

#### CLEANING MACRO INVERT DATA ####

## raw data files
macro_18<-read.csv("00_raw-data/env-data/macroinverts-2018.csv")
macro_rec_18_19_22<-read.csv("00_raw-data/env-data/macroinverts-recipient-2018-2019-2022.csv")

## cleaning
macro_18_clean<-macro_18 %>% 
  separate(col=sample, into=c("country", "lake", "year", "sample"), sep=" ", remove=F) %>% 
  mutate(lake = if_else(lake == "FN", "FG", lake)) %>% 
  mutate(lake = if_else(lake == "G", "GL", lake)) %>% 
  select(lake, year, gammaridae=gamaridae, chironomidae) %>% 
  group_by(lake,year) %>% 
  summarize(across(c(gammaridae, chironomidae), mean, na.rm = T)) %>% 
  filter(lake %in% c(source_lakes, recipient_lakes))

macro_rec_clean<-macro_rec_18_19_22 %>% 
  separate(col=sample, into=c("country", "year", "lake", "sample"), sep="-", remove=F) %>% 
  mutate(lake = if_else(lake == "G", "GL", lake)) %>% 
  select(lake, year, gammaridae = Amphipoda, chironomidae = Chironomidae) %>% 
  group_by(lake,year) %>% 
  summarize(across(c(gammaridae, chironomidae), mean, na.rm = T)) %>% 
  filter(lake %in% c(source_lakes, recipient_lakes))

macro_clean<-bind_rows(macro_18_clean,macro_rec_clean) %>% 
  distinct() %>% 
  mutate(time=if_else(year=="22","post","pre")) %>% 
  select(-year) %>% 
  group_by(lake, time) %>% 
  summarize(across(everything(), mean, na.rm=T))

#### CLEANING ZOOPLANKTON DATA ####

## raw data
zoo_rec_18_19_22<-read.csv("00_raw-data/env-data/zooplankton-recipient-2018-2019-2022.csv")
zoo_19<-read.csv("00_raw-data/env-data/zooplankton-source-2019.csv")
zoo_18<-read.csv("00_raw-data/env-data/zooplankton-2018.csv")

zoo_jen_18_19<-read.csv("00_raw-data/env-data/zooplankton-jen-2018-2019.csv")
zoo_jen_22<-read.csv("00_raw-data/env-data/zooplankton-2022.csv")

## cleaning
zoo_rec_18_19_22[is.na(zoo_rec_18_19_22)]<-0
zoo_rec_clean<-zoo_rec_18_19_22 %>% 
  separate(col=lake, into=c("country", "month_year", "lake"), sep="-", remove=F) %>% 
  mutate(lake = if_else(lake == "G", "GL", lake)) %>% 
  mutate(month_year=str_remove(month_year,"JUL")) %>% 
  mutate(year=as.numeric(str_remove(month_year,"J"))) %>% 
  mutate(calanoid_copepods=CalCopepodids+H..septentrionalis+L.pribilofensis+Epichura.sp.,
         daphnia=D.longiremis+D..middendorffiana+D.mendotae) %>% 
  select(lake,year,calanoid_copepods,daphnia) %>% 
  filter(year != 22)

## cleaning 19
zoo_19_clean<-zoo_19 %>% 
  mutate(lake=case_when(
    lake=="Walby " ~ "WB",
    lake=="Long Lake" ~ "LG",
    lake=="South Rolley" ~ "SR",
    lake=="Wik Lake" ~ "WK",
    lake=="Tern Lake" ~ "TL",
    lake=="Finger" ~ "FG",
    lake=="Watson" ~ "WT",
    lake=="Spirit Lake" ~ "SL")) %>% 
  filter(lake %in% source_lakes) %>% 
  mutate(year=19) %>% 
  mutate(calanoid_copepods=Juvénile+Epischura.sp+Leptodiaptomus.pribilofensis+Unknown.female,
         daphnia=Daphnia.mendotae+Daphnia.longiremus) %>% 
  select(lake,year,calanoid_copepods,daphnia)

## cleaning 18
zoo_18_clean<-zoo_18 %>% 
  mutate(lake=if_else(lake=="G","GL",lake)) %>% 
  filter(lake %in% c(source_lakes, recipient_lakes)) %>% 
  mutate(year=18) %>% 
  mutate(calanoid_copepods=CalCopepedids+Copepe.L.pribilo+Copepe.Epichura+Cal.femelle+Copepe.H.Sept+
           X..H.septentrionalis+X..H..septentrionalis+X..L.pribilo+X..L.pribilofensis+X..Epichura+X..Epichura.sp.,
         daphnia= D..juv+D.longiremis+D..middendorffiana+D.mandotae) %>% 
  select(lake,year,calanoid_copepods,daphnia)

## cleaning jen
zoo_jen_18_19_clean<-zoo_jen_18_19 %>% 
  separate(col=Sample, into=c("country", "lake", "year"), sep="-", remove=F) %>% 
  mutate(lake=if_else(lake=="G","GL",lake)) %>% 
  mutate(year=as.numeric(year)) %>%
  mutate(calanoid_copepods=H..septentrionalis+L.pribilofensis+Epichura.sp...lacustris.,
         daphnia=Daphnia..Juv..unknown.+D.longiremis+D..middendorffiana+D.mendotae) %>% 
  select(lake,year,calanoid_copepods,daphnia)

zoo_jen_22_clean<-zoo_jen_22 %>% 
  separate(col=Sample, into=c("country", "lake", "sample"), sep="-", remove=F) %>% 
  mutate(lake=if_else(lake=="G","GL",lake)) %>% 
  mutate(calanoid_copepods=Calanoid.juvenile+Hetercope.septentrionalis+Epischura.lacrustis+Leptodiaptomus.minutus+Leptodiaptomus.pribilofensis+Calanoid.female..small.+Calanoid.female..large.+Unidentified.Calanoid,
         daphnia=Daphnia.longiremus+Unidentified.Daphnia) %>% 
  select(lake,calanoid_copepods,daphnia) %>% 
  group_by(lake) %>% 
  summarize(across(everything(), mean, na.rm=T)) %>% 
  mutate(year=22)

## combining
zoo_clean<-bind_rows(zoo_18_clean,zoo_19_clean, zoo_rec_clean, zoo_jen_18_19_clean, zoo_jen_22_clean) %>% 
  group_by(lake,year) %>% 
  summarize(across(everything(), mean, na.rm=T)) %>% 
  ungroup() %>%
  mutate(time=if_else(year==22,"post","pre")) %>% 
  select(-year) %>% 
  group_by(lake,time) %>% 
  summarize(across(everything(), mean, na.rm=T))

#### CLEANING WATER CHEM DATA ####

## 2018
water_chem_18<-read.csv("00_raw-data/env-data/water-chemistry-2018.csv") %>% 
  select(lake=Abbreviation,
         region=Region,
         lat=Latitude,
         lon=Longitude,
         surface_area=Surface..ha.,
         max_depth=Max.depth..m.,
         volume=Volume..m3.,
         temp=temperature.C.at.0.5.m,
         do=O2..mg.L..at..0.5.m,
         conductivity=Conductivity..uS.cm..at.0.5.m,
         ph=pH.at.0.5.m,
         secchi_depth=Secchi.depth..m.,
         tp=TP.ug.L,
         tn=TN.mg.L,
         doc=DOC.mg.L,
         chla=chla.ug.L,
         calc=Ca.mg.L) %>% 
  mutate(lake=if_else(lake=="FG or FN","FG",lake)) %>% 
  filter(lake%in%c(source_lakes, recipient_lakes)) %>% 
  mutate(secchi_depth=if_else(secchi_depth=="bottom (5.1)","5.1",secchi_depth)) %>% 
  mutate(region=if_else(region=="Tote Rd - Kenai","Kenai",region)) %>% 
  mutate(across(!c(lake, region), as.numeric)) %>% 
  mutate(year=18)

## 2019
water_chem_19<-read.csv("00_raw-data/env-data/water-chemistry-2019.csv") %>% 
  select(lake,
         temp=Temp..C..at.0.5m,
         do=O2..mg.L..at.0.5m,
         conductivity=Conductivity..uS.cm..at.0.5m,
         ph=pH.at.0.5m,
         secchi_depth=Secchi.depth..m.,
         tp=TP..ug.L.,
         tn=TN..ppm.,
         doc=DOC..mg.L.,
         chla=Chla..ug.L.) %>% 
  filter(lake %in% c(source_lakes,recipient_lakes)) %>% 
  mutate(year=19)

## combining
water_chem<-bind_rows(water_chem_18,water_chem_19) %>% 
  group_by(lake) %>% 
  summarize(across(!c(region), mean, na.rm = T))
 
## old data
old_env<-read.csv("00_raw-data/env-data/old-env-data.csv") %>% 
  select(lake, lake_type=lake.type, ecotype, region, lat, lon=long, area, calc, rla=litt.prop) 

## combining
water_chem_final<-left_join(old_env, select(water_chem,-c(lat,lon,surface_area,volume,calc,year)))

#### COMBINING AND SAVING ####

## combining prey
prey_final<-left_join(macro_clean, zoo_clean, by=c("lake","time")) %>% 
  pivot_wider(values_from=c(gammaridae,chironomidae,calanoid_copepods,daphnia) , names_from = "time") %>% 
  mutate(across(ends_with("_post"), ~ if_else(is.na(.x),
                                              get(sub("_post$", "_pre", cur_column())),
                                              .x))) %>% 
  ungroup() %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.x)))) %>% 
  mutate(macro_total_pre=gammaridae_pre + chironomidae_pre,
         macro_total_post=gammaridae_post + chironomidae_post,
         zoo_total_pre=calanoid_copepods_pre + daphnia_pre,
         zoo_total_post=calanoid_copepods_post + daphnia_post)

## final data
env_final<-left_join(water_chem_final,prey_final)
write.csv(env_final, file = "02_tidy-data/env-data.csv", quote = F, row.names = F)
