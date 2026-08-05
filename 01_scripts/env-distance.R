#### ENVIRONMENTAL DISTANCE ####

## packages
library(tidyverse)

## data
env<-read.csv("02_tidy-data/env-data.csv")

## once I update things
var1<-c("rla","calc")
var2<-c("rla","calc","macro_total_pre","zoo_total_pre")
var3<-c("rla","calc","macro_total_post","zoo_total_post")
var4<-c("rla","calc","max_depth","temp","do","conductivity","ph","secchi_depth","tp","tn","doc","chla","macro_total_pre","zoo_total_pre")
var5<-c("rla","calc","max_depth","temp","do","conductivity","ph","secchi_depth","tp","tn","doc","chla","macro_total_post","zoo_total_post")
var_sets<-list(var1,var2,var3,var4,var5)

## distance matrices
dist_sets<-list()
for (i in 1:length(var_sets)){
  dist_sets[[i]]<-dist(scale(env[var_sets[[i]]])) %>% 
    as.matrix %>% as.data.frame()
  colnames(dist_sets[[i]])<-rownames(dist_sets[[i]])<-env$lake
}

## single df
recipients<-filter(env, lake_type=="Recipient")$lake
sources<-filter(env, lake_type=="Source")$lake
env_dists<-expand.grid(source_lake = sources,  recipient_lake= recipients) %>% as.data.frame()

for (i in seq_along(dist_sets)) {
  dist_mat <- dist_sets[[i]]
  env_dists[[paste0("dist_", i)]] <-
    mapply(function(src, rec) dist_mat[as.character(src), as.character(rec)],
           env_dists$source_lake, env_dists$recipient_lake,
           USE.NAMES = FALSE)
}

## saving
write.csv(env_dists, file="02_tidy-data/env-dists-new.csv",
          quote = F, row.names = F)

