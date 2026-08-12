#### CONSTANTS ####

## packages
library(tidyverse)
library(ggh4x)
library(broom)
library(broom.mixed)
library(purrr)
library(flextable)
library(here)
library(cowplot)
library(patchwork)

## lake factors
sources<-c("SR","WK","LG","SL","WB","FG","WT","TL")
sources_morph<-c("WK","SR","LG","SL","FG","WB","WT","TL")

benthic_sources<-c("FG","TL","WB","WT")
limnetic_sources<-c("LG","SL","SR","WK")

kenai_sources<-c("TL","WT","SL","WK")
matsu_sources<-c("FG","WB","LG","SR")

lakes<-c("CC","LL","LP","CL","FL","HL","RL","GL","LO")

limnetic_lakes<-c("CL","FL","HL","RL")
benthic_lakes<-c("CC","LL","LP")
both_lakes<-c("LO","GL")

all_lake_levels<-c("All",limnetic_lakes,benthic_lakes,both_lakes)

## colors for figures
source_cols_alt<-c("#08519c","#3182bd","#6baed6","#bdd7e7",
                   "#bae4b3","#74c476","#31a354","#006d2c")

all_lake_cols<-c("grey95",rep("#6BAED6",4),rep("#74C476",3),rep("#70B9A6",2))
lake_cols<-all_lake_cols[-1]

## ecotypes
eco_levels<-c("benthic","hybrid","limnetic")
eco_labs<-c("Benthic","Hybrid","Limnetic")
eco_cols<-c("#74C476","#70B9A6","#6BAED6")
limnetic_col<-"#6BAED6"
benthic_col<-"#74C476"
both_col<-"#70B9A6"

## for genotypes
geno_levels<-c("admixed", "pure")
geno_labs<-c("Hybrid", "Pure")
geno_cols<-c("grey55","grey95")

## figure themes
theme_set(
  theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size=10),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size=12),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(size = 10, face = "bold")))

## function for saving figures
save_fig<-function(plot, filename, width, height){
  ggsave(
    plot = plot,
    filename = filename,
    device = "tiff",
    width = width*6.5,
    height = height*6.5,
    units = "in",
    scale = 2,
    dpi = 300)
}
