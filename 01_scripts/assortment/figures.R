#### FIGURES ####

## packages
library(tidyverse)
library(ape)
library(ggtree)
library(cowplot)
library(ggh4x)

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
  ), levels = c("benthic","limnetic","both")))

#### ASSORTMENT FIGS ####

## collapsing genotypes
allele_ecotype <- c(
  setNames(rep("B", length(benthic_sources)), benthic_sources),
  setNames(rep("L", length(limnetic_sources)), limnetic_sources))

allele_region <- c(
  setNames(rep("K", length(benthic_sources)), kenai_sources),
  setNames(rep("M", length(limnetic_sources)), matsu_sources))

## data
freq_plot_data<-freq_data %>% 
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
           TRUE ~ "Hybrid")) %>% 
  left_join(read.csv("02_tidy-data/assortment/assortment-residuals.csv") %>% select(lake, genotype, count_exp_model=exp)) %>% 
  select(lake, genotype, count_obs, count_exp_raw, count_exp_model, starts_with("genotype_")) %>% 
  pivot_longer(cols = starts_with("genotype_"),
               names_to = "group",
               values_to = "type",
               names_prefix = "genotype_") %>% 
  mutate(o_e_raw = case_when(
    group == "source" ~ count_obs - count_exp_raw,
    group != "source" ~ count_obs - count_exp_model
  ))

## qucik ecotype plot
eco_annos<-summary_hierarchical_eco %>% 
  mutate(lake=factor(lake, levels=both_lakes)) %>%
  filter(p_value<0.1) %>%
  mutate(x=2, 
         y=-3.5,
         label = paste0("hat(a)[e] == ", round(e_hat, 3)))

eco_means<-freq_plot_data %>% 
  filter(group=="ecotype_2", lake %in% both_lakes) %>% 
  mutate(lake=factor(lake, levels=both_lakes)) %>%
  group_by(lake, type) %>% 
  summarise(mean_o_e=mean(o_e_raw, na.rm=T), .groups="drop")

eco_assort<-freq_plot_data %>% filter(lake %in% both_lakes) %>% 
  filter(group=="ecotype_2") %>% 
  mutate(lake=factor(lake, levels=both_lakes))%>%
  ggplot(aes(x=type, y=o_e_raw, fill=type))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey50")+
  geom_jitter(shape=21, size=3, alpha=1, width=0.05)+
  geom_point(
    data = eco_means,
    aes(x = type, y = mean_o_e),
    shape = 3, size = 3, stroke = 1,
    show.legend = FALSE) +
   geom_text(
    data = eco_annos,
    aes(x, y, label = label),
    parse = T,
    inherit.aes = FALSE,
    size=4)+
  scale_fill_manual(values = eco_cols)+
  labs(y="Obs. Count - Exp. Count")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill=both_col),
        strip.text = element_text(size=12, face="bold"),
        legend.margin = margin(t = -2, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.margin = margin(t = -2, r = 0, b = 0, l = 0, unit = "pt"))+
  facet_wrap(~lake)

## saving
save_fig(eco_assort,
         filename = "03_figures/assortment/assort-eco-supp-fig.tiff",
         width = 0.5, height=0.4)

## getting rid of that stuff
freq_plot_data<-freq_plot_data %>% 
  filter(group!="ecotype_2")

## means
freq_plot_means<-freq_plot_data %>% 
  group_by(lake, group, type) %>% 
  summarise(mean_o_e=mean(o_e_raw, na.rm=T), .groups="drop")

## all lakes
freq_all <- freq_plot_data %>%
  mutate(lake = "All")

means_all <- freq_plot_data %>% 
  mutate(lake = "All") %>%
  group_by(lake, group, type) %>% 
  summarise(mean_o_e=mean(o_e_raw, na.rm=T), .groups="drop")

freq_plot_data<- bind_rows(freq_plot_data, freq_all) %>% 
  mutate(group=factor(group, levels=c("source","region","ecotype"), labels=c("Source","Region","Ecotype")),
         lake=factor(lake, levels=c("All",limnetic_lakes,benthic_lakes,both_lakes)))

freq_plot_means <- bind_rows(freq_plot_means, means_all) %>%
  mutate(group=factor(group, levels=c("source","region","ecotype"), labels = c("Source","Region","Ecotype")),
         lake=factor(lake, levels=c("All",limnetic_lakes,benthic_lakes,both_lakes)))

## annotations
assortment_df<-read.csv("02_tidy-data/assortment/assortment-summary.csv") %>% 
  mutate(group=factor(group, levels=c("source","region","ecotype"), labels = c("Source","Region","Ecotype")),
         lake=factor(lake, levels=c("All",limnetic_lakes,benthic_lakes,both_lakes))) %>% 
  filter(p_value<0.1) %>% 
  mutate(x=1.5, 
         y=-11,
         label = paste0("F[IS]==", round(Fis, 2)))

assortment_df<-read.csv("02_tidy-data/assortment/assortment-test-results.csv") %>% 
  mutate(subscript=substr(group,1,1)) %>% 
  mutate(group=factor(group, levels=c("source","region","ecotype"), labels = c("Source","Region","Ecotype")),
         lake=factor(lake, levels=c("All",limnetic_lakes,benthic_lakes,both_lakes))) %>% 
  filter(p_value<0.1) %>% 
  mutate(x=1.5, 
         y=-11,
         label = case_when(
           subscript=="s" ~ paste0("hat(a)[s]==", round(estimate, 2)),
           subscript=="e" ~ paste0("hat(a)[e]==", round(estimate, 2)),
           subscript=="r" ~ paste0("hat(a)[r]==", round(estimate, 2))))


## plotting
all_plot<-freq_plot_data %>% 
  filter(group!="Ecotype", lake!="All") %>% 
ggplot(aes(x = type, y = o_e_raw, fill = type)) +
  geom_line(
    data = filter(freq_plot_means, group!="Ecotype", lake!="All"),
    aes(x = type, y = mean_o_e, group = interaction(lake, group)),
    linewidth = 0.5,
    color = "black",
    linetype = "dashed") +
  geom_jitter(shape=21, size=2.5, alpha=0.75, width=0) +
  geom_point(
    data = filter(freq_plot_means, group!="Ecotype", lake!="All"),
    aes(x = type, y = mean_o_e),
    shape = 3, size = 4, stroke = 1,
    show.legend = FALSE) +
  geom_text(
    data = filter(assortment_df, group!="Ecotype", lake!="All"),
    aes(x, y, label = label),
    parse = T,
    inherit.aes = FALSE,
    size=3.5)+
  scale_fill_manual(values = geno_cols,
                    labels = geno_labs) +
  scale_x_discrete(labels = c("", "")) +
  labs(y = "Obs. Count - Exp. Count") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()) +
  guides(
    fill = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1,
        shape = 21)))+
  facet_grid2(
    cols = vars(lake),
    rows = vars(group),
    scales = "free_x",
    strip = strip_themed(
      background_x = elem_list_rect(fill = all_lake_cols[-1])))

all_plot_vert<-all_plot+
  theme(legend.position = "none")+
  facet_grid2(
    rows = vars(lake),
    cols = vars(group),
    scales = "free_x",
    strip = strip_themed(
      background_y = elem_list_rect(fill = all_lake_cols)))
  

## ecotype
eco_plot<-freq_plot_data %>% 
  filter(group=="Ecotype", lake %in% both_lakes) %>% 
  ggplot(aes(x = type, y = o_e_raw, fill = type)) +
  geom_line(
    data = filter(freq_plot_means, group=="Ecotype", lake %in% both_lakes),
    aes(x = type, y = mean_o_e, group = interaction(lake, group)),
    linewidth = 0.5,
    color = "black",
    linetype = "dashed") +
  geom_jitter(shape=21, size=3.5, alpha=0.75, width=0) +
  geom_point(
    data = filter(freq_plot_means, group=="Ecotype", lake %in% both_lakes),
    aes(x = type, y = mean_o_e),
    shape = 3, size = 4, stroke = 1,
    show.legend = FALSE) +
  geom_text(
    data = filter(assortment_df, group=="Ecotype"),
    aes(x, y=-3, label = label),
    parse = T,
    inherit.aes = FALSE,
    size=4.5)+
  scale_fill_manual(values = geno_cols,
                    labels = geno_labs) +
  scale_x_discrete(labels = c("", "")) +
  labs(y = "Obs. Count - Exp. Count") +
  theme(
    legend.text = element_text(size=12, color="white"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()) +
  guides(
    fill = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 0,
        shape = 21,
        color="white")))+
  facet_grid2(
    cols = vars(lake),
    rows = vars(group),
    scales = "free_x",
    strip = strip_themed(
      background_x = elem_list_rect(fill = both_col)))

eco_plot_vert<-eco_plot+
  facet_grid2(
    rows = vars(lake),
    cols = vars(group),
    scales = "free_x",
    strip = strip_themed(
      background_y = elem_list_rect(fill = both_col)))

assort_plot<-plot_grid(all_plot, eco_plot, ncol=2, rel_widths = c(3,1))

save_fig(
  plot = assort_plot,
  filename = "03_figures/assortment/assortment-fig-new.tiff",
  width = 1, height=0.33)

assort_plot_vert<-plot_grid(all_plot_vert, eco_plot_vert, ncol=1, rel_heights = c(4,1))

 
#### DIVERGENCE FIGS ####

## data
div_data<-read.csv("02_tidy-data/assortment/divergence-data.csv") %>% 
  right_join(freq_data, by="genotype") %>% 
  mutate(fst=if_else(is.na(fst),0,fst),
         fst_adapt=if_else(is.na(fst_adapt),0,fst_adapt),
         morph=if_else(is.na(morph),0,morph),
         o_e_raw=count_obs-count_exp_raw)

## plot data
div_data_all<-mutate(div_data, lake="All")
fig_data<-bind_rows(div_data,div_data_all) %>% 
  select(lake, type, o_e_raw, fst, fst_adapt, morph) %>% 
  mutate(lake=factor(lake, levels=c("All",lakes))) %>% 
  pivot_longer(cols = c(fst, fst_adapt, morph), names_to = "divergence", values_to = "value") %>%
  mutate(divergence=factor(divergence, levels=c("fst","fst_adapt","morph"),
                           labels=c("Genetic (genome-wide)","Genetic (adaptive)","Morphological"))) %>%
  filter(type=="admixed")

## quick plot
fig_data %>% 
  ggplot(aes(x=value, y=o_e_raw))+
  geom_smooth(method="lm", se=F, linetype="dashed",color="grey", linewidth=0.75)+
  geom_point(shape=21, size=2, alpha=0.75, fill="grey55")+
  labs(y="Obs. Count - Exp. Count")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        strip.background = element_rect(fill="grey95"),
        strip.text = element_text(face="bold", size=10))+
  facet_grid2(
    rows = vars(lake),
    cols = vars(divergence),
    scales = "free_x",
    strip = strip_themed(
      background_y = elem_list_rect(fill = all_lake_cols)))

## all lakes
div_1<-fig_data %>% 
  filter(lake=="All") %>% 
  ggplot(aes(x=value, y=o_e_raw))+
  geom_smooth(method="lm", se=F, linetype="dashed",color="grey50", linewidth=0.75)+
  geom_point(shape=21, size=2, alpha=0.75, fill="grey45")+
  labs(y="Obs. Count - Exp. Count")+
  theme(axis.title.y = element_text(size=12, color="white"),
        axis.title.x = element_blank())+
  facet_grid2(
    rows = vars(lake),
    cols = vars(divergence),
    scales = "free_x")

## limnetic lakes
div_2<-fig_data %>% 
  filter(lake %in% limnetic_lakes) %>% 
  ggplot(aes(x=value, y=o_e_raw))+
  geom_smooth(method="lm", se=F, linetype="dashed",color="grey50", linewidth=0.75)+
  geom_point(shape=21, size=2, alpha=0.75, fill="grey45")+
  labs(y="Obs. Count - Exp. Count")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, hjust=0.75))+
  facet_grid2(
    rows = vars(lake),
    cols = vars(divergence),
    scales = "free_x",
    strip = strip_themed(
      background_y = elem_list_rect(fill = limnetic_col),
      background_x = element_blank(),
      text_x = element_blank()))

## combining
div_12<-plot_grid(div_1, div_2, ncol=1, rel_heights = c(1.05,3))

## benthic lakes
div_3<-fig_data %>% 
  filter(lake %in% benthic_lakes) %>% 
  ggplot(aes(x=value, y=o_e_raw))+
  geom_smooth(method="lm", se=F, linetype="dashed",color="grey50", linewidth=0.75)+
  geom_point(shape=21, size=2, alpha=0.75, fill="grey45")+
  labs(y="Obs. Count - Exp. Count")+
  theme(axis.title = element_blank())+
  facet_grid2(
    rows = vars(lake),
    cols = vars(divergence),
    scales = "free_x",
    strip = strip_themed(
      background_y = elem_list_rect(fill = benthic_col),
      background_x = elem_list_rect(fill = "grey95")))

## both lakes
div_4<-fig_data %>% 
  filter(lake %in% both_lakes) %>% 
  ggplot(aes(x=value, y=o_e_raw))+
  geom_smooth(method="lm", se=F, linetype="dashed",color="grey50", linewidth=0.75)+
  geom_point(shape=21, size=2, alpha=0.75, fill="grey45")+
  labs(y="Obs. Count - Exp. Count")+
  ylim(-6,12)+
  theme(axis.title = element_blank())+
  facet_grid2(
    rows = vars(lake),
    cols = vars(divergence),
    scales = "free_x",
    strip = strip_themed(
      background_y = elem_list_rect(fill = both_col),
      background_x = element_blank(),
      text_x = element_blank()))

## combining
div_34<-plot_grid(div_3, div_4, ncol=1, rel_heights = c(4,2.5))

## all together
div_plot<-plot_grid(div_12, div_34, ncol=2, rel_widths = c(1,1))

## saving
save_fig(
  plot = div_plot,
  filename = "03_figures/assortment/div-fig-new.tiff",
  width = 1, height=0.5)

#### TREE FIGURES ####

## fst
fst<-read.csv("02_tidy-data/other-data/ak-fst.csv") %>% 
  select(first,second,fst) %>% 
  mutate(pop1 = pmin(first, second),
         pop2 = pmax(first, second)) %>% 
  mutate(genotype=paste0(pop1,pop2))

## adaptive fst
adapt_fst<-read.csv("02_tidy-data/assortment/adaptive-fst-new.csv") %>% 
  # separate pair column into first and second
  separate(pair, into = c("first", "second"), sep = "\\.") %>%
  select(first,second,fst_adapt=weighted_fst) %>% 
  mutate(pop1 = pmin(first, second),
         pop2 = pmax(first, second)) %>% 
  mutate(genotype=paste0(pop1,pop2))

## morph
morph<-read.csv("02_tidy-data/other-data/morph-ld.csv")
delta_ld<-c()
for(i in 1:nrow(fst)){
  pop1<-fst$pop1[i]
  pop2<-fst$pop2[i]
  ld1<-filter(morph,source==pop1)$ld
  ld2<-filter(morph,source==pop2)$ld
  delta_ld[i]<-abs(ld1-ld2)
}

## combining
div_data<-fst %>% mutate(morph=delta_ld) %>% 
  select(genotype,fst,morph) %>% 
  left_join(select(adapt_fst,genotype,fst_adapt), by="genotype")

## tree data
fst_matrix<-matrix(ncol=8,nrow =8)
colnames(fst_matrix)<-rownames(fst_matrix)<-sources
for (i in 1:nrow(fst)) {
  pop_1<-fst$pop1[i]
  pop_2<-fst$pop2[i]
  value<-fst$fst[i]
  fst_matrix[pop_1,pop_2]<-fst_matrix[pop_2,pop_1]<-value
}
fst_matrix[is.na(fst_matrix)] <- 0
fst_dist<-as.dist(fst_matrix)

## making tree
nj_tree<-nj(fst_dist)
rooted_tree<-root(nj_tree, outgroup = "TL", resolve.root = TRUE)
plot(rooted_tree, cex = 1.2)
p<-ggtree(rooted_tree)
p$data <- left_join(p$data, select(morph,source,ld), by = c("label" = "source"))
lim <- max(abs(range(morph$ld, na.rm = TRUE)))

## plotting the tree
tree_plot<-p +
  geom_tippoint(aes(color = ld), size = 5) +
  geom_tiplab(size = 4, hjust = -0.5) +
  scale_color_gradientn(colors = rev(source_cols_alt), limits = c(-lim, lim), name = "Body Shape LD") +
  theme_tree()+
  coord_cartesian(clip = "off")+
  theme(legend.position = "bottom",
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, margin = margin(b=12,r=10)),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.25, "cm"),
        plot.margin = margin(0, 10, 0, 0, unit = "pt"))

ggsave(tree_plot, filename = "03_figures/assortment/tree.tiff",
       device = "tiff",
       width = 700, height=600, units = "px",
       dpi=300, scale=1.5)

## adaptive tree
adapt_fst_matrix<-matrix(ncol=8,nrow =8)
colnames(adapt_fst_matrix)<-rownames(adapt_fst_matrix)<-sources
for (i in 1:nrow(adapt_fst)) {
  pop_1<-adapt_fst$pop1[i]
  pop_2<-adapt_fst$pop2[i]
  value<-adapt_fst$fst_adapt[i]
  adapt_fst_matrix[pop_1,pop_2]<-adapt_fst_matrix[pop_2,pop_1]<-value
}
adapt_fst_matrix[is.na(adapt_fst_matrix)] <- 0
adapt_fst_dist<-as.dist(adapt_fst_matrix)

## making tree
nj_tree<-nj(adapt_fst_dist)
plot(nj_tree, cex = 1.2)
p<-ggtree(nj_tree)
p$data <- left_join(p$data, select(morph,source,ld), by = c("label" = "source"))
lim <- max(abs(range(morph$ld, na.rm = TRUE)))

## plotting the tree
adapt_plot<-p +
  geom_tippoint(aes(color = ld), size = 5) +
  geom_tiplab(size = 4, hjust = -0.5) +
  scale_color_gradientn(colors = rev(source_cols_alt), limits = c(-lim, lim), name = "Body Shape LD") +
  theme_tree()+
  coord_cartesian(clip = "off")+
  theme(legend.position = "bottom",
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, margin = margin(b=12,r=10)),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.25, "cm"),
        plot.margin = margin(0, 10, 0, 0, unit = "pt"))

## div plot
fst_morph<-ggplot(div_data, aes(x=morph, y=fst_adapt))+
  geom_point(shape=21, size=3, fill="grey75", alpha=0.75)+
  geom_smooth(method="lm", se=F, linetype="dashed", color="grey25", linewidth=0.75)+
  labs(x="Morphological Divergence (Δ LD)", y="Genetic Divergence (Adaptive FST)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

plot_grid(fst_morph, adapt_plot, ncol=2, rel_widths = c(1,1))

ggsave(filename = "03_figures/assortment/adaptive-tree-supp.tiff",
       device = "tiff",
       width = 1200, height=600, units = "px",
       dpi=300, scale=2.5)


#### OLD ####

## plotting all
all_plot<-ggplot(freq_data,aes(x=count_exp_raw,y=count_obs,fill=pure))+
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey50")+
  #geom_jitter(shape=21, size=3, alpha=0.6, width=0.05, height=0.05)+
  geom_point(shape=21, size=3, alpha=0.75)+
  scale_fill_manual(values = c("grey45","grey90"), labels=c("Admixed","Pure"))+
  scale_x_sqrt(breaks=c(1,10,20,40,60))+
  scale_y_sqrt(breaks=c(1,10,20,40,60))+
  labs(x="Expected Count", y="Observed Count")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.2,0.75),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.background = element_rect(color="grey90"))

## means
source_means<-freq_data %>% 
  group_by(type) %>% 
  summarise(mean_o_e = mean(o_e, na.rm = TRUE), .groups = "drop")

## plotting diff
all_boxes<-ggplot(freq_data, aes(x = type, y = o_e, fill=type)) +
  geom_line(
    data = source_means,
    aes(x = type, y = mean_o_e, group=1),
    linewidth = 0.7,
    color = "grey50",
    linetype="dashed")+  
  geom_jitter(shape=21, size=3, alpha=0.5, width=0.15)+
  geom_point(data = source_means, aes(x = type, y = mean_o_e),
             shape = 3, size = 6, stroke = 1.25)+
  scale_fill_manual(values = c("grey90","grey25"))+
  scale_x_discrete(labels=c("Admixed","Pure"))+
  labs(y="Obs - Exp", x="Ecotype")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=12))

## source pop fig
source_fig<-plot_grid(all_plot, all_boxes, rel_widths = c(4,3))
ggsave(source_fig, filename = "03_figures/assortment/source-fig.tiff",
       device = "tiff",
       width = 1200, height=600, units = "px",
       dpi=300, scale=2.5)

## by lake
lake_means<-freq_data %>%
  group_by(lake, type) %>%
  summarise(mean_o_e = mean(o_e, na.rm = TRUE), .groups = "drop")

## lake plot
freq_data<-freq_data %>% mutate(lake=factor(lake, levels = lakes))
lake_means<-lake_means %>% mutate(lake=factor(lake, levels = lakes))

lake_points<-ggplot(freq_data, aes(x = type, y = o_e, fill=type)) +
  geom_line(
    data = lake_means,
    aes(x = as.numeric(type), y = mean_o_e, group = lake),
    linewidth = 0.5,
    color = "grey50",
    linetype="dashed")+  
  geom_jitter(shape=21, size=2, alpha=0.5, width=0.15)+
  geom_point(data = lake_means, aes(x = type, y = mean_o_e),
             shape = 3, size = 3, stroke = 1)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  scale_x_discrete(labels=c("",""))+
  labs(y="Obs - Exp", x="Ecotype")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        #legend.position.inside = c(0.9,0.45),
        #legend.title = element_blank(),
        #legend.text = element_text(size=10),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=18, color="white"),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=10),
        strip.background = element_rect(fill = "grey95"))+
  facet_wrap(~lake, ncol=5)
# remove crossbars from legend
#guides(fill = guide_legend(override.aes = list(shape = 21)))

ggsave(lake_points, filename = "03_figures/assortment/lake-fig.tiff",
       device = "tiff",
       width = 900, height=600, units = "px",
       dpi=300, scale=2)


freq_plot_data  <- bind_rows(freq_data, freq_all) %>% 
  mutate(lake = factor(lake, levels = all_plot_levels),
         point_size=if_else(lake=="All","small","big"))
means_plot_data <- bind_rows(lake_means, means_all) %>%
  mutate(lake = factor(lake, levels = all_plot_levels))

all_points<-ggplot(freq_plot_data, aes(x = type, y = o_e, fill = type)) +
  geom_line(
    data = means_plot_data,
    aes(x = as.numeric(type), y = mean_o_e, group = lake),
    linewidth = 0.5,
    color = "black",
    linetype = "dashed"
  ) +
  geom_jitter(aes(size=point_size), shape = 21, alpha = 0.65, width = 0.15) +
  scale_size_manual(values = c(2,1.5)) +
  geom_point(
    data = means_plot_data,
    aes(x = type, y = mean_o_e),
    shape = 3, size = 3, stroke = 1
  ) +
  scale_fill_manual(values = c("grey45", "grey90"),
                    labels = c("Admixed", "Pure")) +
  scale_x_discrete(labels = c("", "")) +
  labs(y = "Obs. Count - Exp. Count", x = "Ecotype") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 18, color = "white"),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 10, face="bold")) +
  facet_wrap2(
    ~ lake,
    ncol = 5,
    strip = strip_themed(
      background_x = elem_list_rect(fill = all_lake_cols),
      text_x = elem_list_text(size = 10)
    )
  )

## all fig
all_fig<-plot_grid(all_plot, all_points, rel_widths = c(5,6), nrow=1)
ggsave(all_fig, filename = "03_figures/assortment/all-fig.tiff",
       device = "tiff",
       width = 1200, height=500, units = "px",
       dpi=300, scale=2.5)


