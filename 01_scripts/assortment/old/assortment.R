#### ASSORTMENT TESTS ####

## packages
library(tidyverse)
library(cowplot)
library(ape)
library(ggtree)

## for factors
sources<-c("SR","WK","LG","SL","WB","FG","WT","TL")
benthic_sources<-c("FG","TL","WB","WT")
limnetic_sources<-c("LG","SL","SR","WK")
lakes<-c("CC","LL","LP","CL","FL","HL","RL","GL","LO")
limnetic_lakes<-c("CL","FL","HL","RL")
benthic_lakes<-c("CC","LL","LP")
both_lakes<-c("LO","GL")
source_cols_alt<-c("#08519c","#3182bd","#6baed6","#bdd7e7",
                   "#bae4b3","#74c476","#31a354","#006d2c")

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

## all f1 data
data_f1<-data_2 %>% 
  filter(lake %in% c("LP","FL","CC") |
           !lake %in% c("LP","FL","CC") & cluster==2) %>% 
  filter(f2==F) %>% 
  bind_rows(data_1)

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

freq_data<-freq_data %>% 
  select(-n) %>% 
  mutate(pure=if_else(substr(genotype,1,2)==substr(genotype,3,4),T,F)) %>% 
  filter(freq_exp>0) %>% 
  mutate(dev_freq=freq_obs-freq_exp) %>% 
  left_join(lake_counts) %>% 
  mutate(
    count_exp=round(freq_exp*n),
    count_obs=freq_obs*n,
    std_resid = ((count_obs+1) - (count_exp+1)) / sqrt((count_exp+1) * (1 - freq_exp)),
    log_dev = log((count_obs + 1) / (count_exp + 1)),
    o_e=count_obs-count_exp)

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

#### tests ####

## tests
t.test(log_dev ~ pure, data = freq_data)
wilcox.test(log_dev ~ pure, data = freq_data)

t.test(o_e ~ pure, data = freq_data)
wilcox.test(o_e ~ pure, data = freq_data)

## means
source_means<-freq_data %>% 
  group_by(pure) %>% 
  summarise(mean_o_e = mean(o_e, na.rm = TRUE), .groups = "drop")

## plotting dev
all_boxes<-ggplot(freq_data, aes(x = pure, y = o_e, fill=pure)) +
  geom_line(
    data = source_means,
    aes(x = pure, y = mean_o_e, group=1),
    linewidth = 0.7,
    color = "grey50",
    linetype="dashed")+  
  #geom_violin(alpha=0.25)+
  geom_jitter(shape=21, size=3, alpha=0.5, width=0.15)+
  geom_point(data = source_means, aes(x = pure, y = mean_o_e),
             shape = 3, size = 6, stroke = 1.25)+
  scale_fill_manual(values = c("grey90","grey25"))+
  scale_x_discrete(labels=c("Admixed","Pure"))+
  labs(y="Obs - Exp", x="Ecotype")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=16))

## source pop fig
source_fig<-plot_grid(all_plot, all_boxes, rel_widths = c(4,3))
ggsave(source_fig, filename = "03_figures/assortment/source-fig.tiff",
       device = "tiff",
       width = 1200, height=600, units = "px",
       dpi=300, scale=2.5)

## by lake
lake_means<-freq_data %>%
  group_by(lake, pure) %>%
  summarise(mean_o_e = mean(o_e, na.rm = TRUE), .groups = "drop")

## lake plot
freq_data<-freq_data %>% mutate(lake=factor(lake, levels = lakes))
lake_means<-lake_means %>% mutate(lake=factor(lake, levels = lakes))
lake_points<-ggplot(freq_data, aes(x = pure, y = o_e, fill=pure)) +
  geom_line(
    data = lake_means,
    aes(x = as.numeric(pure)+1, y = mean_o_e, group = lake),
    linewidth = 0.5,
    color = "grey50",
    linetype="dashed")+  
  geom_point(shape=21, size=2, alpha=0.5)+
  geom_point(data = lake_means, aes(x = pure, y = mean_o_e),
    shape = 3, size = 3, stroke = 1)+
  scale_fill_manual(values = c("grey90","grey25"))+
  scale_x_discrete(labels=c("Admixed","Pure"))+
  labs(y="Obs - Exp", x="Ecotype")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=14),
        strip.text = element_text(size=12, face="bold"),
        strip.background = element_rect(fill = "grey95"))+
  facet_wrap(~lake)

ggsave(lake_points, filename = "03_figures/assortment/lake-fig.tiff",
       device = "tiff",
       width = 1000, height=1000, units = "px",
       dpi=300, scale=2)

## tests
tests_by_lake<-freq_data %>%
  group_by(lake) %>%
  nest() %>%
  mutate(wilcox_test = map(data, ~ wilcox.test(o_e ~ pure, data = .x)),
         p_value = map_dbl(wilcox_test, "p.value"))

dis_lakes<-c("FL","GL")

#### ECOTYPE ####

## assortment by ecotype
eco_data<-freq_data %>% 
  filter(lake %in% c("GL","LO")) %>% 
  mutate(source1=substr(genotype,1,2),
         source2=substr(genotype,3,4)) %>% 
  mutate(eco_pure=if_else(
    (source1 %in% benthic_sources & source2 %in% limnetic_sources) |
      (source2 %in% benthic_sources & source1 %in% limnetic_sources),
    F,T))

## plotting points
eco_points<-ggplot(eco_data, aes(x=sqrt(count_exp),y=sqrt(count_obs),fill=eco_pure))+
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey")+
  geom_jitter(aes(shape=lake), size=3, alpha=0.6, width=0.05, height=0.05)+
  scale_shape_manual(values = c(21,23))+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x="Expected Count (sqrt)", y="Observed Count (sqrt)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))

## eco means
eco_means<-eco_data %>% 
  group_by(lake,eco_pure) %>% 
  summarise(mean_o_e = mean(o_e, na.rm = TRUE), .groups = "drop")

## boxes
eco_boxes<-ggplot(eco_data, aes(x=eco_pure, y=o_e))+
  geom_line(
    data = eco_means,
    aes(x = as.numeric(eco_pure)+1, y = mean_o_e, group = lake),
    linewidth = 0.7,
    color = "grey50",
    linetype="dashed")+
  geom_jitter(aes(shape=lake, fill=eco_pure), size=3, alpha=0.5, width=0.05, height=0)+
  geom_point(data = eco_means, aes(x = eco_pure, y = mean_o_e),
             shape = 3, size = 5, stroke = 1)+
  scale_fill_manual(values = c("grey90","grey25"))+
  scale_shape_manual(values = c(21,23))+
  scale_x_discrete(labels=c("Admixed","Pure"))+
  labs(y="Obs - Exp", x="Ecotype")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=16),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_rect(fill = "grey95"))+
  facet_wrap(~lake)

eco_fig<-plot_grid(eco_points,eco_boxes, rel_widths = c(2,3))
ggsave(eco_fig, filename = "03_figures/assortment/eco-fig.tiff",
       device = "tiff",
       width = 1440, height=600, units = "px",
       dpi=300, scale=2.5)

## tests
wilcox.test(o_e ~ eco_pure, data = filter(eco_data, lake=="GL"))
wilcox.test(o_e ~ eco_pure, data = filter(eco_data, lake=="LO"))

chisq.test(x=filter(freq_data, lake=="FL")$count_obs, p=filter(freq_data, lake=="FL")$freq_exp)
chisq.test(x=filter(freq_data, lake=="HL")$count_obs, p=filter(freq_data, lake=="HL")$freq_exp)


#### DIVERGENCE ####

## fst
fst<-read.csv("02_tidy-data/other-data/ak-fst.csv") %>% 
  select(first,second,fst) %>% 
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
  select(genotype,fst,morph)

## plotting
lm(morph~fst, data = div_data) %>% summary()
ggplot(div_data, aes(x=fst,y=morph))+
  geom_point(alpha=0.5, size=3)+
  labs(x = expression("Genetic Divergence (" * F[ST] * ")"),
       y = expression("Morphological Divergence (" * Delta * "LD" * ")"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12))

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

## adding to all
freq_data<-freq_data %>% 
  left_join(div_data)
freq_data$fst[is.na(freq_data$fst)]<-0
freq_data$morph[is.na(freq_data$morph)]<-0

## plotting morph
ggplot(freq_data, aes(x=morph, y=o_e))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.6, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x=expression("Morphological Divergence (" * Delta * "LD" * ")"),
       y="Obs - Exp")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=14))

## just trendy
lm(o_e~morph, data=filter(freq_data, lake=="GL")) %>% summary()
dis_lakes<-c("GL")
morph_all<-freq_data %>% filter(lake %in% dis_lakes) %>% 
ggplot(aes(x=morph, y=o_e))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.6, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x=expression("Morphological Divergence (" * Delta * "LD" * ")"),
       y="Obs - Exp")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text = element_text(size=14))

## adding to eco data
eco_data<-eco_data %>% 
  left_join(div_data)
eco_data$fst[is.na(eco_data$fst)]<-0
eco_data$morph[is.na(eco_data$morph)]<-0

## eco hybrids
eco_data %>% filter(eco_pure==F) %>% 
ggplot(aes(x=morph, y=log_dev))+
  geom_smooth(aes(group=lake, linetype = lake),method = "lm", se=F, color="grey50", linewidth=0.5)+
  scale_linetype_manual(values = c("solid","dashed"))+
  geom_point(aes(shape=lake), alpha=1, size=3, fill="grey95")+
  scale_shape_manual(values = c(21,23))+
  labs(x=expression("Morphological Divergence (" * Delta * "LD" * ")"),
       y="Log Ratio (obs/exp)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=14))+
  guides(shape = guide_legend(order = 1),
         linetype = guide_legend(order = 2))

lm(log_dev~morph*lake, data=filter(eco_data, eco_pure==F)) %>% summary()
lm(log_dev~morph, data=filter(eco_data, eco_pure==F, lake=="GL")) %>% summary()

## just GL
lm(o_e~morph, data = filter(eco_data, eco_pure==F, lake=="GL")) %>% summary()
morph_high<-eco_data %>% filter(eco_pure==F, lake=="GL") %>% 
  ggplot(aes(x=morph, y=o_e))+
  geom_smooth(method = "lm", se=F, color="grey50", linewidth=0.5, linetype="dashed")+
  geom_point(shape=21, alpha=1, size=3, fill="grey95")+
  labs(x="",
       y="Log Ratio (obs/exp)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text = element_text(size=14))

## within eco
lm(log_dev~morph, data=filter(freq_data, morph<1)) %>% summary()
freq_data %>% filter(morph<1) %>% 
  ggplot(aes(x=morph, y=log_dev))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.6, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x=expression("Morphological Divergence (" * Delta * "LD" * ")"),
       y="Log Ratio (obs/exp)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=14))

## without pures
lm(log_dev~morph, data=filter(freq_data, morph<1, pure==F)) %>% summary()
freq_data %>% filter(morph<1, pure==F) %>% 
  ggplot(aes(x=morph, y=log_dev))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.6, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x=expression("Morphological Divergence (" * Delta * "LD" * ")"),
       y="Log Ratio (obs/exp)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=14))

## just trendy lakes
lm(o_e~morph, data=filter(freq_data, morph<1, lake %in% dis_lakes)) %>% summary()
lm(o_e~morph, data=filter(freq_data, morph<1, lake %in% dis_lakes, pure==F)) %>% summary()
morph_low<-freq_data %>% filter(morph<1, lake %in% dis_lakes) %>% 
  ggplot(aes(x=morph, y=o_e))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.6, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x="",
       y="Obs - Exp")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14))

## MORPH PLOT
morph_plot <- plot_grid(
  morph_low, NULL, morph_all, NULL, morph_high,
  ncol = 5,
  rel_widths = c(1.1, 0.1, 1, 0.1, 1),
  align = 'h'
)

## other approach
all_morph_data<-freq_data %>% filter(lake %in% dis_lakes) %>% 
  select(lake,pure,o_e,morph) %>% 
  mutate(plot="all")

low_morph_data<-freq_data %>% filter(lake %in% dis_lakes, morph<1) %>% 
  select(lake,pure,o_e,morph) %>% 
  mutate(plot="low")

high_morph_data<-freq_data %>% filter(lake %in% dis_lakes, morph>1) %>% 
  select(lake,pure,o_e,morph) %>% 
  mutate(plot="high")

morph_plot_data<-bind_rows(all_morph_data,low_morph_data,high_morph_data) %>% 
  mutate(plot=factor(plot, levels=c("low","all","high"),
                     labels = c("Within-Ecotype","All","Between-Ecotype")))

morph_plot<-ggplot(morph_plot_data, aes(x=morph, y=o_e))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.6, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x=expression("Morphological Divergence (" * Delta * "LD" * ")"),
       y="Obs - Exp")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        strip.text = element_text(size=14, face="bold"),
        strip.background = element_rect(fill = "grey95"))+
  facet_wrap(~plot, scales="free_x")

ggsave(morph_plot, filename = "03_figures/assortment/morph-fig.tiff",
       device = "tiff",
       width = 1440, height=500, units = "px",
       dpi=300, scale=2.5)

## plotting fst
lm(o_e~fst, data = filter(freq_data, lake=="GL")) %>% summary()
lm(o_e~fst+morph, data = filter(freq_data, lake=="GL")) %>% summary()
fst_all<-ggplot(filter(freq_data, lake %in% dis_lakes), aes(x=fst, y=o_e))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.6, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x=expression("Genetic Divergence (" * F[ST] * ")"),
       y="Obs - Exp)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14))

## just admixed
ggplot(filter(freq_data, pure==F), aes(x=fst, y=o_e))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.6, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x=expression("Genetic Divergence (" * F[ST] * ")"),
       y="Log Ratio (obs/exp)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=14))

## just trendy lakes
fst_high<-freq_data %>% filter(pure==F, lake %in% c(dis_lakes)) %>% 
ggplot(aes(x=fst, y=o_e))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.6, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x=expression("Genetic Divergence (" * F[ST] * ")"),
       y="Obs - Exp")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14))

lm(o_e~fst, data=filter(freq_data, pure==F, lake %in% dis_lakes)) %>% summary()

fst_plot<-plot_grid(fst_all, fst_high)

ggsave(fst_plot, filename = "03_figures/assortment/fst-fig.tiff",
       device = "tiff",
       width = 1200, height=550, units = "px",
       dpi=300, scale=2.5)

## both model
lm(o_e~fst+morph, data=filter(freq_data, lake %in% dis_lakes)) %>% summary()
lm(o_e~fst+morph, data=filter(freq_data, lake %in% dis_lakes, pure==F)) %>% summary()

## HL
hl_morph<-ggplot(filter(freq_data, lake=="HL"), aes(x=morph, y=o_e))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.5, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x=expression("Morphological Divergence (" * Delta * "LD" * ")"),
       y="Obs - Exp")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14))

hl_fst<-ggplot(filter(freq_data, lake=="HL"), aes(x=fst, y=o_e))+
  geom_smooth(method = "lm", se=F, linetype="dashed", color="grey50", linewidth=0.5)+
  geom_point(aes(fill=pure), shape=21, alpha=0.5, size=3)+
  scale_fill_manual(values = c("grey90","grey25"), labels=c("Admixed","Pure"))+
  labs(x=expression("Genetic Divergence (" * F[ST] * ")"),
       y="Obs - Exp")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14))

hl_plot<-plot_grid(hl_morph, hl_fst)
ggsave(hl_plot, filename = "03_figures/assortment/hl-fig.tiff",
       device = "tiff",
       width = 1200, height=550, units = "px",
       dpi=300, scale=2.5)

lm(o_e~fst+morph, data=filter(freq_data, lake=="HL")) %>% summary()
filter(freq_data, lake=="HL") %>% 
  ggplot(aes(x=fst,y=morph))+
  geom_point()

#### FITNESS ####

fitness_data <- data_f1 %>%
  filter(year == 2) %>%
  rowwise() %>% 
  mutate(
    is_source = any(c_across(all_of(sources)) == 1),
    pure = if_else(is_source, TRUE, FALSE)
  ) %>%
  ungroup() %>%
  mutate(sum_benthic = FG+TL+WB+WT) %>% 
  mutate(eco_pure=case_when(
    sum_benthic==0 ~ T,
    sum_benthic==1 ~ T,
    T ~ F)) %>% 
  select(fish_id, lake, std_length, pure, eco_pure)

ggplot(fitness_data, aes(x=pure, y=std_length))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~lake, scales="free")

ggplot(fitness_data, aes(x=eco_pure, y=std_length))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~lake, scales="free")

#### RANDOM CODE ####

ecotype_exp<-data_f1 %>% filter(lake %in% c("GL","LO")) %>% 
  mutate(limnetic=LG+SL+SR+WK,
         benthic=FG+TL+WB+WT) %>% 
  group_by(lake) %>% 
  summarize(benthic_sum=sum(benthic), limnetic_sum=sum(limnetic), n=n()) %>% 
  mutate(benthic_prop=benthic_sum/(benthic_sum+limnetic_sum),
         limnetic_prop=limnetic_sum/(benthic_sum+limnetic_sum)) %>% 
  mutate(benthic_freq_exp=benthic_prop^2,
         limnetic_freq_exp=limnetic_prop^2,
         hybrid_freq_exp=benthic_prop*limnetic_prop*2) %>% 
  mutate(benthic_count_exp=round(benthic_freq_exp*n),
         limnetic_count_exp=round(limnetic_freq_exp*n),
         hybrid_count_exp=round(hybrid_freq_exp*n)) %>% 
  select(lake,ends_with("count_exp")) %>% 
  pivot_longer(cols=ends_with("exp"), names_to="genotype", values_to = "count_exp") %>% 
  mutate(genotype=str_remove(genotype,"_count_exp"))

ecotype_obs<-data_f1 %>% filter(lake %in% c("GL","LO")) %>% 
  mutate(limnetic=LG+SL+SR+WK,
         benthic=FG+TL+WB+WT) %>% 
  mutate(genotype=case_when(
    limnetic==0 ~ "benthic",
    limnetic==0.5 ~ "hybrid",
    limnetic==1 ~ "limnetic")) %>% 
  group_by(lake,genotype) %>% 
  summarize(count_obs=n())

ecotype_data<-left_join(ecotype_exp,ecotype_obs)

ggplot(ecotype_data, aes(x=sqrt(count_exp),y=sqrt(count_obs),fill=genotype))+
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey")+
  geom_jitter(aes(shape=lake), size=3, alpha=0.75, width=0.1, height=0.1)+
  scale_shape_manual(values = c(21,23))+
  scale_fill_manual(values = c("#006d2c","#BCDECD","#08519c"))+
  labs(x="Expected Count (sqrt)", y="Observed Count (sqrt)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size=12))


