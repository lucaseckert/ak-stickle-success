#### Environmental Correlations ####

## packages 
library(tidyverse)
library(MuMIn)
library(GGally)
library(cowplot)
library(ggrepel)

## for plots
sources<-c("SR","WK","LG","SL","WB","FG","WT","TL")
source_cols_alt<-c("#08519c","#3182bd","#6baed6","#bdd7e7",
                   "#006d2c","#31a354","#74c476","#bae4b3")

## data
data<-read.csv("02_tidy-data/other-data/morph-ld.csv") %>% 
  left_join(rename(read.csv("02_tidy-data/other-data/env.csv"),source=lake)) %>% 
  mutate(source=factor(source, levels=sources),
         invert_prop=invert/zoo)

## explore
ggpairs(select(data,ld,tp,chla,sp.cond,ph,calc,litt.prop,zoo,invert,invert_prop))

## ggplot the correlation bewtween litt.prop and calc
ggplot(data, aes(x=litt.prop, y=calc))+
  geom_point(aes(fill=source), shape=21, size=5)+
  scale_fill_manual(values=source_cols_alt)+
  labs(x="Relative Littoral Area", y="Calcium Concentration")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=14))

shape_vals<-rep(c(18,20),4)

## rla
lm_rla<-lm(ld~litt.prop, data = data)
rla_anno<-sprintf("adj. r² = %.3f\np = %.3g", summary(lm_rla)$adj.r.squared, summary(lm_rla)$coefficients[2,4])
rla_plot<-ggplot(data, aes(x=litt.prop, y=ld))+
  geom_smooth(method="lm", color="grey", linetype="dashed", se=F, linewidth=1)+
  geom_point(aes(color=source, shape=source), size=7)+
  scale_color_manual(values=source_cols_alt)+
  scale_shape_manual(values=shape_vals)+
  labs(x="Relative Littoral Area", y="Body Shape LD", color="Source\nPopulation", shape="Source\nPopulation")+
  annotate("text", label=rla_anno, x=0.8, y=0.5, size=6)+
  annotate("text", label="A", x=0.9, y=2, size=12)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=14))

## calc
lm_calc<-lm(ld~calc, data = data)
calc_anno<-sprintf("adj. r² = %.3f\np = %.3g", summary(lm_calc)$adj.r.squared, summary(lm_calc)$coefficients[2,4])
calc_plot<-ggplot(data, aes(x=calc, y=ld))+
  geom_smooth(method="lm", color="grey", linetype="dashed", se=F, linewidth=1)+
  geom_point(aes(color=source, shape=source), size=7)+
  scale_color_manual(values=source_cols_alt)+
  scale_shape_manual(values=shape_vals)+
  labs(x="Calcium Concentration", y="Body Shape LD", color="Source\nPopulation", shape="Source\nPopulation")+
  annotate("text", label=calc_anno, x=25, y=0.45, size=6)+
  annotate("text", label="B", x=28, y=2.15, size=12)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=16),
        legend.title = element_text(size = 18, hjust=0.5, face="bold"))

## saving
fig<- plot_grid(rla_plot,calc_plot, ncol = 2, rel_widths = c(1,1.3))
ggsave(fig, filename = "03_figures/env-correlations-fig.tiff",
       device = "tiff",
       width = 1440, height=600, units = "px",
       dpi=300, scale=2.5)

## modelling
global_1<-lm(ld~calc+litt.prop+zoo+daph+invert+gam+chiro,
             data = data,
             na.action = na.fail)
global_2<-lm(ld~calc+litt.prop+tp+zoo+daph.prop+invert+chiro.prop,
             data = data,
             na.action = na.fail)

dredge_1<-dredge(global_1, m.lim=c(2,4))
dredge_2<-dredge(global_2, m.lim=c(2,4))

view(dredge_1)
view(dredge_2)

filter(dredge_1, !is.na(calc), !is.na(litt.prop)) %>% view()
filter(dredge_2, !is.na(calc), !is.na(litt.prop)) %>% view()

## PCA plot
env<-read.csv("02_tidy-data/env-data.csv")
var5<-c("rla","calc","max_depth","temp","do","conductivity","ph","secchi_depth","tp","tn","doc","chla","macro_total_pre","zoo_total_pre")
pca_full<-prcomp(select(env,all_of(var5)), center = T, scale. = T)
env<-cbind(env,pca_full$x)

# Extract loadings
loadings <- as.data.frame(pca_full$rotation[, 1:2]) %>%
  mutate(var = rownames(.))

# Scale loadings to fit on same plot (arrows are often small compared to PCA scores)
loadings <- loadings %>%
  mutate(PC1 = PC1 * 12,   # adjust multiplier as needed
         PC2 = PC2 * 12,
         var = str_remove(var, "_pre"))

pc_lakes<-ggplot(env, aes(x = PC1, y = PC2, fill = lake_type)) +
  #geom_segment(
    #data = loadings, inherit.aes = FALSE,
    #aes(x = 0, y = 0, xend = PC1, yend = PC2),
    #arrow = arrow(length = unit(0.25, "cm")),
    #color = "grey80", linewidth = 0.5) +
  #geom_text(
    #data = loadings, inherit.aes = FALSE,
    #aes(x = PC1 * 1.2, y = PC2 * 1.2, label = var),
    #size = 4, fontface = "italic", color="grey50") +
  stat_ellipse(aes(color = lake_type), linewidth = 1, linetype = "dashed") +
  geom_point(shape = 21, size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = lake),
    vjust = -1, size = 4) +
  scale_fill_manual(values = c("grey5", "white")) +
  scale_color_manual(values = c("grey5", "grey75")) +
  labs(x="PC1", y="PC2") +
  lims(x=c(-8,8), y=c(-6,6))+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
    legend.position = "inside",
    axis.title = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_blank(),
    legend.position.inside = c(0.2,0.15),
    legend.text = element_text(size=14))

pc_loadings<-ggplot(env, aes(x = PC1, y = PC2)) +
  geom_segment(
   data = loadings, inherit.aes = FALSE,
   aes(x = 0, y = 0, xend = PC1, yend = PC2),
   arrow = arrow(length = unit(0.25, "cm")),
   color = "grey", linewidth = 0.5) +
  geom_text_repel(
   data = loadings, inherit.aes = FALSE,
   aes(x = PC1 * 1.1, y = PC2 * 1.1, label = var),
   size = 4, fontface = "italic", color="black") +
  labs(x="PC1", y="PC2") +
  lims(x=c(-8,8), y=c(-6,6))+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.position.inside = c(0.85,0.2),
        legend.text = element_text(size=16))

pc_plot<-plot_grid(pc_lakes, pc_loadings, ncol=2, rel_widths = c(1,1))

ggsave(pc_plot,
       filename = "03_figures/env-pca.tiff",
       device = "tiff",
       width = 1200, height=600, units = "px",
       dpi=300, scale=2.5)

## new data

## data
data<-read.csv("02_tidy-data/other-data/morph-ld.csv") %>% 
  left_join(rename(read.csv("02_tidy-data/env-data.csv"),source=lake)) %>% 
  mutate(source=factor(source, levels=sources))

## forward selection?
null_model<-lm(ld~1, data = data)
full_formula <-ld~rla+calc+max_depth+temp+do+conductivity+ph+secchi_depth+tp+tn+doc+chla+macro_total_pre+zoo_total_pre+macro_total_post+zoo_total_post
forward_fit <- MASS::stepAIC(
  object = null_model,
  scope = list(upper = full_formula, lower = ~ 1),
  direction = "forward",
  trace = TRUE
)

scale_data<-data %>% mutate(across(is.numeric, scale))

global_1<-lm(ld~rla+calc+max_depth+temp+do+conductivity+ph+secchi_depth+tp+tn+doc+chla+macro_total_pre+zoo_total_pre+macro_total_post+zoo_total_post,
             data = scale_data,
             na.action = na.fail)

dredge_1<-MuMIn::dredge(global_1, m.lim=c(2,4))

d_subset<-filter(dredge_1, !is.na(calc), !is.na(rla))
d1 <- d_subset[order(d_subset$delta), ]

first_non_na <- sapply(d1, function(col) {
  out <- which(!is.na(col))
  if (length(out) == 0) return(NA)  # column is all NA
  return(out[1] - 1)  # how many rows before first non-NA
})

