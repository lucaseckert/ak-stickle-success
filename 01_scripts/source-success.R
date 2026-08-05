#### Source Success ####

## Description: quantifying source population success from inferred ancestry data
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## packages 
library(tidyverse)
library(cowplot)
library(patchwork)

## ancestry data
ancestry_data<-read.csv("./02_tidy-data/final-ancestry.csv")

## for factors
sources<-c("LG","SR","WK","SL","WT","WB","FG","TL")
lakes<-c("CL","FL","HL","RL","CC","LL","LP","LO","GL")
limnetic_lakes<-c("CL","FL","HL","RL")
benthic_lakes<-c("CC","LL","LP")
both_lakes<-c("LO","GL")
source_cols<-c("#08519c","#3182bd","#6baed6","#bdd7e7",
               "#fdbe85","#fd8d3c","#e6550d","#a63603")
source_cols_alt<-c("#08519c","#3182bd","#6baed6","#bdd7e7",
                   "#006d2c","#31a354","#74c476","#bae4b3")
## standardize year
ancestry_data<-ancestry_data %>% mutate(year=case_when(
  year==19 ~ 0,
  year==20 ~ 1,
  year==21 ~ 2,
  year==22 & lake!="GL" ~ 3,
  year==22 & lake=="GL" ~ 0.5,
  year==23 & lake=="GL" ~ 1,
  year==24 & lake=="GL" ~ 2
))

## by year, find the proportion of entries where one of the "_round" columns =1
admixture<-ancestry_data %>%
  filter(year %in% c(1,2)) %>% 
  group_by(lake,year) %>%
  summarize(pure = mean(if_any(ends_with("_round"), ~ .x == 1), na.rm = TRUE),
            .groups = "drop") %>% 
  group_by(year) %>%
  summarize(mean(pure))

## source success
round_cols<-c("FG_round","TL_round","WB_round","WT_round","LG_round","SL_round","SR_round","WK_round")
success <- ancestry_data %>%
  group_by(lake, year) %>%
  summarize(across(all_of(round_cols), sum, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = all_of(round_cols), names_to = "source", values_to = "count") %>%
  group_by(lake, year) %>%
  mutate(prop = count / sum(count)) %>%
  mutate(source = factor(substr(source,1,2), levels=sources),
         lake = factor(lake, levels = lakes)) %>% 
  select(lake, year, source, prop)

## mortalities
score_cols<-c("FG_score","TL_score","WB_score","WT_score","LG_score","SL_score","SR_score","WK_score")
morts<-ancestry_data %>% filter(year==0) %>% rowwise() %>% 
  mutate(source = str_remove(score_cols[which.max(c_across(all_of(score_cols)))], "_score")) %>%
  ungroup() %>% select(fish_id,lake,source) %>%
  group_by(lake,source) %>% summarise(dead=n()) %>% 
  mutate(source = factor(source, levels=sources),
         lake = factor(lake, levels=lakes))
morts$dead[is.na(morts$dead)]<-0

## stocking
stocks<-read.csv("./02_tidy-data/other-data/stocking.csv") %>% pivot_longer(cols = FG:WK) %>% 
  select(lake=recipient_lake,source=name,stocked=value) %>% 
  filter(stocked>0) %>% 
  left_join(morts)
stocks$dead[is.na(stocks$dead)]<-0
stocks<-stocks %>% mutate(actual_stocked=stocked-dead) %>% 
  group_by(lake) %>%
  mutate(prop_stocked = stocked / sum(stocked),
         prop_actual_stocked = actual_stocked / sum(actual_stocked)) %>% 
  ungroup() %>% 
  mutate(delta_prop=prop_actual_stocked-prop_stocked) %>% 
  mutate(source = factor(source, levels=sources),
         lake = factor(lake, levels = lakes))

## lethal
lethal<-ancestry_data %>% filter(lake%in%lakes) %>% 
  group_by(lake, year) %>%
  summarize(across(all_of(round_cols), sum, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = all_of(round_cols), names_to = "source", values_to = "count") %>%
  mutate(source=substr(source,1,2)) %>% 
  group_by(lake, year) %>%
  mutate(prop = count / sum(count)) %>%
  mutate(source = factor(source, levels=sources),
         lake = factor(lake, levels=lakes)) %>% 
  select(lake, year, source, prop) %>% 
  pivot_wider(names_from = year, values_from = prop,names_glue = "prop_{year}",values_fill = NA)

## success data to export
ep<-0.01546392
success1<-left_join(stocks,lethal) %>% 
  mutate(logratio_1=log((prop_1+ep)/(prop_actual_stocked+ep)),
         logratio_2=log((prop_2+ep)/(prop_1+ep)),
         logratio_net=log((prop_2+ep)/(prop_actual_stocked+ep))) %>% 
  select(lake,source,logratio_1,logratio_2,logratio_net,prop_actual_stocked,prop_1) %>% 
  pivot_longer(cols = starts_with("logratio"),
               values_to = "logratio",
               names_to = "time",
               names_prefix = "logratio_") %>% 
  mutate(log_init_prop=case_when(
    time %in% c("1","net") ~ log(prop_actual_stocked+ep),
    time == "2" ~ log(prop_1+ep))) %>% 
  select(-starts_with("prop"))

write.csv(success1, file="02_tidy-data/success-data-logratio.csv",
          quote = F, row.names = F)

## props
props<-left_join(stocks,lethal) %>% 
  pivot_longer(cols = starts_with("prop"), names_to = "time", values_to = "prop") %>% 
  mutate(time=factor(time, levels = c("prop_stocked","prop_actual_stocked","prop_1","prop_2","prop_3","prop_0.5"))) %>% 
  mutate(time_numeric=as.numeric(time))

#g_props<-filter(props, lake=="GL", time %in% c("prop_actual_stocked","prop_0.5"))
#saveRDS(g_props, file = "02_tidy-data/g-props.RDS")

## props fig
shape_vals<-rep(c(18,20),4)
props_plot <- props %>%
  filter(time_numeric%in%c(2,3,4), !is.na(prop)) %>%
  mutate(lake = factor(lake, levels = c(benthic_lakes, both_lakes, limnetic_lakes)))
props_fig<-ggplot(props_plot, aes(x=time_numeric,y=prop,color=source))+
  geom_line(linewidth=1.25, alpha=0.7)+
  geom_point(aes(shape = source), size=4.5)+
  scale_shape_manual(values = shape_vals)+
  scale_x_continuous(
    breaks = c(2,3,4),
    labels = c("0", "1", "2"))+
  scale_color_manual(values=source_cols_alt)+
  theme_bw()+
  labs(x="Year", y="Proportion", color="Source\nPopulation", shape="Source\nPopulation")+
  guides(color = guide_legend(
    override.aes = list(
      linetype = 0,
      linewidth = 0)))+
  facet_wrap(~lake, nrow=2, dir="lb")+
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.8),
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 16, hjust=0.5, face="bold"),
        strip.text = element_text(size=14, face="bold"),
        strip.background = element_rect(fill = "grey95"))

ggsave(props_fig,
       filename = "03_figures/props-fig.tiff",
       device = "tiff",
       width = 1440, height=900, units = "px",
       dpi=300, scale=2)

## deltas
deltas<-left_join(stocks,lethal) %>% 
  mutate(delta0=0,
         delta1=(prop_1-prop_actual_stocked),
         delta2=(prop_2-prop_actual_stocked),
         delta3=(prop_3-prop_actual_stocked)) %>% 
  select(-delta_prop) %>% 
  pivot_longer(cols = starts_with("delta"), names_to = "time", values_to = "delta") %>% 
  mutate(time=as.numeric(substr(time,6,6))) %>% 
  mutate(delta_scaled=case_when(
    lake=="LO" ~ delta*(8/4),
    lake=="GL" ~ delta*(7/4),
    TRUE ~ delta
  ))

## this one maybe
deltas_plot <- deltas %>%
  filter(!is.na(delta), time != 3) %>%
  mutate(lake = factor(lake, levels = c(benthic_lakes, both_lakes, limnetic_lakes)))

deltas_fig<- ggplot(deltas_plot, aes(x = time, y = delta, color = source))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey", linewidth=0.75)+
  geom_line(linewidth=1.25, alpha=0.7)+
  geom_point(aes(shape = source), size=4.5)+
  scale_shape_manual(values = shape_vals)+
  scale_x_continuous(
    breaks = c(0, 1, 2),
    labels = c("0", "1", "2"),
    expand = expansion(mult = c(0.1, 0.1)))+
  scale_color_manual(values=source_cols_alt)+
  theme_bw()+
  labs(x="Year", y=expression(Delta~"Proportion"), color="Source\nPopulation", shape="Source\nPopulation")+
  guides(color = guide_legend(
    override.aes = list(
        linetype = 0,
        linewidth = 0)))+
  facet_wrap(~lake, nrow=2, dir="lb")+
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.8),
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 16, hjust=0.5, face="bold"),
        strip.text = element_text(size=14, face="bold"),
        strip.background = element_rect(fill = "grey95"))

ggsave(deltas_fig,
       filename = "03_figures/delta-prop-fig.tiff",
       device = "tiff",
       width = 1440, height=900, units = "px",
       dpi=300, scale=2)

## dans idea
facet_levels <- c(
  "CL", "FL", "HL", "RL",        # row 1
  "CC", "LL", "LP", "BLANK1",    # row 2
  "LO", "GL", "BLANK2", "BLANK3" # row 3
)

deltas_plot3 <- deltas_plot %>%
  mutate(lake_fac = factor(lake, levels = facet_levels)) %>%
  bind_rows(
    tibble(
      lake = NA,
      time = NA,
      delta = NA,
      source = NA,
      lake_fac = factor(
        c("BLANK1", "BLANK2", "BLANK3"),
        levels = facet_levels
      )
    )
  )

rect_df <- data.frame(
  lake_fac = c("BLANK1", "BLANK2", "BLANK3"),
  xmin = 0.5,
  xmax = 1.5,
  ymin = -0.1,
  ymax = 0.1
)

ggplot(deltas_plot3, aes(x = time, y = delta, color = source)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey", linewidth = 0.75) +
  geom_line(linewidth = 1.25, alpha = 0.7) +
  geom_point(aes(shape = source), size = 4.5) +
  scale_shape_manual(values = shape_vals, na.translate = FALSE) +
  scale_color_manual(values = source_cols_alt, na.translate = FALSE) +
  scale_x_continuous(
    breaks = c(0, 1, 2),
    labels = c("0", "1", "2"),
    expand = expansion(mult = c(0.1, 0.1))
  ) +
  facet_wrap(
    ~ lake_fac,
    ncol = 4,
    nrow = 3,
    drop = FALSE,
    labeller = labeller(
      lake_fac = function(x) ifelse(grepl("BLANK", x), "", x)
    )
  ) +
  labs(x="Year", y=expression(Delta~"Proportion"), color="Source\nPopulation", shape="Source\nPopulation") +
  guides(color = guide_legend(
    override.aes = list(
      linetype = 0,
      linewidth = 0
    ),
    nrow = 4
  ),
  shape = guide_legend(nrow = 4)
  )+
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "grey95"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16, hjust = 0.37),
    axis.text = element_text(size = 14),
    legend.position = "inside",
    legend.position.inside = c(0.875, 0.53),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.text = element_text(size=16),
    legend.title = element_text(size = 18, hjust=0.5, face="bold"),
    legend.margin = margin(t = 16, r = 26, b = 20, l = 26)
  )

ggsave(filename = "03_figures/delta-prop-fig-new.tiff",
       device = "tiff",
       width = 1400, height=1000, units = "px",
       dpi=300, scale=2)



#### TIMESPANS ####

## loon
loon<-filter(deltas, lake=="LO", time<4, time>0) %>% 
  mutate(delta = abs(case_when(time==1 ~ prop_1-prop_0,
                               time==2 ~ prop_2-prop_1,
                               time==3 ~ prop_3-prop_2)))

loon_plot<-ggplot(loon, aes(x=time, y=delta, color=source))+
  geom_line(aes(color=source), alpha=0.75, linewidth =1.5)+
  geom_point(aes(shape = source), size=6)+
  scale_shape_manual(values = shape_vals)+
  geom_smooth(method="lm", color="grey50", se=F, linewidth=1, linetype="dashed")+
  scale_color_manual(values=source_cols_alt)+
  scale_x_continuous(breaks=c(1,2,3),
                     labels = c("Year 1","Year 2","Year 3"))+
  labs(x="Year", y=expression("Abs."~Delta~"Proportion"), color="Source\nPopulation", shape="Source\nPopulation")+
  guides(color = guide_legend(
    override.aes = list(
      linetype = 0,
      linewidth = 0)))+
  annotate("text", label="B", x=2.9, y=0.235, size=10)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=18, hjust = 0.65),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=16),
        legend.title = element_text(size = 18, hjust=0.5, face="bold"))

lm(delta~time, data=filter(loon, is.finite(delta))) %>% summary()

## g lake stuff
g_props<-readRDS("02_tidy-data/g-props.RDS") %>%
  select(lake,source,time,prop) %>% 
  pivot_wider(values_from = prop, names_from = time) %>% 
  mutate(delta=prop_0.5-prop_actual_stocked,
         count=prop_0.5*58,
         time=0.5)

g_deltas<-deltas %>% filter(lake=="GL") %>% 
  mutate(delta=prop_2-prop_actual_stocked) %>% 
  select(lake,source,delta) %>% 
  distinct() %>% 
  mutate(time=2) %>% 
  bind_rows(select(g_props,lake,source,delta,time)) %>% 
  bind_rows(
    data.frame(lake="GL", source=g_props$source, time=0, delta=0)
  ) %>% mutate(seg = ifelse(time <= 0.5, "before", "after"))


source_cols_g<-c("#3182bd","#6baed6","#bdd7e7",
                 "#006d2c","#31a354","#74c476","#bae4b3")

g_plot<-ggplot(g_deltas, aes(x=time, y=delta, color=source))+
  geom_vline(xintercept = 0.5, color="grey", linetype="dashed", linewidth=1)+
  geom_line(aes(color=source), alpha=0.75, linewidth =1.5)+
  geom_point(aes(shape = source), size=6)+
  scale_shape_manual(values = shape_vals[-1])+
  scale_color_manual(values=source_cols_alt[-1])+
  scale_x_continuous(breaks = c(0.5,2), labels = c("1 month","2 years"))+
  labs(x="", y=expression(Delta~"Proportion"))+
  annotate("text", label="A", x=0.1, y=0.135, size=10)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=18, hjust=0.8),
        axis.text.y = element_text(size=14))+
  annotate("segment", x = 0.7, xend = 0.8, y = -Inf, yend = -Inf, linewidth=0.1) +
  coord_cartesian(ylim = c(-0.15, 0.15), clip = "off") +
  theme(panel.border = element_rect(color = "black", linewidth = 0.2)) +
  annotate("text", x = 1, y = -0.165, label = "//", size = 6, color="grey20")+
  scale_y_continuous(
    breaks = seq(-0.15, 0.15, by = 0.05), 
    labels = c("-0.15", "-0.10", "-0.05", "0", "0.05", "0.10", "0.15"))

chisq<-chisq.test(x=g_props$count, p=g_props$prop_actual_stocked)
chisq$expected
g_props$p_value_approx <- 2 * pnorm(-abs(chisq$residuals))

time_fig<- plot_grid(g_plot,loon_plot, ncol = 2, rel_widths = c(1,1.3))
ggsave(time_fig, filename = "03_figures/timespans-fig.tiff",
       device = "tiff",
       width = 1440, height=550, units = "px",
       dpi=300, scale=2.5)


#### OLD STUFF #####

ep<-0.01546392
deltas2<-left_join(stocks,lethal) %>% 
  mutate(delta0=0,
         delta1=log((prop_1+ep)/(prop_actual_stocked+ep)),
         delta2=log((prop_2+ep)/(prop_actual_stocked+ep)),
         delta3=log((prop_3+ep)/(prop_actual_stocked+ep))) %>% 
  select(-delta_prop) %>% 
  pivot_longer(cols = starts_with("delta"), names_to = "time", values_to = "delta_logratio") %>% 
  mutate(time=as.numeric(substr(time,6,6)))

## selection for marqus
loon<-props %>% filter(lake=="LO",time%in% c("prop_actual_stocked","prop_3")) %>% 
  mutate(ecotype=factor(if_else(source%in%c("FG","TL","WB","WT"),"benthic","limnetic"))) %>% 
  group_by(time,ecotype) %>% 
  summarize(ecotype_prop=sum(prop))
s<-((loon[4,3]/loon[2,3])^(1/3))-1

## this one maybe
deltas$lake<-factor(deltas$lake, levels = c(benthic_lakes,both_lakes,limnetic_lakes))
ggplot(filter(deltas, !is.na(delta)), aes(x=time,y=delta,color=source))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey", linewidth=1)+
  geom_line(linewidth=1.5)+
  geom_point(size=3)+
  scale_x_continuous(
    breaks = c(0, 1, 2,3),
    labels = c("0", "1", "2","3"))+
  scale_color_manual(values=source_cols_alt)+
  theme_bw()+
  labs(x="Year", y="Change in Frequency", color="Source\nPopulation")+
  facet_wrap(~lake, nrow=2, dir="lb")+
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.8),
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 16, hjust=0.5, face="bold"),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_rect(fill = "grey95"))

ggsave(filename = "03_figures/delta-prop-fig.tiff",
       device = "tiff",
       width = 1440, height=900, units = "px",
       dpi=300, scale=2)

## log ratio plot
deltas2$lake<-factor(deltas2$lake, levels = c(benthic_lakes,both_lakes,limnetic_lakes))
ggplot(filter(deltas2, !is.na(delta_logratio)), aes(x=time,y=delta_logratio,color=source))+
  #geom_hline(yintercept = 0, linetype="dashed", color="grey", linewidth=1)+
  geom_line(linewidth=1.5)+
  geom_point(size=3)+
  scale_x_continuous(
    breaks = c(0, 1, 2,3),
    labels = c("0", "1", "2","3"))+
  scale_color_manual(values=source_cols_alt)+
  theme_bw()+
  labs(x="Year", y="Change in Frequency", color="Source\nPopulation")+
  facet_wrap(~lake, nrow=2, dir="lb")+
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.8),
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 16, hjust=0.5, face="bold"),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_rect(fill = "grey95"))

delta_means<-deltas %>% group_by(lake,time) %>% 
  summarize(delta_mean=mean(delta,na.rm=T),
            delta_sd=sd(delta,na.rm=T))

deltas<-left_join(deltas,delta_means) %>% 
  mutate(z_score=(delta-delta_mean)/delta_sd)

deltas$z_score[deltas$time==0]<-0

## THIS IS THE ONE ##
ggplot(filter(deltas, !is.na(z_score)), aes(x=time,y=z_score,color=source))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey", linewidth=1)+
  geom_line(linewidth=1.5)+
  geom_point(size=3)+
  scale_x_continuous(
    breaks = c(0, 1, 2,3),
    labels = c("0", "1", "2","3"))+
  scale_color_manual(values=source_cols_alt)+
  theme_bw()+
  labs(x="Year", y="Change in Success", color="Source\nPopulation")+
  facet_wrap(~lake, nrow=2, dir="lb")+
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.8),
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 16, hjust=0.5, face="bold"),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_rect(fill = "grey95"))
ggsave(filename = "03_figures/delta-z-fig.tiff",
       device = "tiff",
       width = 1440, height=900, units = "px",
       dpi=300, scale=2)

deltas2<-left_join(stocks,lethal) %>% 
  mutate(delta_prop_0=delta_prop,
         delta_prop_1=prop_20-prop_actual_stocked,
         delta_prop_2=prop_21-prop_20,
         delta_prop_3=prop_22-prop_21) %>% 
  select(-delta_prop) %>% 
  pivot_longer(cols = starts_with("delta"), names_to = "time", values_to = "delta") %>% 
  mutate(time=as.numeric(substr(time,12,12)))

ggplot(filter(deltas2, time>0, time<3), aes(x=source, y=delta, fill=factor(time), color=source))+
  geom_boxplot()+
  scale_color_manual(values=source_cols)

deltas_sd<-deltas %>% mutate(source_ecotype = factor(case_when(
  lake %in% limnetic_lakes ~ "limnetic",
  lake %in% benthic_lakes ~ "benthic",
  lake %in% both_lakes ~ "both"))) %>%
  filter(!is.na(delta)) %>% 
  group_by(source_ecotype,time) %>% 
  summarize(delta_sd=sd(delta))

ggplot(filter(deltas, time%in%c(1,2)), aes(x=source, y=delta, fill=factor(time)))+
  geom_boxplot()

## SUCCESS DATA
deltas2<-left_join(stocks,lethal) %>% 
  mutate(delta_1=(prop_1-prop_actual_stocked),
         delta_2=(prop_2-prop_1),
         delta_net=(prop_2-prop_actual_stocked)) %>% 
  select(lake,source,delta_1,delta_2,delta_net,
         prop_1, prop_2, prop_actual_stocked) %>% 
  pivot_longer(cols = starts_with("delta"), names_to = "time", values_to = "delta") %>% 
  mutate(delta_norm=case_when(
    lake=="LO" ~ delta*2,
    lake=="GL" ~ delta*(7/4),
    TRUE ~ delta)) %>% 
  mutate(initial_prop=case_when(
    time=="delta_1" ~ prop_actual_stocked,
    time=="delta_2" ~ prop_1,
    time== "delta_net" ~ prop_actual_stocked)) %>% 
  mutate(initial_prop_norm=case_when(
    lake=="LO" ~ initial_prop*2,
    lake=="GL" ~ initial_prop*(7/4),
    TRUE ~ initial_prop)) %>% 
  select(lake, source, time, delta, delta_norm, initial_prop_norm)

write.csv(deltas2, file="02_tidy-data/success-data.csv",
          quote = F, row.names = F)

## multinomial
success2<-left_join(stocks,lethal) %>% 
  mutate(prop_0=prop_actual_stocked,
         count_1=round(prop_1*1000),
         count_2=round(prop_2*1000),
         count_net=round(prop_2*1000),
         n_total=1000,
         introduced=1) %>% 
  select(lake,source,introduced,n_total,count_1,count_2,count_net,prop_0,prop_1)

df_full<-expand_grid(lake = lakes, source = sources) %>% 
  left_join(success2, by = c("lake", "source")) %>% 
  mutate(introduced = ifelse(is.na(prop_0), 0, 1),
         across(where(is.numeric), ~replace_na(., 0))) %>% 
  pivot_longer(cols = starts_with("count"),
               values_to = "count",
               names_to = "time",
               names_prefix = "count_") %>% 
  mutate(init_prop=case_when(
    time %in% c("1","net") ~ prop_0,
    time == "2" ~ prop_1)) %>% 
  select(-starts_with("prop"))

write.csv(df_full, file="02_tidy-data/success-data-multinomial.csv",
          quote = F, row.names = F)

## freq dep
deltas2 %>% select(-initial_prop_norm) %>% 
  pivot_wider(names_from = time, values_from = delta_norm) %>% 
  ggplot(aes(x=delta_1, y=delta_2, color=source, fill=source))+
  geom_point(aes(fill=source),size=3, shape=21, alpha=0.75)+
  geom_smooth(aes(fill=source), method="lm", se=T, alpha=0.15)+
  scale_color_manual(values=source_cols_alt)+
  scale_fill_manual(values=source_cols_alt)+
  labs(x="Change in Proportion in Year 1", y="Change in Propotion in Year 2")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))

freq_dep_plot<-deltas2 %>% filter(time %in% c("delta_2")) %>% 
ggplot(aes(x=initial_prop_norm, y=delta_norm, color=source, fill=source))+
  geom_point(aes(fill=source),size=3, shape=21, alpha=0.75)+
  geom_smooth(aes(fill=source), method="lm", se=T, alpha=0.15)+
  scale_color_manual(values=source_cols_alt)+
  scale_fill_manual(values=source_cols_alt)+
  labs(x="Proportion after Year 1", y="Change in Propotion in Year 2")+
  annotate("text", label="B", x=0.025, y=0.35, size=12)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))

saveRDS(freq_dep_plot, file = "03_figures/figure-parts/freq-dep-plot.RDS")

#### OLD CODE ####

## environmental data
env<-read.csv("./02_tidy-data/other-data/env.csv")
pca_full<-prcomp(select(env,-lake,-lat,-long,-lake.type,-ecotype,-region))
pc_scores<-pca_full$x[,1:5]
rownames(pc_scores)<-env$lake
dist_mat<-dist(pc_scores, method = "euclidian") %>% as.matrix()

env_dist<-c()
for (i in 1:nrow(deltas2)) {
 env_dist[i]<-dist_mat[as.character(deltas2$lake[i]),as.character(deltas2$source[i])]
}
deltas2$env_dist<-env_dist

deltas2<-deltas2 %>% 
  mutate(env_dist_scaled=scale(env_dist),
         initial_prop_scaled=scale(initial_prop))
  
ggplot(pc_scores, aes(x=PC1,y=PC2, label=rownames(pc_scores)))+
  geom_text()

## regression
lmm<-lmer(z_score ~ source + env_dist_scaled + initial_prop_scaled + (1 | lake), data = filter(deltas2, time==1))
summary(lmm)

lmm2<-lmer(z_score ~ source + env_dist_scaled + initial_prop_scaled + (1 | lake), data = filter(deltas2, time==2))
summary(lmm2)

## other env data
env2<-read.csv("./02_tidy-data/other-data/env-2022.csv") %>% 
  filter(!is.na(doc)) %>% 
  left_join(select(env,lake,area,z,calc,litt.prop)) %>% 
  bind_rows(filter(env,lake %in% sources)) %>% 
  select(-lat,-long,-lake.type,-ecotype,-region,-secchi,-zoo,-daph,-daph.prop,-invert,-gam,-gam.prop,-chiro,-chiro.prop)

pca2<-prcomp(select(env2,-lake))
pc_scores2<-pca2$x[,1:3]
rownames(pc_scores2)<-env2$lake
dist_mat2<-dist(pc_scores2, method = "euclidian") %>% as.matrix()

env_dist2<-c()
for (i in 1:nrow(deltas2)) {
  env_dist2[i]<-dist_mat2[as.character(deltas2$lake[i]),as.character(deltas2$source[i])]
}
deltas2$env_dist_2<-scale(env_dist2)

## regression
lmm3<-lmer(z_score ~ source + env_dist_2 + initial_prop_scaled + (1 | lake), data = filter(deltas2, time==1))
summary(lmm)

lmm4<-lmer(z_score ~ source + env_dist_2 + initial_prop_scaled + (1 | lake), data = filter(deltas2, time==2))
summary(lmm2)

## freq dep
ggplot(filter(deltas2, time==2), aes(x=initial_prop_scaled, y=z_score))+
  geom_point()+
  geom_smooth(method="lm")


## population level stuff
pop_data<-data.frame(source=sources, effect_size=lmm@beta[1:8]) %>% 
  left_join(select(mort_means,source,overall_rate))

ggplot(pop_data, aes(x=overall_rate,y=effect_size,label=source))+
  geom_text()

lethal<-read.csv("./02_tidy-data/other-data/lethal-2021.csv") %>% 
  mutate(source=substr(Fish_ID,1,2)) %>% 
  filter(source %in% sources, !is.na(Mass_g)) %>% 
  mutate(source=factor(source, levels = sources)) %>% 
  group_by(source) %>% 
  summarize(mean_sl=mean(Std_length_mm, na.rm=T),
            mean_mass=mean(Mass_g, na.rm=T),
            parasite_freq=sum(S_solidus_N>0)/n())

pop_data<-left_join(pop_data,lethal)

ggplot(pop_data, aes(x=mean_sl,y=effect_size,label=source))+
  geom_text()
ggplot(pop_data, aes(x=mean_mass,y=effect_size,label=source))+
  geom_text()
ggplot(pop_data, aes(x=parasite_freq,y=effect_size,label=source))+
  geom_text()

lm1<-lm(effect_size~scale(overall_rate)+scale(mean_sl)+scale(parasite_freq),
        data=pop_data)
summary(lm1)

## loon
loon<-filter(deltas, lake=="LO", time<4, time>0) %>% 
  mutate(delta = abs(case_when(time==1 ~ prop_1-prop_0,
                           time==2 ~ prop_2-prop_1,
                           time==3 ~ prop_3-prop_2)))

loon_plot<-ggplot(loon, aes(x=time,y=delta))+
  geom_line(aes(color=source), alpha=0.5, linewidth =2)+
  geom_point(aes(fill=source), size=5, shape=21, alpha=0.75)+
  geom_smooth(method="lm", color="grey40", se=F, linewidth=1, linetype="dashed")+
  scale_fill_manual(values=source_cols_alt)+
  scale_color_manual(values=source_cols_alt)+
  scale_x_continuous(breaks=c(1,2,3))+
  labs(x="Year", y="Abs. Change in Frequency")+
  annotate("text", label="B", x=2.9, y=0.23, size=10)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=14))

lm(delta~time, data=filter(loon, is.finite(delta))) %>% summary()

## g lake stuff
g_props<-readRDS("02_tidy-data/g-props.RDS") %>%
  select(lake,source,time,prop) %>% 
  pivot_wider(values_from = prop, names_from = time) %>% 
  mutate(delta=prop_0.5-prop_actual_stocked,
         count=prop_0.5*58,
         time=0.5)

g_deltas<-deltas %>% filter(lake=="GL") %>% 
  mutate(delta=prop_2-prop_actual_stocked) %>% 
  select(lake,source,delta) %>% 
  distinct() %>% 
  mutate(time=2) %>% 
  bind_rows(select(g_props,lake,source,delta,time)) %>% 
  bind_rows(
    data.frame(lake="GL", source=g_props$source, time=0, delta=0)
  ) %>% mutate(seg = ifelse(time <= 0.5, "before", "after"))
  

source_cols_g<-c("#08519c","#3182bd","#bdd7e7",
                 "#bae4b3","#74c476","#31a354","#006d2c")

g_plot<-ggplot(g_deltas, aes(x=time, y=delta))+
  #geom_hline(yintercept = 0, color="grey", linetype="dashed", linewidth=1)+
  geom_vline(xintercept = 0.5, color="grey", linetype="dashed", linewidth=1)+
  geom_line(data=filter(g_deltas, time<2), aes(color = source), linewidth = 2, alpha = 0.75) +
  geom_line(data=filter(g_deltas, time>0), aes(color = source), linewidth = 2, alpha = 0.5) +
  geom_point(aes(fill=source), size=5, shape=21, alpha=0.75)+
  scale_color_manual(values=source_cols_g)+
  scale_fill_manual(values=source_cols_g)+
  scale_x_continuous(breaks = c(0.5,2), labels = c("1 month","2 years"))+
  labs(x="", y="Change in Frequency")+
  annotate("text", label="A", x=0.1, y=0.11, size=10)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=18, hjust=0.8),
        axis.text.y = element_text(size=14))


g_plot<-ggplot(g_props, aes(x=source, y=delta, fill=source))+
  geom_bar(stat="identity", color="black", linewidth =0.5, alpha=0.75)+
  geom_hline(yintercept = 0,color="black", linewidth=0.5)+
  scale_fill_manual(values = source_cols_g)+
  labs(x="Source Population", y="Change in Frequency")+
  annotate("text", label="A", x=1, y=0.12, size=10)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=14))

chisq<-chisq.test(x=g_props$count, p=g_props$prop_actual_stocked)

chisq$expected
g_props$p_value_approx <- 2 * pnorm(-abs(chisq$residuals))

fig<- plot_grid(g_plot,loon_plot, ncol = 2)
ggsave(fig, filename = "03_figures/timespans-fig.tiff",
       device = "tiff",
       width = 1440, height=550, units = "px",
       dpi=300, scale=2.5)

deltas3<- deltas %>%
  filter(time %in% c(1,2)) %>% 
  mutate(abs_delta=abs(case_when(time==1 ~ prop_1-prop_actual_stocked,
                                 time==2 ~ prop_2-prop_1))) %>% 
  select(lake,source,time,abs_delta)

ggplot(deltas3, aes(x=time, y=abs_delta))+
  geom_line(aes(color=source, group = interaction(lake,source)), alpha=0.3, linewidth =1)+
  geom_jitter(aes(fill=source), width = 0, shape=21, size=4, alpha=0.75)+
  #geom_boxplot(outliers=F, alpha=0.5, fill="grey90", staplewidth = 0.1)+
  geom_smooth(method = "lm", linewidth=1.25, color="grey50", linetype="dashed", se=F)+
  scale_fill_manual(values=source_cols_alt)+
  scale_color_manual(values=source_cols_alt)+
  scale_x_continuous(breaks=c(1,2), labels = c("Year 1","Year 2"), limits=c(0.85,2.15))+
  labs(y="Abs. Change in Frequency")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank())
ggsave(filename = "03_figures/timespans-sup-fig.tiff",
       device = "tiff",
       width = 800, height=600, units = "px",
       dpi=300, scale=2.5)

deltas4<-deltas3 %>% mutate(time=paste0("time_",time)) %>% 
  pivot_wider(names_from = time, values_from = abs_delta) %>% 
  mutate(delta= time_2 - time_1) %>% 
  mutate(check=delta<0)
mean(deltas4$delta, na.rm=T)
sum(deltas4$check, na.rm=T)/nrow(filter(deltas4,!is.na(check)))

lm(abs_delta~time+lake*source, data = deltas3) %>% summary()
t.test(abs_delta~time, data = deltas3)

ggplot(success, aes(x=lake, y=prop, fill=source))+
  geom_bar(position = "stack", stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = source_cols)+
  xlab("")+
  ylab("Proportional Contribution")+
  guides(fill=guide_legend(title="Source Population:"))+
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        axis.text.y = element_text(size=18)) +
  facet_grid(~year)
