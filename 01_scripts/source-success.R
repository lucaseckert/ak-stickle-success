#### SOURCE SUCCESS ####

## Description: quantifying source population success from inferred ancestry data
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## constants
source(here::here("02_tidy-data/constants.R"))

## ancestry data
ancestry_data<-read.csv("./02_tidy-data/final-ancestry.csv")

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

## summary
ancestry_data %>% count(year,lake)

## percent admixture
admixture<-ancestry_data %>%
  filter(year %in% c(1,2)) %>% 
  group_by(lake,year) %>%
  summarize(pure = mean(if_any(ends_with("_round"), ~ .x == 1), na.rm = TRUE),
            .groups = "drop") %>% 
  group_by(year) %>%
  summarize(mean(pure))

## source success proportions
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

## mortality data
score_cols<-c("FG_score","TL_score","WB_score","WT_score","LG_score","SL_score","SR_score","WK_score")
morts<-ancestry_data %>% filter(year==0) %>% rowwise() %>% 
  mutate(source = str_remove(score_cols[which.max(c_across(all_of(score_cols)))], "_score")) %>%
  ungroup() %>% select(fish_id,lake,source) %>%
  group_by(lake,source) %>% summarise(dead=n()) %>% 
  mutate(source = factor(source, levels=sources),
         lake = factor(lake, levels=lakes))
morts$dead[is.na(morts$dead)]<-0

## stocking data
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

## get counts
count_data<-stocks %>% 
  select(lake, source, count_t0=actual_stocked, prop_t0=prop_actual_stocked)

## lethal data
lethal<-ancestry_data %>% filter(lake%in%lakes) %>% 
  group_by(lake, year) %>%
  summarize(across(all_of(round_cols), sum, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(cols = all_of(round_cols), names_to = "source", values_to = "count") %>%
  mutate(source=substr(source,1,2))

lethal_totals<-lethal %>% 
  filter(year %in% c(1,2)) %>%
  group_by(lake, year) %>% 
  summarize(N = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = year, values_from = N, names_prefix = "N_t")

lethal_counts<-lethal %>% group_by(lake, year) %>% 
  filter(year %in% c(1,2)) %>%
  mutate(prop = count / sum(count)) %>%
  mutate(source = factor(source, levels=sources),
         lake = factor(lake, levels=lakes)) %>% 
  select(lake, year, source, count, prop) %>% 
  pivot_wider(names_from = year, values_from = c(count, prop), names_glue = "{.value}_t{year}", values_fill = NA)

lethal<-lethal %>% group_by(lake, year) %>% 
  mutate(prop = count / sum(count)) %>%
  mutate(source = factor(source, levels=sources),
         lake = factor(lake, levels=lakes)) %>% 
  select(lake, year, source, prop) %>% 
  pivot_wider(names_from = year, values_from = prop,names_glue = "prop_{year}",values_fill = NA)

count_data<-count_data %>% left_join(lethal_counts) %>% 
  left_join(lethal_totals)

## count data for models
write.csv(count_data, file="02_tidy-data/source-success-counts.csv",
          quote = F, row.names = F)

## final proportions
props<-left_join(stocks,lethal) %>% 
  pivot_longer(cols = starts_with("prop"), names_to = "time", values_to = "prop") %>% 
  mutate(time=factor(time, levels = c("prop_stocked","prop_actual_stocked","prop_1","prop_2","prop_3","prop_0.5"))) %>% 
  mutate(time_numeric=as.numeric(time))

## save g proportions for later
g_props<-filter(props, lake=="GL", time %in% c("prop_actual_stocked","prop_0.5"))
saveRDS(g_props, file = "02_tidy-data/g-props.RDS")

## props fig
shape_vals<-rep(c(21,22,23,24),2)
props_plot <- props %>%
  filter(time_numeric%in%c(2,3,4), !is.na(prop)) %>%
  mutate(lake = factor(lake, levels = c(benthic_lakes, both_lakes, limnetic_lakes)))

props_fig<-ggplot(props_plot, aes(x=time_numeric,y=prop,color=source))+
  geom_line(linewidth=1.25, alpha=0.7)+
  geom_point(aes(shape = source, fill = source), size=4)+
  scale_shape_manual(values = shape_vals)+
  scale_x_continuous(
    breaks = c(2,3,4),
    labels = c("0", "1", "2"))+
  scale_color_manual(values=source_cols_alt)+
  scale_fill_manual(values=source_cols_alt)+
  theme_bw()+
  labs(x="Year", y="Proportion", color="Source\nPopulation", shape="Source\nPopulation")+
  guides(
    color = guide_legend(
      override.aes = list(
        shape = shape_vals,
        fill = source_cols_alt,
        linetype = 0,
        linewidth = 0,
        alpha = 1)),
    shape = "none",
    fill = "none") +
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

## data to save for later
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

## this one maybe
deltas_plot <- deltas %>%
  filter(!is.na(delta), time != 3) %>%
  mutate(lake = factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes)))

design1 <- "
 ABCD#
 EFGHI
"

deltas_fig<- ggplot(deltas_plot, aes(x = time, y = delta_scaled, color = source))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey", linewidth=0.75)+
  geom_line(linewidth=1.25, alpha=0.7)+
  geom_point(aes(shape = source, fill = source), size=4)+
  scale_shape_manual(values = shape_vals,labels=paste0("S-",sources))+
  scale_x_continuous(
    breaks = c(0, 1, 2),
    labels = c("0", "1", "2"),
    expand = expansion(mult = c(0.1, 0.1)))+
  scale_color_manual(values=source_cols_alt, labels=paste0("S-",sources))+
  scale_fill_manual(values=source_cols_alt, labels=paste0("S-",sources))+
  theme_bw()+
  labs(x="Year", y=expression(Delta~"Proportion"), color="Source\nPopulation", shape="Source\nPopulation")+
  guides(
    color = guide_legend(
      override.aes = list(
        shape = shape_vals,
        fill = source_cols_alt,
        linetype = 0,
        linewidth = 0,
        alpha = 1)),
    shape = "none",
    fill = "none") +
  ggtitle("Recipient Lake")+
  facet_manual(
    ~factor(lake, labels = paste0("R-", c(limnetic_lakes, benthic_lakes, both_lakes))),
    design=design1)+
  #strip = strip_themed(
  #background_x = elem_list_rect(fill = c(rep(limnetic_col,4),rep(benthic_col,3),rep(both_col,2)))))+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size=16, hjust=0.5, face="bold"),
        legend.position = "inside",
        legend.position.inside = c(0.9,0.825),
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(size=14, face="bold"))

ggsave(deltas_fig,
       filename = "03_figures/delta-prop-fig.tiff",
       device = "tiff",
       width = 1440, height=900, units = "px",
       dpi=300, scale=2)

## other way
design2 <- "
 ABCD
 EFGH
"

## by source population
source_fig<-deltas_plot %>% 
  mutate(ecotype=case_when(
    lake %in% limnetic_lakes ~ "Limnetic",
    lake %in% benthic_lakes ~ "Benthic",
    lake %in% both_lakes ~ "Both"
  )) %>%
  ggplot(aes(x = time, y = delta_scaled, group = lake, color=source))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey", linewidth=0.75)+
  geom_line(linewidth=1.25, alpha=0.7)+
  geom_point(aes(shape=source, fill=source), size=4)+
  scale_shape_manual(values = shape_vals,labels=paste0("S-",sources))+
  scale_fill_manual(values=source_cols_alt, labels=paste0("S-",sources))+
  scale_x_continuous(
    breaks = c(0, 1, 2),
    labels = c("0", "1", "2"),
    expand = expansion(mult = c(0.1, 0.1)))+
  scale_color_manual(values=source_cols_alt)+
  theme_bw()+
  labs(x="Year", y=expression(Delta~"Proportion"), color="Source\nPopulation", shape="Source\nPopulation")+
  guides(
    color = guide_legend(
      override.aes = list(
        shape = shape_vals,
        fill = source_cols_alt,
        linetype = 0,
        linewidth = 0,
        alpha = 1)),
    shape = "none",
    fill = "none") +
  ggtitle("Source Population")+
  facet_manual(
    ~factor(source, labels = paste0("S-", sources)),
    design = design2)+
  #strip = strip_themed(
  #background_x = elem_list_rect(fill = source_cols_alt),
  #text_x = elem_list_text(color = c(rep("white",2),rep("black",4),rep("white",2)))))+
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size=16, hjust=0.5, face="bold"),
        legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size=14),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(size=14, face="bold"))

## final figure
fig2<-cowplot::plot_grid(deltas_fig, source_fig,
                         nrow = 2,
                         rel_heights = c(1,0.9))
## save
save_fig(
  plot = fig2,
  filename = "03_figures/fig2-new.tiff",
  width = 0.75, height=0.75)

#### TIMESPANS ####

## loon
loon<-filter(deltas, lake=="LO", time<4, time>0) %>% 
  mutate(delta = abs(case_when(time==1 ~ prop_1-prop_0,
                               time==2 ~ prop_2-prop_1,
                               time==3 ~ prop_3-prop_2)))

loon_plot<-ggplot(loon, aes(x=time, y=delta, color=source))+
  geom_line(aes(color=source), alpha=0.75, linewidth =1.5)+
  geom_point(aes(shape = source, fill=source), size=5.5)+
  scale_shape_manual(values = shape_vals, labels = paste0("S-",sources))+
  scale_fill_manual(values=source_cols_alt, labels = paste0("S-",sources))+
  geom_smooth(method="lm", color="grey50", se=F, linewidth=1, linetype="dashed")+
  scale_color_manual(values=source_cols_alt, labels = paste0("S-",sources))+
  scale_x_continuous(breaks=c(1,2,3),
                     labels = c("Year 1","Year 2","Year 3"))+
  labs(x="Year", y=expression("Abs."~Delta~"Proportion"), color="Source\nPopulation", shape="Source\nPopulation")+
  guides(
    color = guide_legend(
      override.aes = list(
        shape = shape_vals,
        fill = source_cols_alt,
        linetype = 0,
        linewidth = 0,
        alpha = 1)),
    shape = "none",
    fill = "none") +
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
  ) %>% mutate(seg = ifelse(time <= 0.5, "before", "after")) %>% 
  mutate(source = factor(source, levels=sources),
         lake = factor(lake, levels = lakes))


source_cols_g<-source_cols_alt[-1]
shape_vals_g<-shape_vals[-1]

g_plot<-ggplot(g_deltas, aes(x=time, y=delta, color=source))+
  geom_vline(xintercept = 0.5, color="grey", linetype="dashed", linewidth=1)+
  geom_line(aes(color=source), alpha=0.75, linewidth =1.5)+
  geom_point(aes(shape = source, fill=source), size=5.5)+
  scale_shape_manual(values = shape_vals[-1])+
  scale_color_manual(values=source_cols_g)+
  scale_fill_manual(values=source_cols_g)+
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

time_fig<- cowplot::plot_grid(g_plot,loon_plot, ncol = 2, rel_widths = c(1,1.3))
ggsave(time_fig, filename = "03_figures/timespans-fig.tiff",
       device = "tiff",
       width = 1440, height=550, units = "px",
       dpi=300, scale=2.5)
