#### DIVERGENCE TESTS ####

## packages
library(tidyverse)
library(lme4)
library(lmerTest)

## constants
source(here::here("02_tidy-data/constants.R"))

## freq data
freq_data<-read.csv("02_tidy-data/assortment/freq-data.csv") %>% 
  mutate(type=factor(if_else(pure==T,"pure","admixed")),
         o_e_raw=count_obs-count_exp_raw)

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
div_data<-fst %>% 
  left_join(select(adapt_fst,genotype,fst_adapt), by="genotype") %>% 
  mutate(morph=delta_ld) %>% 
  select(genotype,fst,fst_adapt,morph)

## saving
write.csv(div_data, file="02_tidy-data/assortment/divergence-data.csv",
          quote=F, row.names = F)

## export to word doc table
div_data %>%
  mutate(pop1 = substr(genotype, 1, 2),
         pop2 = substr(genotype, 3, 4)) %>%
  select(pop1, pop2, morph, fst, fst_adapt)  %>% 
  arrange(pop1, pop2) %>%
  mutate(across(c(morph, fst, fst_adapt), round, 3)) %>% 
  flextable::flextable() %>% 
  print(preview = "docx")

## merging with freq data
data<-freq_data %>% 
  left_join(div_data, by="genotype") %>% 
  mutate(fst=if_else(is.na(fst),0,fst),
         fst_adapt=if_else(is.na(fst_adapt),0,fst_adapt),
         morph=if_else(is.na(morph),0,morph))

## plotting
data_all<-mutate(data, lake="All")
fig_data<-bind_rows(data,data_all) %>% 
  select(lake, type, count_obs, freq_exp, count_exp_raw, o_e_raw, fst, fst_adapt, morph) %>% 
  mutate(lake=factor(lake, levels=c("All",lakes))) %>% 
  pivot_longer(cols = c(fst, fst_adapt, morph), names_to = "divergence", values_to = "value") %>%
  mutate(divergence=factor(divergence, levels=c("fst","fst_adapt","morph"),
                           labels=c("Genetic (genome-wide)","Genetic (adaptive)","Morphological"))) %>%
  filter(type=="admixed")

fig_data %>% 
  ggplot(aes(x=value, y=o_e_raw))+
  geom_smooth(method="lm", se=F, linetype="dashed",color="grey", linewidth=0.75)+
  geom_point(shape=21, size=2, alpha=0.75, fill="grey45")+
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
      background_y = elem_list_rect(fill = all_lake_cols)
    )
  )

ggsave(filename = "03_figures/assortment/div-fig.tiff",
       device = "tiff",
       width = 600, height=800, units = "px",
       dpi=300, scale=3)


## plotting
data_all<-mutate(data, lake="All")
fig_data<-bind_rows(data,data_all) %>%
  filter(substr(genotype,1,2)=="LG" | substr(genotype,3,4)=="LG") %>% 
  select(lake, type, o_e_raw, fst, fst_adapt, morph) %>% 
  mutate(lake=factor(lake, levels=c("All",lakes))) %>% 
  pivot_longer(cols = c(fst, fst_adapt, morph), names_to = "divergence", values_to = "value") %>%
  mutate(divergence=factor(divergence, levels=c("fst","fst_adapt","morph"),
                           labels=c("Genetic (neutral)","Genetic (adaptive)","Morphological"))) %>%
  filter(type=="admixed")

fig_data %>% 
  ggplot(aes(x=value, y=o_e_raw))+
  geom_smooth(method="lm", se=F, linetype="dashed",color="grey", linewidth=0.75)+
  geom_point(shape=21, size=2, alpha=0.75, fill="grey45")+
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
      background_y = elem_list_rect(fill = all_lake_cols)
    )
  )


## linear models 
lm_results<-fig_data %>% 
  group_by(lake, divergence) %>% 
  group_modify(~ {
    model <- lm(o_e_raw ~ scale(value), data = .x)
    tidy_model <- broom::tidy(model)
    glance_model <- broom::glance(model)
    results <- tidy_model %>%
      filter(term == "scale(value)") %>%
      select(estimate, std.error, statistic, p.value) %>%
      bind_cols(glance_model %>% select(r.squared, adj.r.squared, sigma, df.residual))
    results
  }) %>% 
  ungroup() %>%
  group_by(divergence) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

## export
lm_results %>% 
  filter(lake != "All") %>% 
  group_by(divergence) %>% 
  mutate(p_adj = p.adjust(p.value, method = "BH"),
         lake=factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes))) %>%
  ungroup() %>% 
  arrange(lake, divergence) %>%
  select(lake, divergence, estimate, std.error, statistic, p.value, p_adj, adj.r.squared) %>%
  mutate(across(c(estimate, std.error, statistic, p.value, p_adj, adj.r.squared), round, digits = 3)) %>%
  flextable::flextable() %>% 
  print(preview = "docx")
  
## all models
sig_lakes<-c("GL","FL","RL","LP")

lmer(o_e_raw~scale(fst)+(1|lake), data=filter(data, type=="admixed")) %>% r2()
lmer(o_e_raw~scale(fst)+(1|lake), data=filter(data, type=="admixed", lake %in% sig_lakes)) %>% r2()

lmer(o_e_raw~scale(fst_adapt)+(1|lake), data=filter(data, type=="admixed")) %>% r2()
lmer(o_e_raw~scale(fst_adapt)+(1|lake), data=filter(data, type=="admixed", lake %in% sig_lakes)) %>% r2()

lmer(o_e_raw~scale(morph)+(1|lake), data=filter(data, type=="admixed")) %>% r2()
lmer(o_e_raw~scale(morph)+(1|lake), data=filter(data, type=="admixed", lake %in% sig_lakes)) %>% r2()

## model
fst_mod <- lmer(o_e~scale(fst)+(1|lake), data=filter(data, type=="admixed", lake %in% c("RL","GL","FL","LP")))
summary(fst_mod)

morph_mod <- lmer(o_e~scale(morph)+(1|lake), data=filter(data, type=="admixed", lake %in% c("RL","GL","FL","LP")))
summary(morph_mod)

lm(o_e~morph, data=filter(data, type=="admixed", lake=="FL")) %>% summary()

ggplot(data, aes(x=fst_adapt, y=morph))+
  geom_point()+
  geom_smooth(method="lm", se=F, linetype="dashed",color="grey", linewidth=0.75)+
  labs(x="Adaptive Fst", y="Morphological Divergence")+
  theme_bw()

#### NEW MODELS ####

## all lakes data
all_data<-data %>% 
  filter(type=="admixed") %>% 
  pivot_longer(cols = c(fst, fst_adapt, morph), names_to = "divergence", values_to = "value") %>%
  mutate(divergence=factor(divergence, levels=c("fst","fst_adapt","morph"),
                           labels=c("Genetic (genome-wide)","Genetic (adaptive)","Morphological")))

## glmm models
divergence_types <- unique(all_data$divergence)  # or whatever column labels the 3 divergence measures
all_models <- map(divergence_types, function(div) {
  glmer(
    count_obs ~ value + offset(log(count_exp_raw)) + (1 | lake),
    family = poisson,
    data = filter(all_data, divergence == div)
  )
})

names(all_models) <- divergence_types

## glm models
glm_results <- fig_data %>%
  filter(lake!="All") %>%
  group_by(lake, divergence) %>%
  group_modify(~ {
    
    model <- glm(
      count_obs ~ scale(value) +
        offset(log(count_exp_raw)),
      family = poisson,
      data = .x
    )
    
    tidy_model  <- broom::tidy(model)
    glance_model <- broom::glance(model)
    
    tidy_model %>%
      filter(term == "scale(value)") %>%
      select(estimate, std.error, statistic, p.value) %>%
      bind_cols(
        glance_model %>%
          select(deviance, df.residual, AIC)
      )
  }) %>%
  ungroup() %>%
  group_by(divergence) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()


## new plot
lake_div_models <- fig_data %>%
  filter(lake != "All") %>%
  group_by(lake, divergence) %>%
  group_map(~ glm(
    count_obs ~ scale(value) + offset(log(count_exp_raw)),
    family = poisson,
    data = .x
  ))

# Save names for referencing
lake_div_names <- fig_data %>%
  filter(lake != "All") %>%
  distinct(lake, divergence) %>%
  mutate(model_index = row_number())

names(lake_div_models) <- paste(lake_div_names$lake, lake_div_names$divergence, sep = "_")

## prediction grid
pred_data <- fig_data %>%
  filter(lake != "All") %>%
  group_by(lake, divergence) %>%
  summarise(
    value = seq(min(value), max(value), length.out = 100),
    count_exp_raw = mean(count_exp_raw),
    .groups = "drop"
  )

pred_data <- pred_data %>%
  rowwise() %>%
  mutate(
    pred = predict(
      lake_div_models[[paste(lake, divergence, sep = "_")]],
      newdata = data.frame(value = value, count_exp_raw = count_exp_raw),
      type = "response"
    )
  ) %>%
  ungroup() %>%
  mutate(ratio_pred = pred / count_exp_raw)

plot_data <- fig_data %>%
  filter(lake != "All") %>%
  mutate(ratio = count_obs / count_exp_raw)

ggplot(plot_data, aes(x = value, y = ratio)) +
  geom_point(alpha = 0.6) +
  geom_line(
    data = pred_data,
    aes(x = value, y = ratio_pred, colour = divergence),
    size = 1
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid2(
    rows = vars(lake),
    cols = vars(divergence),
    scales = "free_x",
    strip = strip_themed(
      background_y = elem_list_rect(fill = all_lake_cols[-1])
    )
  )

## plotting all
plot_data_all <- data %>%
  filter(type == "admixed") %>%
  mutate(ratio = count_obs / count_exp_raw) 

pred_data_all<- map_df(names(all_models), function(div) {
  
  df_sub <- filter(all_data, divergence == div)
  
  tibble(
    value = seq(min(df_sub$value), max(df_sub$value), length.out = 100),
    divergence = div,
    count_exp_raw = mean(df_sub$count_exp_raw)
  ) %>%
    rowwise() %>%
    mutate(
      pred = predict(
        all_models[[div]],
        newdata = data.frame(
          value = value,
          count_exp_raw = count_exp_raw
        ),
        re.form = NA,  # removes random effect → population-average line
        type = "response"
      ),
      ratio_pred = pred / count_exp_raw
    ) %>%
    ungroup()
}) %>% 
  mutate(divergence = factor(divergence, levels = c("Genetic (genome-wide)", "Genetic (adaptive)", "Morphological")))

div_1<-all_data %>% 
  mutate(ratio=count_obs/count_exp_raw,
         lake="All") %>% 
  ggplot(aes(x=value, y=ratio))+
  geom_line(
    data = mutate(pred_data_all, lake="All"),
    aes(x = value, y = ratio_pred),
    size = 0.75, color = "grey45", linetype = "dashed"
  ) +
  geom_point(shape=21, size=2, alpha=0.75, fill="grey45")+
  labs(y="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, hjust=0.75))+
  facet_grid2(
    rows = vars(lake),
    cols = vars(divergence),
    scales = "free_x",
    strip = strip_themed(
      background_y = elem_list_rect(fill = "grey95"),
      background_x = elem_list_rect(fill = "grey95")))
  
div_2<-plot_data %>% 
  filter(lake %in% limnetic_lakes) %>% 
  ggplot(aes(x=value, y=ratio))+
  geom_line(
    data = filter(pred_data, lake %in% limnetic_lakes),
    aes(x = value, y = ratio_pred),
    size = 0.75, color = "grey45", linetype = "dashed"
  ) +
  geom_point(shape=21, size=2, alpha=0.75, fill="grey45")+
  labs(y="Obs/Exp Counts")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12, hjust=0.70))+
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

div_3<-plot_data %>% 
  filter(lake %in% benthic_lakes) %>% 
  ggplot(aes(x=value, y=ratio))+
  geom_line(
    data = filter(pred_data, lake %in% benthic_lakes),
    aes(x = value, y = ratio_pred),
    size = 0.75, color = "grey45", linetype = "dashed"
  ) +
  geom_point(shape=21, size=2, alpha=0.75, fill="grey45")+
  theme(axis.title = element_blank())+
  facet_grid2(
    rows = vars(lake),
    cols = vars(divergence),
    scales = "free_x",
    strip = strip_themed(
      background_y = elem_list_rect(fill = benthic_col),
      background_x = elem_list_rect(fill = "grey95")))

## both lakes
div_4<-plot_data %>% 
  filter(lake %in% both_lakes) %>% 
  ggplot(aes(x=value, y=ratio))+
  geom_line(
    data = filter(pred_data, lake %in% both_lakes),
    aes(x = value, y = ratio_pred),
    size = 0.75, color = "grey45", linetype = "dashed"
  ) +
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
  filename = "03_figures/assortment/div-fig-new-new.tiff",
  width = 1, height=0.5)


# ── Data prep ────────────────────────────────────────────────────────────────

all_data <- data %>%
  filter(type == "admixed") %>%
  pivot_longer(cols = c(fst, fst_adapt, morph),
               names_to  = "divergence",
               values_to = "value") %>%
  mutate(divergence = factor(divergence,
                             levels = c("fst", "fst_adapt", "morph"),
                             labels = c("Genetic (genome-wide)",
                                        "Genetic (adaptive)",
                                        "Morphological")))

# ── GLMMs (one per divergence type, lake as random effect) ───────────────────

all_models <- map(levels(all_data$divergence), function(div) {
  glmer(count_obs ~ value + offset(log(count_exp_raw)) + (1 | lake),
        family = poisson,
        data   = filter(all_data, divergence == div))
}) %>% set_names(levels(all_data$divergence))

## GLM subsets
sub_data<-all_data %>% 
  filter(lake %in% c("GL","FL","RL","LP"))

sub_models <- map(levels(sub_data$divergence), function(div) {
  glmer(count_obs ~ value + offset(log(count_exp_raw)) + (1 | lake),
        family = poisson,
        data   = filter(sub_data, divergence == div))
}) %>% set_names(levels(sub_data$divergence))

## GLMM results
tidy_models <- function(model_list) {
  imap_dfr(model_list, ~ tidy(.x, conf.int = TRUE) %>% 
             mutate(model = .y))
}

# Combine all results
summary_all <- tidy_models(all_models) %>% mutate(group = "All")
summary_sub <- tidy_models(sub_models) %>% mutate(group = "Subset")

# Final combined dataframe
final_df <- bind_rows(summary_all, summary_sub) %>%
  filter(effect == "fixed", term=="value") %>%
  select(group, model, estimate, std.error, statistic, p.value)
  
## exporting
final_df %>% 
  mutate(across(c(estimate, std.error, statistic), round, 3),
         p.value=if_else(p.value<0.001, "<0.001", as.character(round(p.value, 3)))) %>%
  flextable::flextable() %>% 
  print(preview = "docx")

## per lake GLMs

glm_results <- all_data %>%
  group_by(lake, divergence) %>%
  group_modify(~ {
    model <- glm(count_obs ~ scale(value) + offset(log(count_exp_raw)),
                 family = poisson, data = .x)
    broom::tidy(model) %>%
      filter(term == "scale(value)") %>%
      select(estimate, std.error, statistic, p.value) %>%
      bind_cols(broom::glance(model) %>% select(deviance, df.residual, AIC))
  }) %>%
  ungroup() %>%
  group_by(divergence) %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

## export
glm_results %>% 
  mutate(lake=factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes))) %>%
  mutate(across(c(estimate, std.error, statistic), round, 3),
         p.value=if_else(p.value<0.001, "<0.001", as.character(round(p.value, 3))),
         p_adj=if_else(p_adj<0.001, "<0.001", as.character(round(p_adj, 3)))) %>%
  arrange(lake, divergence) %>%
  select(lake, divergence, estimate, std.error, statistic, p.value, p_adj) %>%
  flextable::flextable() %>% 
  print(preview = "docx")

#### PREDICTION GRIDS ####

# Per-lake predictions (re-use the GLMs from glm_results via group_map)
lake_div_models <- all_data %>%
  group_by(lake, divergence) %>%
  group_map(~ glm(count_obs ~ scale(value) + offset(log(count_exp_raw)),
                  family = poisson, data = .x)) %>%
  set_names(all_data %>%
              distinct(lake, divergence) %>%
              mutate(nm = paste(lake, divergence, sep = "_")) %>%
              pull(nm))

pred_data <- all_data %>%
  group_by(lake, divergence) %>%
  summarise(value         = seq(min(value), max(value), length.out = 100),
            count_exp_raw = mean(count_exp_raw),
            .groups = "drop") %>%
  rowwise() %>%
  mutate(ratio_pred = predict(
    lake_div_models[[paste(lake, divergence, sep = "_")]],
    newdata = data.frame(value = value, count_exp_raw = count_exp_raw),
    type = "response") / count_exp_raw) %>%
  ungroup()

# All-lakes predictions (population-average, no random effect)
pred_data_all <- map_df(levels(all_data$divergence), function(div) {
  df_sub <- filter(all_data, divergence == div)
  tibble(value         = seq(min(df_sub$value), max(df_sub$value), length.out = 100),
         divergence    = div,
         count_exp_raw = mean(df_sub$count_exp_raw)) %>%
    rowwise() %>%
    mutate(ratio_pred = predict(
      all_models[[div]],
      newdata  = data.frame(value = value, count_exp_raw = count_exp_raw),
      re.form  = NA,
      type     = "response") / count_exp_raw) %>%
    ungroup()
}) %>%
  mutate(divergence = factor(divergence, levels = levels(all_data$divergence)))

plot_data <- all_data %>%
  mutate(ratio = count_obs / count_exp_raw)

#### PLOTTING ####

point_layer  <- geom_point(shape = 21, size = 2, alpha = 0.75, fill = "grey45")
line_layer   <- function(df) geom_line(data = df, aes(x = value, y = ratio_pred),
                                       size = 0.75, color = "grey45", linetype = "dashed")
no_x_title   <- theme(axis.title.x = element_blank())

## all lakes
div_1 <- plot_data %>%
  mutate(lake = "All") %>%
  ggplot(aes(x = value, y = ratio)) +
  line_layer(mutate(pred_data_all, lake = "All")) +
  point_layer +
  labs(y = "") +
  no_x_title +
  theme(axis.title.y = element_text(size = 12, hjust = 0.75)) +
  facet_grid2(rows = vars(lake), cols = vars(divergence), scales = "free_x",
              strip = strip_themed(background_y = elem_list_rect(fill = "grey95"),
                                   background_x = elem_list_rect(fill = "grey95")))

## limnetic lakes
div_2 <- plot_data %>%
  filter(lake %in% limnetic_lakes) %>%
  ggplot(aes(x = value, y = ratio)) +
  line_layer(filter(pred_data, lake %in% limnetic_lakes)) +
  point_layer +
  labs(y = "Obs/Exp Counts") +
  no_x_title +
  theme(axis.title.y = element_text(size = 12, hjust = 0.70)) +
  facet_grid2(rows = vars(lake), cols = vars(divergence), scales = "free_x",
              strip = strip_themed(background_y = elem_list_rect(fill = limnetic_col),
                                   background_x = element_blank(),
                                   text_x       = element_blank()))

## benthic lakes
div_3 <- plot_data %>%
  filter(lake %in% benthic_lakes) %>%
  ggplot(aes(x = value, y = ratio)) +
  line_layer(filter(pred_data, lake %in% benthic_lakes)) +
  point_layer +
  theme(axis.title = element_blank()) +
  facet_grid2(rows = vars(lake), cols = vars(divergence), scales = "free_x",
              strip = strip_themed(background_y = elem_list_rect(fill = benthic_col),
                                   background_x = elem_list_rect(fill = "grey95")))

## both lakes 
div_4 <- plot_data %>%
  filter(lake %in% both_lakes) %>%
  ggplot(aes(x = value, y = ratio)) +
  line_layer(filter(pred_data, lake %in% both_lakes)) +
  point_layer +
  theme(axis.title = element_blank()) +
  facet_grid2(rows = vars(lake), cols = vars(divergence), scales = "free_x",
              strip = strip_themed(background_y = elem_list_rect(fill = both_col),
                                   background_x = element_blank(),
                                   text_x       = element_blank()))

## combining
div_plot <- plot_grid(
  plot_grid(div_1, div_2, ncol = 1, rel_heights = c(1.05, 3)),
  plot_grid(div_3, div_4, ncol = 1, rel_heights = c(4, 2.5)),
  ncol = 2, rel_widths = c(1, 1)
)

## saving
save_fig(plot     = div_plot,
         filename = "03_figures/assortment/div-fig-new-new.tiff",
         width    = 1,
         height   = 0.5)

