#### TRAIT MODELS ####

## Description: correlations between traits and general population success
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## packages
library(tidyverse)
library(purrr)

## for factors
sources<-c("LG","SR","WK","SL","WT","WB","FG","TL")
lakes<-c("CL","FL","HL","RL","CC","LL","LP","LO")
source_cols_alt<-c("#08519c","#3182bd","#6baed6","#bdd7e7",
                   "#006d2c","#31a354","#74c476","#bae4b3")

#### data import ####

## success data
success<-read.csv("02_tidy-data/source-effect-means-new.csv")

## genetic diversity
diversity<-read.csv("02_tidy-data/other-data/ak_diversity.csv") %>%
  select(ends_with("pi"), ends_with("watterson")) %>% 
  pivot_longer(cols=c(ends_with("pi"),ends_with("watterson")),names_to = "source") %>% 
  separate(source, into = c("source", "metric"), sep = "\\.", remove = FALSE)

ggplot(diversity, aes(x=source, y=value, fill=metric))+
  geom_bar(stat="identity", position="dodge")

diversity %>% pivot_wider(names_from = metric) %>% 
  ggplot(aes(x=theta_watterson, y=theta_pi,label=source))+
  geom_text()

## mortality rates
morts<-read.csv("02_tidy-data/other-data/mortality-rates.csv") %>% 
  select(source,mortality_rate=overall_rate) %>% 
  unique()

## max body size
lethal<-read.csv("02_tidy-data/other-data/lethal-data.csv") %>% 
  mutate(source=factor(substr(fishID,1,2)),
         year=factor(year)) %>% 
  filter(source %in% sources)

body_size<-lethal %>% group_by(source,year) %>% 
  summarize(max_body_size=quantile(std_length,0.95,na.rm=T),
            median_body_size=median(std_length, na.rm=T))

ggplot(body_size, aes(x=source, y=max_body_size, fill=year))+
  geom_bar(stat="identity", position="dodge")

body_size %>% 
  filter(year=="2019") %>% 
  ggplot(aes(x=median_body_size, y=max_body_size))+
  geom_point()

## infection rate
parasites<-lethal %>%
  mutate(parasite=if_else(parasite_no>0,1,0)) %>% 
  group_by(source,year) %>% 
  summarize(parasite_rate=sum(parasite)/n()) %>% 
  mutate(parasite_rate_sqrt=sqrt(parasite_rate))

ggplot(parasites, aes(x=source, y=parasite_rate, fill=year))+
  geom_bar(stat="identity", position="dodge")

## gravidity
gravidity<-lethal %>%
  group_by(source,year) %>% 
  summarize(gravid_prop=sum(gravid, na.rm=T)/n())

ggplot(gravidity, aes(x=source, y=gravid_prop, fill=year))+
  geom_bar(stat="identity", position="dodge")

## combining
data<-success %>%
  left_join(pivot_wider(diversity,names_from = metric)) %>% 
  left_join(morts) %>% 
  left_join(filter(body_size, year=="2019")) %>% 
  left_join(filter(parasites, year=="2019")) %>% 
  left_join(filter(gravidity, year=="2019")) %>% 
  select(-year, -theta_watterson, -parasite_rate, -gravid_prop, -max_body_size)

#### graphs #####
graph_data<-data %>% pivot_longer(cols = c(mean_effect_1, mean_effect_net),
                      names_to = "timespan",
                      names_prefix = "mean_effect_",
                      values_to = "response_value") %>% 
  pivot_longer(cols = c(theta_pi,mortality_rate,median_body_size,parasite_rate_sqrt),
               names_to = "predictor_var",
               values_to = "predictor_value")

predictor_labels <- c(
  "median_body_size" = "Median Body Size (mm)",
  "mortality_rate" = "Mortality Rate",
  "parasite_rate_sqrt" = "Parasite Prevelance (sqrt)",
  "theta_pi" = "Genetic Diversity (pi)")

timespan_labels <- c(
  "1" = "Year 1 Change",
  "net" = "Net Change")

shape_vals<-rep(c(18,20),4)

graph_data %>%
  mutate(source=factor(source, levels=sources)) %>%
  ggplot(aes(x=predictor_value,y=response_value, color=source))+
  geom_smooth(aes(group = 1), method = "lm", color = "grey", linetype="dashed", se = FALSE) +
  geom_point(aes(shape=source), size=7)+
  scale_shape_manual(values = shape_vals)+
  scale_color_manual(values = source_cols_alt)+
  facet_grid(rows = vars(timespan), cols = vars(predictor_var), scales = "free",
             labeller = labeller(
               timespan = timespan_labels,
               predictor_var = predictor_labels
             )) +
  labs(y="Mean Source Effect Size", color="Source\nPopulation", shape="Source\nPopulation")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16, face="bold"),
        strip.background = element_rect(fill = "grey95"),
        legend.text = element_text(size=16),
        legend.title = element_text(size = 18, hjust=0.5, face="bold"))

ggsave(filename = "03_figures/trait-success-new.tiff",
       device = "tiff",
       width = 1440, height=800, units = "px",
       dpi=300, scale=3)

#### models ####
combs<-expand.grid(
  timespan = unique(graph_data$timespan),
  predictor_var = unique(graph_data$predictor_var))

models<-list()
for (i in 1:nrow(combs)){
  time=combs$timespan[i]
  predictor<-combs$predictor_var[[i]]
  response <- paste0("mean_effect_", as.character(time))
  
  formula_str <- paste0("scale(", response, ") ~ scale(", predictor, ")")
  
  lm<-lm(formula = formula_str, data=data) %>% summary()
  models[[i]]<-lm
  names(models)[i]<-formula_str
}

## summary table
summary_table <- map_df(models, ~{
  tibble(
    effect = .$coefficients[2,1],
    r2 = .$r.squared,
    adj_r2 = .$adj.r.squared,
    p_value = pf(.$fstatistic[1], .$fstatistic[2], .$fstatistic[3], lower.tail = FALSE)
  )
}, .id = "formula") %>% 
  mutate(time=str_extract(formula, "mean_effect_[^ ]+")) %>%
  group_by(time) %>%
  mutate(adj_p_value=p.adjust(p_value, method = "BH")) %>%
  ungroup()
  
## export
summary_table %>% 
  select(time, formula, effect, r2, adj_r2, p_value, adj_p_value) %>% 
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  arrange(time) %>%
  flextable::flextable() %>% 
  flextable::save_as_docx(path="03_figures/tables/trait-lm-results.docx")


lm(scale(mean_effect_1) ~ scale(mortality_rate) + scale(theta_pi), data=data) %>% summary()
