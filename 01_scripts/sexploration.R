
ancestry_data<-read.csv("02_tidy-data/final-ancestry.csv")

lethal<-rename(read.csv("02_tidy-data/other-data/lethal-data.csv"),fish_id=fishID) %>% 
  select(-year)

data<-ancestry_data %>% 
  left_join(lethal) %>% 
  #filter(year==2) %>% 
  mutate(sex=factor(sex)) %>% 
  filter(!is.na(sex)) %>% 
  select(fish_id,year,lake,sex,ends_with("round"))

loon_data<-filter(data, lake=="LO") %>% 
  mutate(benthic_ancestry=FG_round+TL_round+WB_round+WT_round)

g_data<-filter(data, lake=="GL") %>% 
  mutate(benthic_ancestry=FG_round+TL_round+WB_round+WT_round)

ggplot(loon_data, aes(x=sex, y=benthic_ancestry))+
  geom_violin()+
  geom_jitter(width=0.1, height = 0)+
  facet_wrap(~year, scales="free")

loon_data %>% count(year,sex,benthic_ancestry)
loon_data %>% count(year,sex)
loon_data %>% group_by(year,sex) %>% 
  summarize(mean(benthic_ancestry))

round_cols<-colnames(select(data, ends_with("round")))

sex_data<-loon_data %>% 
  group_by(year,sex) %>% 
  summarize(across(all_of(round_cols), sum, na.rm = TRUE), .groups = "drop") %>% 
  pivot_longer(cols = all_of(round_cols), names_to = "source", values_to = "count") %>%
  group_by(year, sex) %>%
  mutate(prop = count / sum(count)) %>%
  mutate(source = factor(substr(source,1,2), levels=sources))

sex_data<-g_data %>% 
  group_by(sex) %>% 
  summarize(across(all_of(round_cols), sum, na.rm = TRUE), .groups = "drop") %>% 
  pivot_longer(cols = all_of(round_cols), names_to = "source", values_to = "count") %>%
  group_by(sex) %>%
  mutate(prop = count / sum(count)) %>%
  mutate(source = factor(substr(source,1,2), levels=sources))

ggplot(sex_data, aes(x=source, y=prop, fill=sex)) +
  geom_bar(stat = "identity", position="dodge")

  