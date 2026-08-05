#### SELECTION NEW ####

## constants
source(here::here("02_tidy-data/constants.R"))

## freq data
freq_data<-read.csv("02_tidy-data/assortment/freq-data.csv") %>%
  select(lake, genotype, count_exp=count_exp_raw, count_obs)

## rawer data
data_f1<-read.csv("02_tidy-data/assortment/final-f1-data.csv") %>% 
  select(lake, year, ends_with("round"))

colnames(data_f1)<-str_remove(colnames(data_f1),"_round")

## hwe data
hwe_data<-data_f1 %>%
  rowwise() %>%
  mutate(
    alleles = list(names(pick(-c(lake,year)))[pick(-c(lake,year)) > 0]),
    geno = if(length(alleles) == 2) {
      paste0(sort(alleles), collapse = "")
    } else if(length(alleles) == 1) {
      paste0(alleles, alleles)
    } else {
      NA_character_
    }
  ) %>%
  ungroup() %>% 
  mutate(pop=paste(lake,year,sep="_")) %>%
  filter(!is.na(geno)) %>% 
  mutate(geno=paste0(substr(geno,1,2),"/",substr(geno,3,4)))

## genind data
genind_df<-hwe_data %>%
  select(geno) %>% 
  df2genind(
    sep="/",
    ncode=2,
    pop=hwe_data$pop,
    ploidy=2)

## f stats
hf_data <- genind2hierfstat(genind_df)
basic_stats<-basic.stats(hf_data)
fis<- basic_stats$Fis[1,]

## tests
hw.test(genind_df, B = 10000)
pop_split<-seppop(genind_df)
hwe_by_lake<-lapply(pop_split, function(x) hw.test(x, B = 10000))

## summary
pop_summary<-data.frame(
  lake_year = names(basic_stats$Fis[1,]),
  n = sapply(pop_split, nInd),
  He = basic_stats$Hs[1,],
  Ho = basic_stats$Ho[1,],
  Fis = basic_stats$Fis[1,],
  p_value = sapply(hwe_by_lake, function(x) x[, "Pr.exact"]),
  chi_sq   = sapply(hwe_by_lake, function(x) x[, "chi^2"]),
  df      = sapply(hwe_by_lake, function(x) x[, "df"])) %>% 
  mutate(p_adj=p.adjust(p_value, method = "BH"),
         group="source") %>% 
  mutate(lake=substr(lake_year,1,2),
         year=substr(lake_year,4,4))

#### props ####
props<-hwe_data %>% 
  group_by(pop) %>% 
  summarize(
    total = n(),
    prop_pure = sum(substr(geno, 1,2) == substr(geno, 4,5)) / n(),
    prop_admixed = sum(substr(geno, 1,2) != substr(geno, 4,5)) / n()
  ) %>% 
  ungroup() %>% 
  mutate(lake=substr(pop,1,2),
         year=substr(pop,4,4)) %>% 
  select(-pop) %>% 
  mutate(year=as.numeric(year))

prop_fig<-props %>% mutate(lake=factor(lake, levels=c(limnetic_lakes,benthic_lakes,both_lakes))) %>%
ggplot(aes(x=year, y=prop_admixed))+
  geom_point(aes(size=total))+
  geom_line(linewidth = 0.5,color = "black",linetype = "dashed")+
  scale_x_continuous(breaks=c(1,2), labels=c("1","2"), limits = c(0.5,2.5))+
  labs(x="Year", y="Proportion Admixed", size="Sample Size")+
  theme(legend.position = "top",
        legend.title = element_text(size=12, hjust=0.5))+
  facet_grid2(
    cols = vars(lake),
    strip = strip_themed(
      background_x = elem_list_rect(fill = all_lake_cols[-1])))

save_fig(
  plot = prop_fig,
  filename = "03_figures/assortment/selection-fig-2.tiff",
  width = 1, height=0.33)

fisher_results <- props %>%
  filter(lake!="GL") %>% 
  mutate(
    n_admixed = round(total * prop_admixed),
    n_pure = total - n_admixed
  ) %>%
  group_by(lake) %>%
  summarize(
    # We create the matrix inside the function
    p_val = fisher.test(matrix(c(n_admixed, n_pure), nrow = 2))$p.value,
    # Calculate the change in proportion for context
    prop_diff = prop_admixed[year == "2"] - prop_admixed[year == "1"]
  ) %>%
  mutate(p_adj = p.adjust(p_val, method = "BH"))



