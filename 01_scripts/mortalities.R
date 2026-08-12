#### MORTALITY ANALYSIS ####

## Description: exploring patterns in pre-transplant mortalities
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## constants
source(here::here("02_tidy-data/constants.R"))

## data
ancestry_data<-read.csv("./02_tidy-data/final-ancestry.csv")
data19<-filter(ancestry_data,year=="19")

## mortalities
score_cols<-c("FG_score","TL_score","WB_score","WT_score","LG_score","SL_score","SR_score","WK_score")
mortalities<-data19 %>% rowwise() %>% mutate(source = str_remove(score_cols[which.max(c_across(all_of(score_cols)))], "_score")) %>%
  ungroup() %>% select(fish_id,lake,source) %>% mutate(fish_id=substr(fish_id,1,8)) %>% unique() %>% 
  group_by(lake,source) %>% summarise(fish=n()) %>% 
  mutate(source = factor(source, levels=sources),
         lake = factor(lake, levels = lakes))
mortalities$fish[is.na(mortalities$fish)]<-0

## stocking data
stocking<-read.csv("./02_tidy-data/other-data/stocking.csv") %>% pivot_longer(cols = FG:WK) %>% 
  select(lake=recipient_lake,source=name,stocked=value) %>% 
  filter(stocked>0)

## mortality rates
mort_df<-mortalities %>% rename(dead=fish) %>% full_join(stocking) %>% filter(!is.na(lake), lake!="GL")
mort_df$dead[is.na(mort_df$dead)]<-0
mort_df<-mort_df %>% mutate(mort_rate=dead/stocked,
                            actual_stocked=stocked-dead)

## stocking
stock_df<-mort_df %>%
  select(lake,source,actual_stocked) %>% 
  pivot_wider(names_from = source, values_from = actual_stocked) %>% 
  mutate(total = rowSums(across(where(is.numeric)), na.rm = TRUE))
stock_df[is.na(stock_df)]<-0

## adding g
g_stocks <- data.frame(
  lake = "GL",
  TL = 500,
  WB = 500,
  LG = 0,
  SL = 500,
  SR = 500,
  WK = 495,
  FG = 500,
  WT = 500,
  total = 3495)

stock_df <- stock_df %>% bind_rows(g_stocks)

## saving
write.csv(stock_df, file = "02_tidy-data/other-data/actual-stocks.csv",
          row.names = F, quote = F)

## means
mort_means<-mort_df %>% group_by(source) %>% summarize(dead=sum(dead),stocked=sum(stocked)) %>% 
  mutate(overall_rate=dead/stocked)
mort_df<-left_join(mort_df,select(mort_means,source,overall_rate)) %>% 
  mutate(source = factor(source, levels=sources),
         lake = factor(lake, levels = lakes))

## saving
write.csv(mort_df, file = "02_tidy-data/other-data/mortality-rates.csv",
          quote = F, row.names = F)

## plotting 
stocked_labels <- mort_df %>%
  group_by(lake) %>%
  summarise(stocked = sum(stocked))
dead_totals <- mort_df %>%
  group_by(lake) %>%
  summarise(total_dead = sum(dead)) %>%
  left_join(stocked_labels, by = "lake")

ggplot(mort_df, aes(x=lake,y=dead,fill=source))+
  geom_bar(stat = "identity")+
  geom_text(data = dead_totals,
            aes(x = lake, y = total_dead + 5, label = paste0(stocked," total")),  # +5 for spacing
            inherit.aes = FALSE,
            size = 6) +
  theme_bw()+
  scale_fill_manual(values = source_cols_alt)+
  labs(x="Lake",y="Dead",fill="Source")+
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size = 22),
        legend.position = "top")

## plotting
mean_rate<-sum(mort_df$dead)/sum(mort_df$stocked)
ggplot(mort_df)+
  geom_hline(yintercept = mean_rate,color="grey",linetype="dashed", linewidth=1)+
  geom_point(aes(x=source,y=mort_rate,size=dead,fill=source),shape=21, alpha=0.75)+
  scale_fill_manual(values = source_cols_alt)+
  scale_size_continuous(range = c(1.5,8))+
  geom_point(aes(x=source,y=overall_rate),color="black",size=5,shape=3,stroke=1)+
  labs(x="",y="Mortality Rate")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=14))

ggsave(filename = "03_figures/mortality-fig.tiff",
       device = "tiff",
       width = 1200, height=600, units = "px",
       dpi=300, scale=2.5)
