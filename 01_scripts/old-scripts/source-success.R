## SOURCE SUCCESS ## 

## packages
library(tidyverse)

## loading data
data<-read_csv("./02_tidy-data/final-ancestry-allsnps.csv")
data$lake<-factor(data$lake)

data<-data %>% mutate(year=factor(substr(fish,3,4)))

## grouping by pop and computing proportional success
popData<- as.data.frame(data %>% group_by(lake,year) %>% 
                          summarise(FG_total=sum(FG, na.rm = T),
                                    TL_total=sum(TL, na.rm = T),
                                    WB_total=sum(WB, na.rm = T),
                                    WT_total=sum(WT, na.rm = T),
                                    LG_total=sum(LG, na.rm = T),
                                    SL_total=sum(SL, na.rm = T),
                                    SR_total=sum(SR, na.rm = T),
                                    WK_total=sum(WK, na.rm = T),)) %>%
  mutate(FG_prop=(FG_total/(FG_total+TL_total+WB_total+WT_total+LG_total+SL_total+SR_total+WK_total)),
         TL_prop=(TL_total/(FG_total+TL_total+WB_total+WT_total+LG_total+SL_total+SR_total+WK_total)),
         WB_prop=(WB_total/(FG_total+TL_total+WB_total+WT_total+LG_total+SL_total+SR_total+WK_total)),
         WT_prop=(WT_total/(FG_total+TL_total+WB_total+WT_total+LG_total+SL_total+SR_total+WK_total)),
         LG_prop=(LG_total/(FG_total+TL_total+WB_total+WT_total+LG_total+SL_total+SR_total+WK_total)),
         SL_prop=(SL_total/(FG_total+TL_total+WB_total+WT_total+LG_total+SL_total+SR_total+WK_total)),
         SR_prop=(SR_total/(FG_total+TL_total+WB_total+WT_total+LG_total+SL_total+SR_total+WK_total)),
         WK_prop=(WK_total/(FG_total+TL_total+WB_total+WT_total+LG_total+SL_total+SR_total+WK_total)),) %>% 
  gather(source,prop,FG_prop:WK_prop) %>% 
  mutate(source=as.factor(substr(source,1,2)))

## ordering levels
popData$source<-factor(popData$source, levels = c("LG","SL","SR","WK","FG","TL","WB","WT"))
popData$lake<-factor(popData$lake, levels = c("CL","FL","HL","RL","CC","LL","LP","LO"))

## plotting success
ggplot(popData, aes(x=lake, y=prop, fill=source))+
  geom_bar(position = "stack", stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = c("#eff3ff","#bdd7e7","#6baed6","#2171b5",
                               "#feedde","#fdbe85","#fd8d3c","#d94701"))+
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

## alternative colour scheme
ggplot(popData, aes(x=lake, y=prop, fill=source))+
  geom_bar(position = "stack", stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = c("#c1cad0","#8697a3","#4d6778","#0f3b50",
                               "#d1cebe","#a39f81","#777248","#4b4910"))+
  xlab("")+
  ylab("Proportional Contribution")+
  guides(fill=guide_legend(title="Source Population:"))+
  theme(legend.position = "top", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        axis.text.y = element_text(size=18))


## need to make this a change in proportion using the stocking data

## individual level plots
LO <- data %>% filter(lake=="LO") %>% 
  select(fish, FG_round, TL_round, WB_round, WT_round, LG_round, SL_round, SR_round, WK_round) %>% 
  mutate(type=factor(case_when(FG_round == 1 ~ "pure",
                               TL_round == 1 ~ "pure",
                               WB_round == 1 ~ "pure",
                               WT_round == 1 ~ "pure",
                               LG_round == 1 ~ "pure",
                               SL_round == 1 ~ "pure",
                               SR_round == 1 ~ "pure",
                               WK_round == 1 ~ "pure",
                               FG_round == 0.25 ~ "F2",
                               TL_round == 0.25 ~ "F2",
                               WB_round == 0.25 ~ "F2",
                               WT_round == 0.25 ~ "F2",
                               LG_round == 0.25 ~ "F2",
                               SL_round == 0.25 ~ "F2",
                               SR_round == 0.25 ~ "F2",
                               WK_round == 0.25 ~ "F2",
                               FG_round == 0.5 & TL_round != 0.25 & WB_round != 0.25 & WT_round != 0.25 & LG_round != 0.25 & SL_round != 0.25 & SR_round != 0.25 & WK_round != 0.25  ~ "F1",
                               TL_round == 0.5 & FG_round != 0.25 & WB_round != 0.25 & WT_round != 0.25 & LG_round != 0.25 & SL_round != 0.25 & SR_round != 0.25 & WK_round != 0.25  ~ "F1",
                               WB_round == 0.5 & TL_round != 0.25 & FG_round != 0.25 & WT_round != 0.25 & LG_round != 0.25 & SL_round != 0.25 & SR_round != 0.25 & WK_round != 0.25  ~ "F1",
                               WT_round == 0.5 & TL_round != 0.25 & WB_round != 0.25 & FG_round != 0.25 & LG_round != 0.25 & SL_round != 0.25 & SR_round != 0.25 & WK_round != 0.25  ~ "F1",
                               LG_round == 0.5 & TL_round != 0.25 & WB_round != 0.25 & WT_round != 0.25 & FG_round != 0.25 & SL_round != 0.25 & SR_round != 0.25 & WK_round != 0.25  ~ "F1",
                               SL_round == 0.5 & TL_round != 0.25 & WB_round != 0.25 & WT_round != 0.25 & LG_round != 0.25 & FG_round != 0.25 & SR_round != 0.25 & WK_round != 0.25  ~ "F1",
                               SR_round == 0.5 & TL_round != 0.25 & WB_round != 0.25 & WT_round != 0.25 & LG_round != 0.25 & SL_round != 0.25 & FG_round != 0.25 & WK_round != 0.25  ~ "F1",
                               WK_round == 0.5 & TL_round != 0.25 & WB_round != 0.25 & WT_round != 0.25 & LG_round != 0.25 & SL_round != 0.25 & SR_round != 0.25 & FG_round != 0.25  ~ "F1")))

indLO<- LO %>% gather(source, prop, FG_round:WK_round) %>%
  mutate(source=substr(source,1,2), fish=as.factor(fish)) 
indLO$source<-factor(indLO$source, levels = c("LG","SL","SR","WK","FG","TL","WB","WT"))

orderLO<- LO %>% rowwise() %>%
  mutate(zeros = sum(c_across(FG_round:WK_round) == 0)) %>%
  as.data.frame()
orderLO$type<-factor(orderLO$type, levels = c("pure","F1","F2"))
orderLO <- arrange(orderLO, type, desc(zeros),desc(FG_round),desc(TL_round),desc(WB_round),desc(WT_round),desc(LG_round),desc(SR_round), desc(SL_round),desc(WK_round))
LOorder<-c(orderLO$fish)


ggplot(indLO, aes(x=factor(fish, levels = LOorder), y=prop, fill=source))+
  geom_bar(position = "fill", stat = "identity")+
  theme_void()+
  scale_fill_manual(values = c("#eff3ff","#bdd7e7","#6baed6","#2171b5","#feedde","#fdbe85","#fd8d3c","#d94701"))+
  xlab("")+
  ylab("Proportion")+
  ggtitle("Loon Lake")+
  guides(fill=guide_legend(title=""))+
  theme(legend.position = "left",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(size = 28, hjust = 0.5, color = "darkorange2"))

#### G LAKE ####

## data
g<-read.csv("./02_tidy-data/g-ancestry.csv", stringsAsFactors = T)
g_props<-g %>% summarise(FG_prop=mean(FG),
                         TL_prop=mean(TL),
                         WB_prop=mean(WB),
                         WT_prop=mean(WT),
                         SL_prop=mean(SL),
                         SR_prop=mean(SR),
                         WK_prop=mean(WK))

## stocking
g_stocks<-read.csv("./02_tidy-data/g-stocking.csv", stringsAsFactors = T) %>% 
  mutate(prop0 = stock/sum(g_stocks$stock),
         prop1 = t(g_props),
         delta = ((prop1-prop0)/prop0)*100,
         source = factor(source, levels = c("FG","TL","WB","WT","SL","SR","WK")))

## graphing success
ggplot(g_stocks, aes(x=source, y=delta))+
  geom_hline(yintercept = 68, linetype="dashed")+
  geom_hline(yintercept = -52, linetype="dashed")+
  geom_bar(stat = "identity",
           fill=c("#d1cebe","#a39f81","#777248","#4b4910","#8697a3","#4d6778","#0f3b50"),
           color="black", lwd=0.5)+
  geom_hline(yintercept = 0, color="black")+
  theme_bw()+
  ylab("% Change in Abundance")+
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth =1),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12))+
  annotate(geom = "text", x=0.88, y=84, label="(B)", size=8)

## hlines
(((4/58)-g_stocks[1,3])/g_stocks[1,3])*100
(((14/58)-g_stocks[1,3])/g_stocks[1,3])*100



## simulating
g_sim<-function(stocks){
  pop_i<-c(rep("FG",stocks[1,2]),
           rep("TL",stocks[2,2]),
           rep("WB",stocks[3,2]),
           rep("WT",stocks[4,2]),
           rep("SL",stocks[5,2]),
           rep("SR",stocks[6,2]),
           rep("WK",stocks[7,2]))
  surviving<-round(runif(1, min = 0.1, max=1)*sum(stocks$stock))
  pop_f<-sample(pop_i, surviving, replace = F)
  pop_sample<-sample(pop_f, 58, replace = F)
  counts<-c(sum(pop_sample == "FG"),
            sum(pop_sample == "TL"),
            sum(pop_sample == "WB"),
            sum(pop_sample == "WT"),
            sum(pop_sample == "SL"),
            sum(pop_sample == "SR"),
            sum(pop_sample == "WK"))
  counts
}

gsims<-data.frame(t(replicate(1000, g_sim(g_stocks))))
colnames(gsims)<-c(g_stocks$source)

FG_lower<-quantile(gsims$FG, 0.025)
FG_upper<-quantile(gsims$FG, 0.975)

TL_lower<-quantile(gsims$TL, 0.025)
TL_upper<-quantile(gsims$TL, 0.975)

WB_lower<-quantile(gsims$WB, 0.025)
WB_upper<-quantile(gsims$WB, 0.975)

WT_lower<-quantile(gsims$WT, 0.025)
WT_upper<-quantile(gsims$WT, 0.975)

SL_lower<-quantile(gsims$SL, 0.025)
SL_upper<-quantile(gsims$SL, 0.975)

SR_lower<-quantile(gsims$SR, 0.025)
SR_upper<-quantile(gsims$SR, 0.975)

WK_lower<-quantile(gsims$WK, 0.025)
WK_upper<-quantile(gsims$WK, 0.975)

g_props*58





