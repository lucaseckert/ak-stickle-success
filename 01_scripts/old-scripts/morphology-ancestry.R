## morphology + ancestry ## 

#### packages and data #### 

## packages
library(tidyverse)
library(GGally)
library(gridExtra)

## loading data
ancestry<-read_csv("./02_tidy-data/final-ancestry.csv")[,-1]
morphology<-read_csv("./02_tidy-data/linear-measurements-2021.csv")
lethal<-read_csv("./02_tidy-data/lethal-2021.csv")
lethal$sex<-factor(lethal$Sex_f_m_u)
lethal$fish<-lethal$Fish_ID

## combining
data<-full_join(morphology,ancestry, by=c("fish","lake"))
data$lake<-factor(data$lake)
data<-left_join(data,
                select(lethal, fish, sex, mass=Mass_g, std_length=Std_length_mm, schisto=S_solidus_N),
                by="fish")

## adjusting
data$adj.body.depth<-data$body.depth/data$st.length
data$adj.eye.diameter<-data$eye.diameter/data$st.length
data$adj.head.length<-data$head.length/data$st.length
data$adj.snout.length<-data$snout.length/data$st.length
data$adj.cp.depth<-data$cp.depth/data$st.length

data %>% select(lake,st.length,adj.body.depth,adj.cp.depth,adj.eye.diameter,adj.head.length,adj.snout.length) %>% 
  ggpairs(aes(color=lake, alpha=0.25))

#### sexual dimorphism ####

ggplot(data, aes(st.length, std_length))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  facet_wrap(~lake)

sexdata<-filter(data, sex!="u")

ggplot(sexdata, aes(sex, adj.body.depth))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~lake)
ggplot(sexdata, aes(sex, adj.eye.diameter))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~lake)
ggplot(sexdata, aes(sex, adj.cp.depth))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~lake)
ggplot(sexdata, aes(sex, adj.head.length))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~lake)
ggplot(sexdata, aes(sex, adj.snout.length))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~lake)

#### benthic hybrids ####

## LP
lp<-data %>% filter(lake == "LP") %>% 
  filter(FG_round == 0, TL_round==0)
summary(lp)

## bd
lp_bd <- ggplot(lp, aes(WB_round, adj.body.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="red", fill="red", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.body.depth ~ WB_round, lp))

## cpd
lp_cpd<- ggplot(lp, aes(WB_round, adj.cp.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="red", fill="red", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.cp.depth ~ WB_round, lp))

# ed
lp_ed <- ggplot(lp, aes(WB_round, adj.eye.diameter))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.eye.diameter ~ WB_round, lp))

## hl
lp_hl<-ggplot(lp, aes(WB_round, adj.head.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.head.length ~ WB_round, lp))

## sn
lp_sn<-ggplot(lp, aes(WB_round, adj.snout.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(panel.grid = element_blank())
summary(lm(adj.snout.length ~ WB_round, lp))

## CC
cc<-data %>% filter(lake == "CC") %>% 
  filter(FG_round == 0, TL_round==0)
summary(cc)

## bd
cc_bd <- ggplot(cc, aes(WB_round, adj.body.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.body.depth ~ WB_round, cc))

## cpd
cc_cpd<- ggplot(cc, aes(WB_round, adj.cp.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.cp.depth ~ WB_round, cc))

# ed
cc_ed <- ggplot(cc, aes(WB_round, adj.eye.diameter))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.eye.diameter ~ WB_round, cc))

## hl
cc_hl<-ggplot(cc, aes(WB_round, adj.head.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.head.length ~ WB_round, cc))

## sn
cc_sn<-ggplot(cc, aes(WB_round, adj.snout.length))+
  geom_point()+
  geom_smooth(method = "lm", color="yellowgreen", fill="yellowgreen", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.snout.length ~ WB_round, cc))

## LL
ll<-data %>% filter(lake == "LL") %>% 
  filter(FG_round == 0, TL_round==0)
summary(ll)

## bd
ll_bd <- ggplot(ll, aes(WB_round, adj.body.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.body.depth ~ WB_round, ll))

## cpd
ll_cpd<- ggplot(ll, aes(WB_round, adj.cp.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.cp.depth ~ WB_round, ll))

# ed
ll_ed <- ggplot(ll, aes(WB_round, adj.eye.diameter))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.eye.diameter ~ WB_round, ll))

## hl
ll_hl<-ggplot(ll, aes(WB_round, adj.head.length))+
  geom_point()+
  geom_smooth(method = "lm", color="yellowgreen", fill="yellowgreen", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.head.length ~ WB_round, ll))

## sn
ll_sn<-ggplot(ll, aes(WB_round, adj.snout.length))+
  geom_point()+
  geom_smooth(method = "lm", color="yellowgreen", fill="yellowgreen", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.snout.length ~ WB_round, ll))

## plotting all
grid.arrange(cc_bd, cc_cpd, cc_ed, cc_sn, cc_hl,
             ll_bd, ll_cpd, ll_ed, ll_sn, ll_hl,
             lp_bd, lp_cpd, lp_ed, lp_sn, lp_hl,
             nrow = 3)

#### limnetic hybrids ####

## CL
cl<-data %>% filter(lake == "CL") %>% 
  filter(WK_round == 0, SL_round==0)
summary(cl)

## bd
cl_bd <- ggplot(cl, aes(LG_round, adj.body.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.body.depth ~ LG_round, cl))

## cpd
cl_cpd<- ggplot(cl, aes(LG_round, adj.cp.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.cp.depth ~ LG_round, cl))

# ed
cl_ed <- ggplot(cl, aes(LG_round, adj.eye.diameter))+
  geom_point()+
  geom_smooth(method = "lm", color="red", fill="red", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.eye.diameter ~ LG_round, cl))

## hl
cl_hl<-ggplot(cl, aes(LG_round, adj.head.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.head.length ~ LG_round, cl))

## sn
cl_sn<-ggplot(cl, aes(LG_round, adj.snout.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(panel.grid = element_blank())
summary(lm(adj.snout.length ~ LG_round, cl))

## FL
fl<-data %>% filter(lake == "FL") %>% 
  filter(SL_round == 0, WK_round==0)
summary(fl)

## bd
fl_bd <- ggplot(fl, aes(LG_round, adj.body.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="yellowgreen", fill="yellowgreen", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.body.depth ~ LG_round, fl))

## cpd
fl_cpd<- ggplot(fl, aes(LG_round, adj.cp.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.cp.depth ~ LG_round, fl))

# ed
fl_ed <- ggplot(fl, aes(LG_round, adj.eye.diameter))+
  geom_point()+
  geom_smooth(method = "lm", color="red", fill="red", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.eye.diameter ~ LG_round, fl))

## hl
fl_hl<-ggplot(fl, aes(LG_round, adj.head.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.head.length ~ LG_round, fl))

## sn
fl_sn<-ggplot(fl, aes(LG_round, adj.snout.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.snout.length ~ LG_round, fl))

## HL
hl<-data %>% filter(lake == "HL") %>% 
  filter(SL_round == 0, WK_round==0)
summary(hl)

## bd
hl_bd <- ggplot(hl, aes(LG_round, adj.body.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.body.depth ~ LG_round, hl))

## cpd
hl_cpd<- ggplot(hl, aes(LG_round, adj.cp.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.cp.depth ~ LG_round, hl))

# ed
hl_ed <- ggplot(hl, aes(LG_round, adj.eye.diameter))+
  geom_point()+
  geom_smooth(method = "lm", color="red", fill="red", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.eye.diameter ~ LG_round, hl))

## hl
hl_hl<-ggplot(hl, aes(LG_round, adj.head.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.head.length ~ LG_round, hl))

## sn
hl_sn<-ggplot(hl, aes(LG_round, adj.snout.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.snout.length ~ LG_round, hl))

## RL
rl<-data %>% filter(lake == "RL") %>% 
  filter(SL_round == 0, WK_round==0)
summary(rl)

## bd
rl_bd <- ggplot(rl, aes(LG_round, adj.body.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.body.depth ~ LG_round, rl))

## cpd
rl_cpd<- ggplot(rl, aes(LG_round, adj.cp.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.cp.depth ~ LG_round, rl))

# ed
rl_ed <- ggplot(rl, aes(LG_round, adj.eye.diameter))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.eye.diameter ~ LG_round, rl))

## hl
rl_hl<-ggplot(rl, aes(LG_round, adj.head.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.head.length ~ LG_round, rl))

## sn
rl_sn<-ggplot(rl, aes(LG_round, adj.snout.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.snout.length ~ LG_round, rl))

## plotting all
grid.arrange(cl_bd, cl_cpd, cl_ed, cl_sn, cl_hl,
             fl_bd, fl_cpd, fl_ed, fl_sn, fl_hl,
             hl_bd, hl_cpd, hl_ed, hl_sn, hl_hl,
             rl_bd, rl_cpd, rl_ed, rl_sn, rl_hl,
             nrow = 4)


#### loon hybrids ####
loon <- data %>% filter(lake=="LO") %>% 
  mutate(benthic_prop = FG_round+TL_round+WB_round+WT_round)

loon_bd<- ggplot(loon, aes(benthic_prop, adj.body.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="yellowgreen", fill="yellowgreen", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.body.depth~benthic_prop, loon))

loon_cpd<- ggplot(loon, aes(benthic_prop, adj.cp.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="yellowgreen", fill="yellowgreen", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.cp.depth~benthic_prop, loon))

loon_ed<- ggplot(loon, aes(benthic_prop, adj.eye.diameter))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.eye.diameter~benthic_prop, loon))

loon_sn<- ggplot(loon, aes(benthic_prop, adj.snout.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.snout.length~benthic_prop, loon))

loon_hl<- ggplot(loon, aes(benthic_prop, adj.head.length))+
  geom_point()+
  geom_smooth(method = "lm", color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank())
summary(lm(adj.head.length~benthic_prop, loon))

grid.arrange(loon_bd, loon_cpd, loon_ed, loon_sn, loon_hl, nrow=1)

#### individual ancestry graphs ####

## plotting individuals
ind1<- data %>% filter(lake=="LP") %>% 
  select(fish, lake, FG_round, TL_round, WB_round, WT_round, LG_round, SL_round, SR_round, WK_round)
ind2<- ind1 %>% gather(source, prop, FG_round:WK_round) %>%
  mutate(source=substr(source,1,2), fish=as.factor(fish))

ind2$source<-factor(ind$source, levels = c("FG","TL","WB","WT","LG","SL","SR","WK"))

ggplot(ind2, aes(x=factor(fish), y=prop, fill=source))+
  geom_bar(position = "fill", stat = "identity")+
  theme_void()+
  scale_fill_manual(values = c("#feedde","#fdbe85","#fd8d3c","#d94701",
                               "#eff3ff","#bdd7e7","#6baed6","#2171b5"))+
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

pureFG<-filter(ind1,TL_round==0,WB_round==0,WT_round==0)
nrow(pureFG)

pureTL<-filter(ind1,FG_round==0,WB_round==0,WT_round==0)
nrow(pureTL)

pureWB<-filter(ind1,TL_round==0,FG_round==0,WT_round==0)
nrow(pureWB)

pureWT<-filter(ind1,TL_round==0,WB_round==0,FG_round==0)
nrow(pureWT)

#### GxE ####

## WB
wb<-data %>% filter(WB_round ==1)
ggplot(wb, aes(lake,adj.body.depth))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

## WT
wt<-data %>% filter(WT_round ==1)
ggplot(wt, aes(lake,adj.body.depth))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

## FG
fg<-data %>% filter(FG_round ==1)
ggplot(fg, aes(lake,adj.body.depth))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

## LG
lg<-data %>% filter(LG_round ==1)
ggplot(lg, aes(lake,adj.body.depth))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

## SR
sr<-data %>% filter(SR_round ==1)
ggplot(sr, aes(lake,adj.body.depth))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

## not enough pure SL or WK

## GxE in body depth between CC and LL
fg$source<-"FG"
wb$source<-"WB"
wt$source<-"WT"
ccll <- rbind(fg,wb,wt) %>% 
  filter(lake %in% c("CC","LL"))%>%
  mutate(lake2 = case_when(lake == "CC" ~ 0,
                           lake == "LL" ~ 1))

ggplot(ccll, aes(x=lake, y=adj.body.depth))+
  geom_boxplot(color="grey60")+
  geom_point()+
  theme_bw()+
  xlab("")+
  ylab("Adj. Body Depth")+
  theme(panel.grid = element_blank())+
  facet_grid(~source)

ggplot(ccll, aes(x=lake2, y=adj.body.depth))+
  geom_point()+
  geom_smooth(method = "lm", color="purple", fill="purple", alpha=0.25 )+
  theme_bw()+
  xlab("")+
  ylab("Adj. Body Depth")+
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12))+
  scale_x_continuous(breaks=c(0,1),
                   labels=c("CC", "LL"),
                   limits = c(-0.25,1.25))+
  facet_grid(~source)

ggplot(ccll, aes(x=lake2, y=adj.eye.diameter))+
  geom_point()+
  geom_smooth(method = "lm", color="purple", fill="purple", alpha=0.25 )+
  theme_bw()+
  xlab("")+
  ylab("Adj. Eye Diameter")+
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12))+
  scale_x_continuous(breaks=c(0,1),
                     labels=c("CC", "LL"),
                     limits = c(-0.25,1.25))+
  facet_grid(~source)
  
ggplot(ccll, aes(x=lake2, y=st.length))+
  geom_point()+
  geom_smooth(method = "lm", color="purple", fill="purple", alpha=0.25 )+
  theme_bw()+
  xlab("")+
  ylab("Standard Length")+
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=12))+
  scale_x_continuous(breaks=c(0,1),
                     labels=c("CC", "LL"),
                     limits = c(-0.25,1.25))+
  facet_grid(~source)

#### isotope data ####

isotopes<-read_csv("./02_tidy-data/isotopes-2021.csv") %>% 
  select(fish=SampleID, N15, C13, O18)

data<- left_join(data, isotopes, by="fish")
summary(data)

select(data, lake, adj.body.depth, adj.eye.diameter, N15, C13, O18) %>% 
  ggpairs(aes(color=lake))

ggplot(data, aes(adj.body.depth, N15))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

ggplot(data, aes(adj.body.depth, C13))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

ggplot(data, aes(adj.body.depth, O18))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

ggplot(data, aes(LG_round, N15))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

ggplot(data, aes(WB_round, N15))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

## PCA for ancestry

## benthic PCA
benthics<-c("CC","LP","LL")
benthic<-filter(data, lake %in% benthics)
pcab<-princomp(~WB_round+WT_round+FG_round+TL_round,
               data = benthic, cor = T, scores = T)
benthic.pca<-predict(pcab, newdata = benthic)

benthic<-cbind(benthic, benthic.pca)
ggplot(benthic, aes(Comp.1, N15))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

ggplot(benthic, aes(Comp.2, N15))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

## limnetic PCA
limnetics<-c("CL","FL","HL","RL")
limnetic<-filter(data, lake %in% limnetics)
pcal<-princomp(~WK_round+LG_round+SR_round+SL_round,
               data = limnetic, cor = T, scores = T)
limnetic.pca<-predict(pcal, newdata = limnetic)

limnetic<-cbind(limnetic, limnetic.pca)
ggplot(limnetic, aes(Comp.1, N15))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

ggplot(limnetic, aes(Comp.2, N15))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

## loon PCA
pcalo<-princomp(~WB_round+WT_round+FG_round+TL_round+LG_round+SR_round+SL_round,
               data = loon, cor = T, scores = T)
loon.pca<-predict(pcalo, newdata = loon)

loon<-cbind(loon, loon.pca)
ggplot(loon, aes(Comp.1, C13))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

ggplot(loon, aes(Comp.2, N15))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  facet_grid(~lake)

ggplot(loon, aes(benthic_prop, N15))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()


#### benthic hybrids + sex ####

## LP
lp<-sexdata %>% filter(lake == "LP") %>% 
  filter(FG_round == 0, TL_round==0)
summary(lp)

## bd
lp_bd <- ggplot(lp, aes(WB_round, adj.body.depth, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.body.depth ~ WB_round, lp))

## cpd
lp_cpd<- ggplot(lp, aes(WB_round, adj.cp.depth, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.cp.depth ~ WB_round, lp))

# ed
lp_ed <- ggplot(lp, aes(WB_round, adj.eye.diameter, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.eye.diameter ~ WB_round, lp))

## hl
lp_hl<-ggplot(lp, aes(WB_round, adj.head.length, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.head.length ~ WB_round, lp))

## sn
lp_sn<-ggplot(lp, aes(WB_round, adj.snout.length, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.snout.length ~ WB_round, lp))

## CC
cc<-sexdata %>% filter(lake == "CC") %>% 
  filter(FG_round == 0, TL_round==0)
summary(cc)

## bd
cc_bd <- ggplot(cc, aes(WB_round, adj.body.depth, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.body.depth ~ WB_round, cc))

## cpd
cc_cpd<- ggplot(cc, aes(WB_round, adj.cp.depth, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.cp.depth ~ WB_round, cc))

# ed
cc_ed <- ggplot(cc, aes(WB_round, adj.eye.diameter, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.eye.diameter ~ WB_round, cc))

## hl
cc_hl<-ggplot(cc, aes(WB_round, adj.head.length, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.head.length ~ WB_round, cc))

## sn
cc_sn<-ggplot(cc, aes(WB_round, adj.snout.length, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.snout.length ~ WB_round, cc))

## LL
ll<-sexdata %>% filter(lake == "LL") %>% 
  filter(FG_round == 0, TL_round==0)
summary(ll)

## bd
ll_bd <- ggplot(ll, aes(WB_round, adj.body.depth, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(plot.title = element_text(hjust=0.5, size=16),
        axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.body.depth ~ WB_round, ll))

## cpd
ll_cpd<- ggplot(ll, aes(WB_round, adj.cp.depth, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.cp.depth ~ WB_round, ll))

# ed
ll_ed <- ggplot(ll, aes(WB_round, adj.eye.diameter, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.eye.diameter ~ WB_round, ll))

## hl
ll_hl<-ggplot(ll, aes(WB_round, adj.head.length, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.head.length ~ WB_round, ll))

## sn
ll_sn<-ggplot(ll, aes(WB_round, adj.snout.length, color=sex))+
  geom_point()+
  geom_smooth(method = "lm", alpha=0.25, aes(group=sex))+
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.title = element_text(size=12),
        panel.grid = element_blank(),
        legend.position = "none")
summary(lm(adj.snout.length ~ WB_round, ll))

## plotting all
grid.arrange(cc_bd, cc_cpd, cc_ed, cc_sn, cc_hl,
             ll_bd, ll_cpd, ll_ed, ll_sn, ll_hl,
             lp_bd, lp_cpd, lp_ed, lp_sn, lp_hl,
             nrow = 3)
