## simulations ##

## packages
library(tidyverse)

#### for simulations ####
a<-matrix(rep(0.25,16), nrow = 4)

## loading data
data<-read_csv("./02_tidy-data/final-ancestry.csv")[,-1]
data$lake<-factor(data$lake)

## CC F1
CCF1<- data %>% filter(lake=="CC") %>% 
  select(fish, lake, FG_round, TL_round, WB_round, WT_round, LG_round, SL_round, SR_round, WK_round) %>% 
  filter(FG_round!=0.25,TL_round!=0.25,WB_round!=0.25,WT_round!=0.25,)
nrow(CCF1)

CC_FG<-sum(CCF1$FG_round)/nrow(CCF1)
CC_TL<-sum(CCF1$TL_round)/nrow(CCF1)
CC_WB<-sum(CCF1$WB_round)/nrow(CCF1)
CC_WT<-sum(CCF1$WT_round)/nrow(CCF1)

pCC<-c(CC_FG,CC_TL,CC_WB,CC_WT)
mCC<-a*pCC

CC_AA<-mCC[1,1]
CC_AB<-mCC[1,2]+mCC[2,1]
CC_AC<-mCC[1,3]+mCC[3,1]
CC_AD<-mCC[1,4]+mCC[4,1]
CC_BB<-mCC[2,2]
CC_BC<-mCC[2,3]+mCC[3,2]
CC_BD<-mCC[2,4]+mCC[4,2]
CC_CC<-mCC[3,3]
CC_CD<-mCC[3,4]+mCC[4,3]
CC_DD<-mCC[4,4]

CC_AA<-CC_FG*CC_FG
CC_AB<-CC_FG*CC_TL*2
CC_AC<-CC_FG*CC_WB*2
CC_AD<-CC_FG*CC_WT*2
CC_BB<-CC_TL*CC_TL
CC_BC<-CC_TL*CC_WB*2
CC_BD<-CC_TL*CC_WT*2
CC_CC<-CC_WB*CC_WB
CC_CD<-CC_WB*CC_WT*2
CC_DD<-CC_WT*CC_WT

p2CC<-c(CC_AA,CC_AB,CC_AC,CC_AD,CC_BB,CC_BC,CC_BD,CC_CC,CC_CD,CC_DD)
eCC<-round(p2CC*nrow(CCF1),0)
names(eCC)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
eCC

oCC_FGp<-nrow(filter(CCF1,FG_round==1,TL_round==0,WB_round==0,WT_round==0))
oCC_FGxTL<-nrow(filter(CCF1,WB_round==0,WT_round==0,FG_round!=1,TL_round!=1))
oCC_FGxWB<-nrow(filter(CCF1,FG_round!=1,TL_round==0,WB_round!=1,WT_round==0))
oCC_FGxWT<-nrow(filter(CCF1,FG_round!=1,TL_round==0,WB_round==0,WT_round!=1))
oCC_TLp<-nrow(filter(CCF1,FG_round==0,TL_round==1,WB_round==0,WT_round==0))
oCC_TLxWB<-nrow(filter(CCF1,FG_round==0,TL_round!=1,WB_round!=1,WT_round==0))
oCC_TLxWT<-nrow(filter(CCF1,FG_round==0,TL_round!=1,WB_round==0,WT_round!=1))
oCC_WBp<-nrow(filter(CCF1,FG_round==0,TL_round==0,WB_round==1,WT_round==0))
oCC_WBxWT<-nrow(filter(CCF1,FG_round==0,TL_round==0,WB_round!=1,WT_round!=1))
oCC_WTp<-nrow(filter(CCF1,FG_round==0,TL_round==0,WB_round==0,WT_round==1))

oCC<-c(oCC_FGp,oCC_FGxTL,oCC_FGxWB,oCC_FGxWT,oCC_TLp,oCC_TLxWB,oCC_TLxWT,oCC_WBp,oCC_WBxWT,oCC_WTp)
names(oCC)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
oCC

rCC<-oCC-eCC
rCC

pCC<-rCC*100/eCC

CC_success<-c(CC_FG, (CC_FG+CC_TL)/2, (CC_FG+CC_WB)/2, (CC_FG+CC_WT)/2, CC_TL, (CC_TL+CC_WB)/2, (CC_TL+CC_WT)/2, CC_WB, (CC_WB+CC_WT)/2, CC_WT)
names(CC_success)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
CC_success

CC_success2<-c(CC_FG*CC_FG, CC_FG*CC_TL, CC_FG*CC_WB, CC_FG*CC_WT, CC_TL*CC_TL, CC_TL*CC_WB, CC_TL*CC_WT, CC_WB*CC_WB, CC_WB*CC_WT, CC_WT*CC_WT)
names(CC_success2)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
CC_success2

CC_df<-data.frame(lake="CC", hybrid=names(rCC), ex=eCC, ob=oCC, residual=rCC, success=CC_success, success2=CC_success2)

ggplot(CC_df, aes(x=success, y=residual))+
  geom_point()+
  theme_bw()

## LL F1
LLF1<- data %>% filter(lake=="LL") %>% 
  select(fish, lake, FG_round, TL_round, WB_round, WT_round, LG_round, SL_round, SR_round, WK_round) %>% 
  filter(FG_round!=0.25,TL_round!=0.25,WB_round!=0.25,WT_round!=0.25,)
nrow(LLF1)

LL_FG<-sum(LLF1$FG_round)/nrow(LLF1)
LL_TL<-sum(LLF1$TL_round)/nrow(LLF1)
LL_WB<-sum(LLF1$WB_round)/nrow(LLF1)
LL_WT<-sum(LLF1$WT_round)/nrow(LLF1)

pLL<-c(LL_FG,LL_TL,LL_WB,LL_WT)
mLL<-a*pLL

LL_AA<-mLL[1,1]
LL_AB<-mLL[1,2]+mLL[2,1]
LL_AC<-mLL[1,3]+mLL[3,1]
LL_AD<-mLL[1,4]+mLL[4,1]
LL_BB<-mLL[2,2]
LL_BC<-mLL[2,3]+mLL[3,2]
LL_BD<-mLL[2,4]+mLL[4,2]
LL_CC<-mLL[3,3]
LL_CD<-mLL[3,4]+mLL[4,3]
LL_DD<-mLL[4,4]

LL_AA<-LL_FG*LL_FG
LL_AB<-LL_FG*LL_TL*2
LL_AC<-LL_FG*LL_WB*2
LL_AD<-LL_FG*LL_WT*2
LL_BB<-LL_TL*LL_TL
LL_BC<-LL_TL*LL_WB*2
LL_BD<-LL_TL*LL_WT*2
LL_CC<-LL_WB*LL_WB
LL_CD<-LL_WB*LL_WT*2
LL_DD<-LL_WT*LL_WT

p2LL<-c(LL_AA,LL_AB,LL_AC,LL_AD,LL_BB,LL_BC,LL_BD,LL_CC,LL_CD,LL_DD)
eLL<-round(p2LL*nrow(LLF1),0)
names(eLL)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
eLL

oLL_FGp<-nrow(filter(LLF1,FG_round==1,TL_round==0,WB_round==0,WT_round==0))
oLL_FGxTL<-nrow(filter(LLF1,WB_round==0,WT_round==0,FG_round!=1,TL_round!=1))
oLL_FGxWB<-nrow(filter(LLF1,FG_round!=1,TL_round==0,WB_round!=1,WT_round==0))
oLL_FGxWT<-nrow(filter(LLF1,FG_round!=1,TL_round==0,WB_round==0,WT_round!=1))
oLL_TLp<-nrow(filter(LLF1,FG_round==0,TL_round==1,WB_round==0,WT_round==0))
oLL_TLxWB<-nrow(filter(LLF1,FG_round==0,TL_round!=1,WB_round!=1,WT_round==0))
oLL_TLxWT<-nrow(filter(LLF1,FG_round==0,TL_round!=1,WB_round==0,WT_round!=1))
oLL_WBp<-nrow(filter(LLF1,FG_round==0,TL_round==0,WB_round==1,WT_round==0))
oLL_WBxWT<-nrow(filter(LLF1,FG_round==0,TL_round==0,WB_round!=1,WT_round!=1))
oLL_WTp<-nrow(filter(LLF1,FG_round==0,TL_round==0,WB_round==0,WT_round==1))

oLL<-c(oLL_FGp,oLL_FGxTL,oLL_FGxWB,oLL_FGxWT,oLL_TLp,oLL_TLxWB,oLL_TLxWT,oLL_WBp,oLL_WBxWT,oLL_WTp)
names(oLL)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
oLL

rLL<-oLL-eLL
rLL

LL_success<-c(LL_FG, (LL_FG+LL_TL)/2, (LL_FG+LL_WB)/2, (LL_FG+LL_WT)/2, LL_TL, (LL_TL+LL_WB)/2, (LL_TL+LL_WT)/2, LL_WB, (LL_WB+LL_WT)/2, LL_WT)
names(LL_success)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
LL_success

LL_success2<-c(LL_FG*LL_FG, LL_FG*LL_TL, LL_FG*LL_WB, LL_FG*LL_WT, LL_TL*LL_TL, LL_TL*LL_WB, LL_TL*LL_WT, LL_WB*LL_WB, LL_WB*LL_WT, LL_WT*LL_WT)
names(LL_success2)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
LL_success2

LL_df<-data.frame(lake="LL", hybrid=names(rLL), ex=eLL, ob=oLL, residual=rLL, success=LL_success, success2=LL_success2)

ggplot(LL_df, aes(x=success, y=residual))+
  geom_point()+
  theme_bw()

## LP F1
LPF1<- data %>% filter(lake=="LP") %>% 
  select(fish, lake, FG_round, TL_round, WB_round, WT_round, LG_round, SL_round, SR_round, WK_round) %>% 
  filter(FG_round!=0.25,TL_round!=0.25,WB_round!=0.25,WT_round!=0.25,)
nrow(LPF1)

LP_FG<-sum(LPF1$FG_round)/nrow(LPF1)
LP_TL<-sum(LPF1$TL_round)/nrow(LPF1)
LP_WB<-sum(LPF1$WB_round)/nrow(LPF1)
LP_WT<-sum(LPF1$WT_round)/nrow(LPF1)

pLP<-c(LP_FG,LP_TL,LP_WB,LP_WT)
mLP<-a*pLP

LP_AA<-mLP[1,1]
LP_AB<-mLP[1,2]+mLP[2,1]
LP_AC<-mLP[1,3]+mLP[3,1]
LP_AD<-mLP[1,4]+mLP[4,1]
LP_BB<-mLP[2,2]
LP_BC<-mLP[2,3]+mLP[3,2]
LP_BD<-mLP[2,4]+mLP[4,2]
LP_CC<-mLP[3,3]
LP_CD<-mLP[3,4]+mLP[4,3]
LP_DD<-mLP[4,4]

LP_AA<-LP_FG*LP_FG
LP_AB<-LP_FG*LP_TL*2
LP_AC<-LP_FG*LP_WB*2
LP_AD<-LP_FG*LP_WT*2
LP_BB<-LP_TL*LP_TL
LP_BC<-LP_TL*LP_WB*2
LP_BD<-LP_TL*LP_WT*2
LP_CC<-LP_WB*LP_WB
LP_CD<-LP_WB*LP_WT*2
LP_DD<-LP_WT*LP_WT

p2LP<-c(LP_AA,LP_AB,LP_AC,LP_AD,LP_BB,LP_BC,LP_BD,LP_CC,LP_CD,LP_DD)
eLP<-round(p2LP*nrow(LPF1),0)
names(eLP)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
eLP

oLP_FGp<-nrow(filter(LPF1,TL_round==0,WB_round==0,WT_round==0))
oLP_FGxTL<-nrow(filter(LPF1,WB_round==0,WT_round==0,FG_round!=1,TL_round!=1))
oLP_FGxWB<-nrow(filter(LPF1,FG_round!=1,TL_round==0,WB_round!=1,WT_round==0))
oLP_FGxWT<-nrow(filter(LPF1,FG_round!=1,TL_round==0,WB_round==0,WT_round!=1))
oLP_TLp<-nrow(filter(LPF1,FG_round==1,WB_round==0,WT_round==0))
oLP_TLxWB<-nrow(filter(LPF1,FG_round==0,TL_round!=1,WB_round!=1,WT_round==0))
oLP_TLxWT<-nrow(filter(LPF1,FG_round==0,TL_round!=1,WB_round==0,WT_round!=1))
oLP_WBp<-nrow(filter(LPF1,FG_round==0,TL_round==0,WT_round==0))
oLP_WBxWT<-nrow(filter(LPF1,FG_round==0,TL_round==0,WB_round!=1,WT_round!=1))
oLP_WTp<-nrow(filter(LPF1,FG_round==0,TL_round==0,WB_round==0))

oLP<-c(oLP_FGp,oLP_FGxTL,oLP_FGxWB,oLP_FGxWT,oLP_TLp,oLP_TLxWB,oLP_TLxWT,oLP_WBp,oLP_WBxWT,oLP_WTp)
names(oLP)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
oLP

rLP<-oLP-eLP
rLP

LP_success<-c(LP_FG, (LP_FG+LP_TL)/2, (LP_FG+LP_WB)/2, (LP_FG+LP_WT)/2, LP_TL, (LP_TL+LP_WB)/2, (LP_TL+LP_WT)/2, LP_WB, (LP_WB+LP_WT)/2, LP_WT)
names(LP_success)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
LP_success

LP_success2<-c(LP_FG*LP_FG, LP_FG*LP_TL, LP_FG*LP_WB, LP_FG*LP_WT, LP_TL*LP_TL, LP_TL*LP_WB, LP_TL*LP_WT, LP_WB*LP_WB, LP_WB*LP_WT, LP_WT*LP_WT)
names(LP_success2)<-c("FGp","FGxTL","FGxWB","FGxWT","TLp","TLxWB","TLxWT","WBp","WBxWT","WTp")
LP_success2

LP_df<-data.frame(lake="LP", hybrid=names(rLP), ex=eLP, ob=oLP, residual=rLP, success=LP_success, success2=LP_success2)

ggplot(LP_df, aes(x=success, y=residual))+
  geom_point()+
  theme_bw()

## CL F1
CLF1<- data %>% filter(lake=="CL") %>% 
  select(fish, lake, FG_round, TL_round, WB_round, WT_round, LG_round, SL_round, SR_round, WK_round) %>% 
  filter(LG_round!=0.25,SL_round!=0.25,SR_round!=0.25,WK_round!=0.25,)
nrow(CLF1)

CL_LG<-sum(CLF1$LG_round)/nrow(CLF1)
CL_SL<-sum(CLF1$SL_round)/nrow(CLF1)
CL_SR<-sum(CLF1$SR_round)/nrow(CLF1)
CL_WK<-sum(CLF1$WK_round)/nrow(CLF1)

pCL<-c(CL_LG,CL_SL,CL_SR,CL_WK)
mCL<-a*pCL

CL_AA<-mCL[1,1]
CL_AB<-mCL[1,2]+mCL[2,1]
CL_AC<-mCL[1,3]+mCL[3,1]
CL_AD<-mCL[1,4]+mCL[4,1]
CL_BB<-mCL[2,2]
CL_BC<-mCL[2,3]+mCL[3,2]
CL_BD<-mCL[2,4]+mCL[4,2]
CL_CC<-mCL[3,3]
CL_CD<-mCL[3,4]+mCL[4,3]
CL_DD<-mCL[4,4]

CL_AA<-CL_LG*CL_LG
CL_AB<-CL_LG*CL_SL*2
CL_AC<-CL_LG*CL_SR*2
CL_AD<-CL_LG*CL_WK*2
CL_BB<-CL_SL*CL_SL
CL_BC<-CL_SL*CL_SR*2
CL_BD<-CL_SL*CL_WK*2
CL_CC<-CL_SR*CL_SR
CL_CD<-CL_SR*CL_WK*2
CL_DD<-CL_WK*CL_WK

p2CL<-c(CL_AA,CL_AB,CL_AC,CL_AD,CL_BB,CL_BC,CL_BD,CL_CC,CL_CD,CL_DD)
eCL<-round(p2CL*nrow(CLF1),0)
names(eCL)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
eCL

oCL_LGp<-nrow(filter(CLF1,LG_round==1,SL_round==0,SR_round==0,WK_round==0))
oCL_LGxSL<-nrow(filter(CLF1,SR_round==0,WK_round==0,LG_round!=1,SL_round!=1))
oCL_LGxSR<-nrow(filter(CLF1,LG_round!=1,SL_round==0,SR_round!=1,WK_round==0))
oCL_LGxWK<-nrow(filter(CLF1,LG_round!=1,SL_round==0,SR_round==0,WK_round!=1))
oCL_SLp<-nrow(filter(CLF1,LG_round==0,SL_round==1,SR_round==0,WK_round==0))
oCL_SLxSR<-nrow(filter(CLF1,LG_round==0,SL_round!=1,SR_round!=1,WK_round==0))
oCL_SLxWK<-nrow(filter(CLF1,LG_round==0,SL_round!=1,SR_round==0,WK_round!=1))
oCL_SRp<-nrow(filter(CLF1,LG_round==0,SL_round==0,SR_round==1,WK_round==0))
oCL_SRxWK<-nrow(filter(CLF1,LG_round==0,SL_round==0,SR_round!=1,WK_round!=1))
oCL_WKp<-nrow(filter(CLF1,LG_round==0,SL_round==0,SR_round==0,WK_round==1))

oCL<-c(oCL_LGp,oCL_LGxSL,oCL_LGxSR,oCL_LGxWK,oCL_SLp,oCL_SLxSR,oCL_SLxWK,oCL_SRp,oCL_SRxWK,oCL_WKp)
names(oCL)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
oCL

rCL<-oCL-eCL
rCL

CL_success<-c(CL_LG, (CL_LG+CL_SL)/2, (CL_LG+CL_SR)/2, (CL_LG+CL_WK)/2, CL_SL, (CL_SL+CL_SR)/2, (CL_SL+CL_WK)/2, CL_SR, (CL_SR+CL_WK)/2, CL_WK)
names(CL_success)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
CL_success

CL_success2<-c(CL_LG*CL_LG, CL_LG*CL_SL, CL_LG*CL_SR, CL_LG*CL_WK, CL_SL*CL_SL, CL_SL*CL_SR, CL_SL*CL_WK, CL_SR*CL_SR, CL_SR*CL_WK, CL_WK*CL_WK)
names(CL_success2)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
CL_success2

CL_df<-data.frame(lake="CL", hybrid=names(rCL), ex=eCL, ob=oCL, residual=rCL, success=CL_success, success2=CL_success2)

ggplot(CL_df, aes(x=success, y=residual))+
  geom_point()+
  theme_bw()

## FL F1
FLF1<- data %>% filter(lake=="FL") %>% 
  select(fish, lake, FG_round, TL_round, WB_round, WT_round, LG_round, SL_round, SR_round, WK_round) %>% 
  filter(LG_round!=0.25,SL_round!=0.25,SR_round!=0.25,WK_round!=0.25,)
nrow(FLF1)

FL_LG<-sum(FLF1$LG_round)/nrow(FLF1)
FL_SL<-sum(FLF1$SL_round)/nrow(FLF1)
FL_SR<-sum(FLF1$SR_round)/nrow(FLF1)
FL_WK<-sum(FLF1$WK_round)/nrow(FLF1)

pFL<-c(FL_LG,FL_SL,FL_SR,FL_WK)
mFL<-a*pFL

FL_AA<-mFL[1,1]
FL_AB<-mFL[1,2]+mFL[2,1]
FL_AC<-mFL[1,3]+mFL[3,1]
FL_AD<-mFL[1,4]+mFL[4,1]
FL_BB<-mFL[2,2]
FL_BC<-mFL[2,3]+mFL[3,2]
FL_BD<-mFL[2,4]+mFL[4,2]
FL_CC<-mFL[3,3]
FL_CD<-mFL[3,4]+mFL[4,3]
FL_DD<-mFL[4,4]

FL_AA<-FL_LG*FL_LG
FL_AB<-FL_LG*FL_SL*2
FL_AC<-FL_LG*FL_SR*2
FL_AD<-FL_LG*FL_WK*2
FL_BB<-FL_SL*FL_SL
FL_BC<-FL_SL*FL_SR*2
FL_BD<-FL_SL*FL_WK*2
FL_CC<-FL_SR*FL_SR
FL_CD<-FL_SR*FL_WK*2
FL_DD<-FL_WK*FL_WK

p2FL<-c(FL_AA,FL_AB,FL_AC,FL_AD,FL_BB,FL_BC,FL_BD,FL_CC,FL_CD,FL_DD)
eFL<-round(p2FL*nrow(FLF1),0)
names(eFL)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
eFL

oFL_LGp<-nrow(filter(FLF1,LG_round==1,SL_round==0,SR_round==0,WK_round==0))
oFL_LGxSL<-nrow(filter(FLF1,SR_round==0,WK_round==0,LG_round!=1,SL_round!=1))
oFL_LGxSR<-nrow(filter(FLF1,LG_round!=1,SL_round==0,SR_round!=1,WK_round==0))
oFL_LGxWK<-nrow(filter(FLF1,LG_round!=1,SL_round==0,SR_round==0,WK_round!=1))
oFL_SLp<-nrow(filter(FLF1,LG_round==0,SL_round==1,SR_round==0,WK_round==0))
oFL_SLxSR<-nrow(filter(FLF1,LG_round==0,SL_round!=1,SR_round!=1,WK_round==0))
oFL_SLxWK<-nrow(filter(FLF1,LG_round==0,SL_round!=1,SR_round==0,WK_round!=1))
oFL_SRp<-nrow(filter(FLF1,LG_round==0,SL_round==0,SR_round==1,WK_round==0))
oFL_SRxWK<-nrow(filter(FLF1,LG_round==0,SL_round==0,SR_round!=1,WK_round!=1))
oFL_WKp<-nrow(filter(FLF1,LG_round==0,SL_round==0,SR_round==0,WK_round==1))

oFL<-c(oFL_LGp,oFL_LGxSL,oFL_LGxSR,oFL_LGxWK,oFL_SLp,oFL_SLxSR,oFL_SLxWK,oFL_SRp,oFL_SRxWK,oFL_WKp)
names(oFL)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
oFL

rFL<-oFL-eFL
rFL

FL_success<-c(FL_LG, (FL_LG+FL_SL)/2, (FL_LG+FL_SR)/2, (FL_LG+FL_WK)/2, FL_SL, (FL_SL+FL_SR)/2, (FL_SL+FL_WK)/2, FL_SR, (FL_SR+FL_WK)/2, FL_WK)
names(FL_success)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
FL_success

FL_success2<-c(FL_LG*FL_LG, FL_LG*FL_SL, FL_LG*FL_SR, FL_LG*FL_WK, FL_SL*FL_SL, FL_SL*FL_SR, FL_SL*FL_WK, FL_SR*FL_SR, FL_SR*FL_WK, FL_WK*FL_WK)
names(FL_success2)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
FL_success2

FL_df<-data.frame(lake="FL", hybrid=names(rFL), ex=eFL, ob=oFL, residual=rFL, success=FL_success, success2=FL_success2)

ggplot(FL_df, aes(x=success, y=residual))+
  geom_point()+
  theme_bw()

## HL F1
HLF1<- data %>% filter(lake=="HL") %>% 
  select(fish, lake, FG_round, TL_round, WB_round, WT_round, LG_round, SL_round, SR_round, WK_round) %>% 
  filter(LG_round!=0.25,SL_round!=0.25,SR_round!=0.25,WK_round!=0.25,)
nrow(HLF1)

HL_LG<-sum(HLF1$LG_round)/nrow(HLF1)
HL_SL<-sum(HLF1$SL_round)/nrow(HLF1)
HL_SR<-sum(HLF1$SR_round)/nrow(HLF1)
HL_WK<-sum(HLF1$WK_round)/nrow(HLF1)

pHL<-c(HL_LG,HL_SL,HL_SR,HL_WK)
mHL<-a*pHL

HL_AA<-mHL[1,1]
HL_AB<-mHL[1,2]+mHL[2,1]
HL_AC<-mHL[1,3]+mHL[3,1]
HL_AD<-mHL[1,4]+mHL[4,1]
HL_BB<-mHL[2,2]
HL_BC<-mHL[2,3]+mHL[3,2]
HL_BD<-mHL[2,4]+mHL[4,2]
HL_CC<-mHL[3,3]
HL_CD<-mHL[3,4]+mHL[4,3]
HL_DD<-mHL[4,4]

HL_AA<-HL_LG*HL_LG
HL_AB<-HL_LG*HL_SL*2
HL_AC<-HL_LG*HL_SR*2
HL_AD<-HL_LG*HL_WK*2
HL_BB<-HL_SL*HL_SL
HL_BC<-HL_SL*HL_SR*2
HL_BD<-HL_SL*HL_WK*2
HL_CC<-HL_SR*HL_SR
HL_CD<-HL_SR*HL_WK*2
HL_DD<-HL_WK*HL_WK

p2HL<-c(HL_AA,HL_AB,HL_AC,HL_AD,HL_BB,HL_BC,HL_BD,HL_CC,HL_CD,HL_DD)
eHL<-round(p2HL*nrow(HLF1),0)
names(eHL)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
eHL

oHL_LGp<-nrow(filter(HLF1,LG_round==1,SL_round==0,SR_round==0,WK_round==0))
oHL_LGxSL<-nrow(filter(HLF1,SR_round==0,WK_round==0,LG_round!=1,SL_round!=1))
oHL_LGxSR<-nrow(filter(HLF1,LG_round!=1,SL_round==0,SR_round!=1,WK_round==0))
oHL_LGxWK<-nrow(filter(HLF1,LG_round!=1,SL_round==0,SR_round==0,WK_round!=1))
oHL_SLp<-nrow(filter(HLF1,LG_round==0,SL_round==1,SR_round==0,WK_round==0))
oHL_SLxSR<-nrow(filter(HLF1,LG_round==0,SL_round!=1,SR_round!=1,WK_round==0))
oHL_SLxWK<-nrow(filter(HLF1,LG_round==0,SL_round!=1,SR_round==0,WK_round!=1))
oHL_SRp<-nrow(filter(HLF1,LG_round==0,SL_round==0,SR_round==1,WK_round==0))
oHL_SRxWK<-nrow(filter(HLF1,LG_round==0,SL_round==0,SR_round!=1,WK_round!=1))
oHL_WKp<-nrow(filter(HLF1,LG_round==0,SL_round==0,SR_round==0,WK_round==1))

oHL<-c(oHL_LGp,oHL_LGxSL,oHL_LGxSR,oHL_LGxWK,oHL_SLp,oHL_SLxSR,oHL_SLxWK,oHL_SRp,oHL_SRxWK,oHL_WKp)
names(oHL)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
oHL

rHL<-oHL-eHL
rHL

HL_success<-c(HL_LG, (HL_LG+HL_SL)/2, (HL_LG+HL_SR)/2, (HL_LG+HL_WK)/2, HL_SL, (HL_SL+HL_SR)/2, (HL_SL+HL_WK)/2, HL_SR, (HL_SR+HL_WK)/2, HL_WK)
names(HL_success)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
HL_success

HL_success2<-c(HL_LG*HL_LG, HL_LG*HL_SL, HL_LG*HL_SR, HL_LG*HL_WK, HL_SL*HL_SL, HL_SL*HL_SR, HL_SL*HL_WK, HL_SR*HL_SR, HL_SR*HL_WK, HL_WK*HL_WK)
names(HL_success2)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
HL_success2

HL_df<-data.frame(lake="HL", hybrid=names(rHL), ex=eHL, ob=oHL, residual=rHL, success=HL_success, success2=HL_success2)

ggplot(HL_df, aes(x=success, y=residual))+
  geom_point()+
  theme_bw()

## RL F1
RLF1<- data %>% filter(lake=="RL") %>% 
  select(fish, lake, FG_round, TL_round, WB_round, WT_round, LG_round, SL_round, SR_round, WK_round) %>% 
  filter(LG_round!=0.25,SL_round!=0.25,SR_round!=0.25,WK_round!=0.25,)
nrow(RLF1)

RL_LG<-sum(RLF1$LG_round)/nrow(RLF1)
RL_SL<-sum(RLF1$SL_round)/nrow(RLF1)
RL_SR<-sum(RLF1$SR_round)/nrow(RLF1)
RL_WK<-sum(RLF1$WK_round)/nrow(RLF1)

pRL<-c(RL_LG,RL_SL,RL_SR,RL_WK)
mRL<-a*pRL

RL_AA<-mRL[1,1]
RL_AB<-mRL[1,2]+mRL[2,1]
RL_AC<-mRL[1,3]+mRL[3,1]
RL_AD<-mRL[1,4]+mRL[4,1]
RL_BB<-mRL[2,2]
RL_BC<-mRL[2,3]+mRL[3,2]
RL_BD<-mRL[2,4]+mRL[4,2]
RL_CC<-mRL[3,3]
RL_CD<-mRL[3,4]+mRL[4,3]
RL_DD<-mRL[4,4]

RL_AA<-RL_LG*RL_LG
RL_AB<-RL_LG*RL_SL*2
RL_AC<-RL_LG*RL_SR*2
RL_AD<-RL_LG*RL_WK*2
RL_BB<-RL_SL*RL_SL
RL_BC<-RL_SL*RL_SR*2
RL_BD<-RL_SL*RL_WK*2
RL_CC<-RL_SR*RL_SR
RL_CD<-RL_SR*RL_WK*2
RL_DD<-RL_WK*RL_WK

p2RL<-c(RL_AA,RL_AB,RL_AC,RL_AD,RL_BB,RL_BC,RL_BD,RL_CC,RL_CD,RL_DD)
eRL<-round(p2RL*nrow(RLF1),0)
names(eRL)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
eRL

oRL_LGp<-nrow(filter(RLF1,LG_round==1,SL_round==0,SR_round==0,WK_round==0))
oRL_LGxSL<-nrow(filter(RLF1,SR_round==0,WK_round==0,LG_round!=1,SL_round!=1))
oRL_LGxSR<-nrow(filter(RLF1,LG_round!=1,SL_round==0,SR_round!=1,WK_round==0))
oRL_LGxWK<-nrow(filter(RLF1,LG_round!=1,SL_round==0,SR_round==0,WK_round!=1))
oRL_SLp<-nrow(filter(RLF1,LG_round==0,SL_round==1,SR_round==0,WK_round==0))
oRL_SLxSR<-nrow(filter(RLF1,LG_round==0,SL_round!=1,SR_round!=1,WK_round==0))
oRL_SLxWK<-nrow(filter(RLF1,LG_round==0,SL_round!=1,SR_round==0,WK_round!=1))
oRL_SRp<-nrow(filter(RLF1,LG_round==0,SL_round==0,SR_round==1,WK_round==0))
oRL_SRxWK<-nrow(filter(RLF1,LG_round==0,SL_round==0,SR_round!=1,WK_round!=1))
oRL_WKp<-nrow(filter(RLF1,LG_round==0,SL_round==0,SR_round==0,WK_round==1))

oRL<-c(oRL_LGp,oRL_LGxSL,oRL_LGxSR,oRL_LGxWK,oRL_SLp,oRL_SLxSR,oRL_SLxWK,oRL_SRp,oRL_SRxWK,oRL_WKp)
names(oRL)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
oRL

rRL<-oRL-eRL
rRL

RL_success<-c(RL_LG, (RL_LG+RL_SL)/2, (RL_LG+RL_SR)/2, (RL_LG+RL_WK)/2, RL_SL, (RL_SL+RL_SR)/2, (RL_SL+RL_WK)/2, RL_SR, (RL_SR+RL_WK)/2, RL_WK)
names(RL_success)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
RL_success

RL_success2<-c(RL_LG*RL_LG, RL_LG*RL_SL, RL_LG*RL_SR, RL_LG*RL_WK, RL_SL*RL_SL, RL_SL*RL_SR, RL_SL*RL_WK, RL_SR*RL_SR, RL_SR*RL_WK, RL_WK*RL_WK)
names(RL_success2)<-c("LGp","LGxSL","LGxSR","LGxWK","SLp","SLxSR","SLxWK","SRp","SRxWK","WKp")
RL_success2

RL_df<-data.frame(lake="RL",hybrid=names(rRL), ex=eRL, ob=oRL, residual=rRL, success=RL_success, success2=RL_success2)

ggplot(RL_df, aes(x=success2, y=residual))+
  geom_point()+
  theme_bw()

## summary
rBen<-rCC+rLL+rLP
rLim<-rCL+rFL+rHL+rRL

rBen_pure<-rBen[1]+rBen[5]+rBen[8]+rBen[10]
rLim_pure<-rLim[1]+rLim[5]+rLim[8]+rLim[10]

## full data
full_data<-rbind(CC_df, LL_df, LP_df, CL_df, FL_df, HL_df, RL_df)

ggplot(full_data, aes(x=success, y=residual))+
  geom_point()+
  geom_smooth(method = "lm",color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("Observed-Expected Individuals")+
  xlab("Success of Source Population(s)")+
  theme(panel.grid.minor = element_blank())

ggplot(full_data, aes(x=success2, y=residual))+
  geom_point()+
  geom_smooth(method = "lm",color="grey", fill="grey", alpha=0.25)+
  theme_bw()+
  ylab("Observed-Expected Individuals")+
  xlab("Success of Source Population(s)")+
  theme(panel.grid.minor = element_blank())
  

## just pure types 
pures<-c("FGp","TLp","WBp","WTp","LGp","SLp","SRp","WKp")
pure_data<-filter(full_data, hybrid %in% pures)

ggplot(pure_data, aes(x=success, y=residual))+
  geom_point()+
  geom_smooth(method = "lm",color="grey", fill="grey", alpha=0.25)+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw()+
  ylab("Observed-Expected Pure Individuals")+
  xlab("Source Success")+
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size=14))

ggplot(pure_data, aes(x=success2, y=residual))+
  geom_point()+
  geom_smooth(method = "lm",color="grey", fill="grey", alpha=0.25)+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw()+
  ylab("Observed-Expected Pure Individuals")+
  xlab("Source Success")+
  theme(panel.grid.minor = element_blank(),
        axis.title = element_text(size=14))

ggplot(pure_data, aes(x=hybrid,y=residual))+
  geom_hline(yintercept = 0, linetype="dashed", color="grey", lwd=1)+
  geom_violin()+
  geom_point()+
  xlab("Source")+
  ylab("Residual")+
  ylim(c(-5,5))+
  scale_x_discrete(labels=c("FG","LG","SL","SR","TL","WB","WK","WT"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size=12))+
  annotate("text",x=4.5,y=4,label="Over-Represented", size=5)+
  annotate("text",x=4.5,y=-4,label="Under-Represented", size=5)

pure_data$hybrid<-factor(pure_data$hybrid, levels = c("LGp","SLp","SRp","WKp","FGp","TLp","WBp","WTp"))
ggplot(pure_data, aes(x=ex, y=ob, fill=hybrid))+
  geom_abline(slope = 1, linetype="dashed", color="grey", lwd=1)+
  geom_point(size=2, shape=21)+
  theme_bw()+
  scale_fill_manual(values = c("#eff3ff","#bdd7e7","#6baed6","#2171b5",
                               "#feedde","#fdbe85","#fd8d3c","#d94701"),
                    labels = c("LG","SL","SR","WK","FG","TL","WB","WT"),
                    name= "Source")+
  ylab("Observed Frequency")+
  xlab("Expected Frequency")+
  annotate("text",x=7.5,y=17.5,label="Over-Represented", size=5)+
  annotate("text",x=15,y=5,label="Under-Represented", size=5)+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size=10),
        legend.text = element_text(size = 10))


## just hybrids 
hybrid_data<-filter(full_data, !(hybrid %in% pures))

ggplot(hybrid_data, aes(x=success, y=residual))+
  geom_point()+
  geom_smooth(method = "lm",color="grey", fill="grey", alpha=0.25)+
  geom_hline(yintercept = 0, linetype="dashed")+
  theme_bw()+
  ylab("Observed-Expected Pure Individuals")+
  xlab("Source Success")+
  theme(panel.grid.minor = element_blank())

## testing size of lake
areas<-read.csv("./02_tidy-data/env.csv", stringsAsFactors = T) %>% 
  filter(lake.type=="Recipient") %>% 
  select(lake,area)

CC_pures<-filter(CC_df, hybrid %in% pures)
CC_lm<-lm(residual~success, data = CC_pures)
CC_slope<-data.frame(lake="CC", ecotype="benthic",slope=CC_lm$coefficients[2])

LL_pures<-filter(LL_df, hybrid %in% pures)
LL_lm<-lm(residual~success, data = LL_pures)
LL_slope<-data.frame(lake="LL", ecotype="benthic",slope=LL_lm$coefficients[2])

LP_pures<-filter(LP_df, hybrid %in% pures)
LP_lm<-lm(residual~success, data = LP_pures)
LP_slope<-data.frame(lake="LP", ecotype="benthic",slope=LP_lm$coefficients[2])

CL_pures<-filter(CL_df, hybrid %in% pures)
CL_lm<-lm(residual~success, data = CL_pures)
CL_slope<-data.frame(lake="CL", ecotype="limnetic",slope=CL_lm$coefficients[2])

FL_pures<-filter(FL_df, hybrid %in% pures)
FL_lm<-lm(residual~success, data = FL_pures)
FL_slope<-data.frame(lake="FL", ecotype="limnetic",slope=FL_lm$coefficients[2])

HL_pures<-filter(HL_df, hybrid %in% pures)
HL_lm<-lm(residual~success, data = HL_pures)
HL_slope<-data.frame(lake="HL", ecotype="limnetic",slope=HL_lm$coefficients[2])

RL_pures<-filter(RL_df, hybrid %in% pures)
RL_lm<-lm(residual~success, data = RL_pures)
RL_slope<-data.frame(lake="RL", ecotype="limnetic",slope=RL_lm$coefficients[2])

slopes<-rbind(CC_slope, LL_slope, LP_slope, CL_slope, FL_slope, HL_slope, RL_slope)
slope_df<-left_join(slopes, areas) 

ggplot(slope_df, aes(x=area, y=slope, color=ecotype))+
  geom_point()+
  theme_bw()

## redoing the slopes

CC_lm<-lm(residual~success2, data = CC_df)
CC_slope<-data.frame(lake="CC", ecotype="benthic",slope=CC_lm$coefficients[2])

LL_lm<-lm(residual~success2, data = LL_df)
LL_slope<-data.frame(lake="LL", ecotype="benthic",slope=LL_lm$coefficients[2])

LP_lm<-lm(residual~success2, data = LP_df)
LP_slope<-data.frame(lake="LP", ecotype="benthic",slope=LP_lm$coefficients[2])

CL_lm<-lm(residual~success2, data = CL_df)
CL_slope<-data.frame(lake="CL", ecotype="limnetic",slope=CL_lm$coefficients[2])

FL_lm<-lm(residual~success2, data = FL_df)
FL_slope<-data.frame(lake="FL", ecotype="limnetic",slope=FL_lm$coefficients[2])

HL_lm<-lm(residual~success2, data = HL_df)
HL_slope<-data.frame(lake="HL", ecotype="limnetic",slope=HL_lm$coefficients[2])

RL_lm<-lm(residual~success2, data = RL_df)
RL_slope<-data.frame(lake="RL", ecotype="limnetic",slope=RL_lm$coefficients[2])

slopes<-rbind(CC_slope, LL_slope, LP_slope, CL_slope, FL_slope, HL_slope, RL_slope)
slope_df<-left_join(slopes, areas) 

ggplot(slope_df, aes(x=area, y=slope, color=ecotype))+
  geom_point()+
  theme_bw()

## percent over-representation
eCC
oCC
CC_po<-((20-10)/10)*100
CC_po_df<-data.frame(lake="CC",po=CC_po)

eLL
oLL
LL_po<-((24-10)/10)*100
LL_po_df<-data.frame(lake="LL",po=LL_po)

eLP
oLP
LP_po<-((19-10)/10)*100
LP_po_df<-data.frame(lake="LP",po=LP_po)

eCL
oCL
CL_po<-((18-10)/10)*100
CL_po_df<-data.frame(lake="CL",po=CL_po)

eFL
oFL
FL_po<-((16-9)/9)*100
FL_po_df<-data.frame(lake="FL",po=FL_po)

eHL
oHL
HL_po<-((16-7)/7)*100
HL_po_df<-data.frame(lake="HL",po=HL_po)

eRL
oRL
RL_po<-((7-6)/6)*100
RL_po_df<-data.frame(lake="RL",po=RL_po)


benthics<-c("CC","LL","LP")
po_df<-rbind(CC_po_df, LL_po_df, LP_po_df, CL_po_df, FL_po_df, HL_po_df, RL_po_df) %>% 
  left_join(areas) %>% 
  mutate(ecotype = ifelse(lake%in%benthics,"Benthic Sources","Limnetic Sources"))


ggplot(po_df, aes(x=area, y=po, color=ecotype))+
  geom_point(size=4)+
  scale_color_manual(values = c("#d94701","#2171b5"))+
  geom_smooth(method = "lm", alpha=0.15, linetype="dashed", se=F)+
  theme_bw()+
  ylab("% Over-representation")+
  xlab("Lake Area")+
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title = element_text(size=14))

## testing for condition dependence
lethal<-read.csv("./02_tidy-data/lethal-2021.csv", stringsAsFactors = T) %>%
  select(fish=Fish_ID, mass=Mass_g, length=Std_length_mm, parasite=S_solidus_N) %>% 
  right_join(data) %>% 
  mutate(condition= 100000*(mass/(length^3))) %>% 
  filter(FG_round==1 | TL_round==1 | WB_round==1 | WT_round==1 | LG_round==1 | SL_round==1 | SR_round==1 | WK_round==1) %>% 
  mutate(pop=factor(case_when(
    FG_round==1 ~ "FG",
    TL_round==1 ~ "TL",
    WB_round==1 ~ "WB",
    WT_round==1 ~ "WT",
    LG_round==1 ~ "LG",
    SL_round==1 ~ "SL",
    SR_round==1 ~ "SR",
    WK_round==1 ~ "WK"))) %>% 
  mutate(infected = factor(ifelse(parasite>0,1,0)))

ggplot(lethal, aes(x=pop, y=condition, fill=pop))+
  geom_boxplot()+
  theme_bw()+
  facet_grid(~lake)

lethal_lim<-filter(lethal, pop %in% c("LG","SR"), lake!="LO")
ggplot(lethal_lim, aes(x=lake, y=condition, fill=pop))+
  geom_boxplot()+
  theme_bw()

ggplot(lethal_lim, aes(x=pop, fill=infected))+
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  facet_grid(~lake)

ggplot(lethal_lim, aes(x=lake, y=length, fill=pop))+
  geom_boxplot()+
  theme_bw()

lethal_ben<-filter(lethal, pop %in% c("WB","WT"), lake!="LO")
ggplot(lethal_ben, aes(x=lake, y=condition, fill=pop))+
  geom_boxplot()+
  theme_bw()

ggplot(lethal_ben, aes(x=lake, y=length, fill=pop))+
  geom_boxplot()+
  theme_bw()




## other condition method
lethal2<-read.csv("./02_tidy-data/lethal-2021.csv", stringsAsFactors = T) %>%
  select(fish=Fish_ID, mass=Mass_g, length=Std_length_mm, parasite=S_solidus_N) %>% 
  right_join(data) %>% 
  mutate(condition= 100000*(mass/(length^3))) %>% 
  mutate(infected = factor(ifelse(parasite>0,1,0)))

CC2<-filter(lethal2, lake=="CC") %>% 
  mutate(fitness = FG_round*CC_FG + TL_round*CC_TL + WB_round*CC_WB + WT_round*CC_WT)

LL2<-filter(lethal2, lake=="LL") %>% 
  mutate(fitness = FG_round*LL_FG + TL_round*LL_TL + WB_round*LL_WB + WT_round*LL_WT)

LP2<-filter(lethal2, lake=="LP") %>% 
  mutate(fitness = FG_round*LP_FG + TL_round*LP_TL + WB_round*LP_WB + WT_round*LP_WT)

CL2<-filter(lethal2, lake=="CL") %>% 
  mutate(fitness = LG_round*CL_LG + SL_round*CL_SL + SR_round*CL_SR + WK_round*CL_WK)

FL2<-filter(lethal2, lake=="FL") %>% 
  mutate(fitness = LG_round*FL_LG + SL_round*FL_SL + SR_round*FL_SR + WK_round*FL_WK)

HL2<-filter(lethal2, lake=="HL") %>% 
  mutate(fitness = LG_round*HL_LG + SL_round*HL_SL + SR_round*HL_SR + WK_round*HL_WK)

RL2<-filter(lethal2, lake=="RL") %>% 
  mutate(fitness = LG_round*RL_LG + SL_round*RL_SL + SR_round*RL_SR + WK_round*RL_WK)

lethal3<-rbind(CC2,LL2,LP2,CL2,FL2,HL2,RL2) %>% 
summary(lethal3)

ggplot(lethal3, aes(x=fitness, y=condition))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

ggplot(lethal3, aes(x=fitness, y=condition))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()

ggplot(lethal3, aes(x=fitness, y=length))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

## checking condition correlates
morpho<-read.csv("./02_tidy-data/linear-measurements-2021.csv") %>% 
  select(fish, body.depth) %>% 
  right_join(lethal2) %>% 
  mutate(adj.bd = body.depth/length)

ggplot(morpho, aes(x=body.depth, y=condition))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

ggplot(morpho, aes(x=length, y=condition))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

ggplot(morpho, aes(x=adj.bd, y=condition))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

## linear model for condition
condition_data<-morpho %>% select(fish, lake, mass, length, body.depth) %>% 
  na.omit()
summary(condition_data)

condition_lm<-lm(mass~length+body.depth, data=condition_data)
summary(condition_lm)

condition_data$residual<-condition_lm$residuals

condition_data2<-right_join(condition_data, select(lethal, fish, pop))

ggplot(condition_data2, aes(x=pop, y=residual, fill=pop))+
  geom_boxplot()+
  theme_bw()+
  facet_grid(~lake)

condition_data3<-right_join(condition_data, select(lethal3, fish, fitness))

ggplot(condition_data3, aes(x=fitness, y=residual))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)


## with PCA
condition_data4<-full_join(condition_data, select(data, fish, FG_round:WK_round))

condition_ben<-filter(condition_data4, lake %in% benthics) %>% 
  select(fish, lake, mass, length, body.depth, residual, FG_round:WT_round) %>% 
  na.omit() %>% 
  mutate(condition = mass/(length*body.depth))

cb_pca<-prcomp(~FG_round+TL_round+WB_round+WT_round, data=condition_ben,
               center=T, scale=T)
cb_pcs<-predict(cb_pca, newdata = condition_ben)

condition_ben$PC1<-cb_pcs[,1]
condition_ben$PC2<-cb_pcs[,2]

ggplot(condition_ben, aes(x=PC1, y=residual))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

ggplot(condition_ben, aes(x=PC1, y=condition))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

ggplot(condition_ben, aes(x=PC2, y=residual))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

ggplot(condition_ben, aes(x=PC2, y=condition))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

limnetics<-c("CL","FL","HL","RL")
condition_lim<-filter(condition_data4, lake %in% limnetics) %>% 
  select(fish, lake,mass, length, body.depth, residual, LG_round:WK_round) %>% 
  na.omit() %>% 
  mutate(condition = mass/(length*body.depth))

cl_pca<-prcomp(~LG_round+SL_round+SR_round+WK_round, data=condition_lim,
               center=T, scale=T)
cl_pcs<-predict(cl_pca, newdata = condition_lim)

condition_lim$PC1<-cl_pcs[,1]
condition_lim$PC2<-cl_pcs[,2]

ggplot(condition_lim, aes(x=PC1, y=residual))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

ggplot(condition_lim, aes(x=PC1, y=condition))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

ggplot(condition_lim, aes(x=PC2, y=residual))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

ggplot(condition_lim, aes(x=PC2, y=condition))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

## with weird fitness
morpho2 <- morpho %>% mutate(condition2 = mass/(length*body.depth)) %>% 
  full_join(lethal3)

ggplot(morpho2, aes(x=fitness, y=condition2))+
  geom_point()+
  geom_smooth(method = "lm", se=F)+
  theme_bw()+
  facet_grid(~lake)

#### 8 source scenario ####

## loading data
data<-read_csv("./02_tidy-data/final-ancestry.csv")[,-1]
data$lake<-factor(data$lake)

## LO F1
LOF1<- data %>% filter(lake=="LO") %>% 
  select(fish, lake, FG_round, TL_round, WB_round, WT_round, LG_round, SL_round, SR_round, WK_round) %>% 
  filter(FG_round!=0.25,TL_round!=0.25,WB_round!=0.25,WT_round!=0.25,LG_round!=0.25,SL_round!=0.25,SR_round!=0.25,WK_round!=0.25)
nrow(LOF1)

LO_FG<-sum(LOF1$FG_round)/nrow(LOF1)
LO_TL<-sum(LOF1$TL_round)/nrow(LOF1)
LO_WB<-sum(LOF1$WB_round)/nrow(LOF1)
LO_WT<-sum(LOF1$WT_round)/nrow(LOF1)
LO_LG<-sum(LOF1$LG_round)/nrow(LOF1)
LO_SL<-sum(LOF1$SL_round)/nrow(LOF1)
LO_SR<-sum(LOF1$SR_round)/nrow(LOF1)
LO_WK<-sum(LOF1$WK_round)/nrow(LOF1)

## simplify by ecotype
LO_ben<-LO_FG+LO_TL+LO_WB+LO_WT
LO_lim<-LO_LG+LO_SL+LO_SR+LO_WK

LO_BB<-LO_ben*LO_ben
LO_BL<-2*LO_ben*LO_lim
LO_LL<-LO_lim*LO_lim

eLO_BB<-round(LO_BB*nrow(LOF1),0)
eLO_BL<-round(LO_BL*nrow(LOF1),0)
eLO_LL<-round(LO_LL*nrow(LOF1),0)
eLO<-c(eLO_BB, eLO_BL, eLO_LL)
names(eLO)<-c("B*B","B*L","L*L")

oLO_BB<-nrow(filter(LOF1,LG_round==0,SL_round==0,SR_round==0,WK_round==0))
oLO_LL<-nrow(filter(LOF1,FG_round==0,TL_round==0,WB_round==0,WT_round==0))
oLO_BL<-nrow(LOF1)-(oLO_BB+oLO_LL)
oLO<-c(oLO_BB, oLO_BL, oLO_LL)
names(oLO)<-c("B*B","B*L","L*L")

rLO<-oLO-eLO
rLO

## by region
LO_K<-LO_TL+LO_WT+LO_SL+LO_WK
LO_M<-LO_LG+LO_SR+LO_FG+LO_WB

LO_KK<-LO_K*LO_K
LO_KM<-2*LO_K*LO_M
LO_MM<-LO_M*LO_M

eLO_KK<-round(LO_KK*nrow(LOF1),0)
eLO_KM<-round(LO_KM*nrow(LOF1),0)
eLO_MM<-round(LO_MM*nrow(LOF1),0)
eLO<-c(eLO_KK, eLO_KM, eLO_MM)
names(eLO)<-c("K*K","K*M","M*M")

oLO_KK<-nrow(filter(LOF1,LG_round==0,SR_round==0,FG_round==0,WB_round==0))
oLO_MM<-nrow(filter(LOF1,TL_round==0,WT_round==0,SL_round==0,WK_round==0))
oLO_KM<-nrow(LOF1)-(oLO_KK+oLO_MM)
oLO<-c(oLO_KK, oLO_KM, oLO_MM)
names(oLO)<-c("K*K","K*M","M*M")

rLO<-oLO-eLO
rLO

## interaction?
LO_BK<-LO_TL+LO_WT
LO_BM<-LO_FG+LO_WB
LO_LK<-LO_SL+LO_WK
LO_LM<-LO_LG+LO_SR

LO_AA<-LO_BK*LO_BK
LO_AB<-LO_BK*LO_BM*2
LO_AC<-LO_BK*LO_LK*2
LO_AD<-LO_BK*LO_LM*2
LO_BB<-LO_BM*LO_BM
LO_BC<-LO_BM*LO_LK*2
LO_BD<-LO_BM*LO_LM*2
LO_CC<-LO_LK*LO_LK
LO_CD<-LO_LK*LO_LM*2
LO_DD<-LO_LM*LO_LM

p2LO<-c(LO_AA,LO_BB,LO_CC,LO_DD)
eLO<-round(p2LO*nrow(LOF1),0)
names(eLO)<-c("BKxBK","BMxBM","LKxLK","LMxLM")
eLO

oLO_BKxBK<-nrow(filter(LOF1,LG_round==0,SL_round==0,SR_round==0,WK_round==0,FG_round==0,WB_round==0))
oLO_BMxBM<-nrow(filter(LOF1,LG_round==0,SL_round==0,SR_round==0,WK_round==0,TL_round==0,WT_round==0))
oLO_LKxLK<-nrow(filter(LOF1,LG_round==0,FG_round==0,SR_round==0,WB_round==0,TL_round==0,WT_round==0))
oLO_LMxLM<-nrow(filter(LOF1,SL_round==0,FG_round==0,WK_round==0,WB_round==0,TL_round==0,WT_round==0))
oLO<-c(oLO_BKxBK,oLO_BMxBM,oLO_LKxLK,oLO_LMxLM)

rLO<-oLO-eLO
rLO



## 
q<-read.csv("./02_tidy-data/assortment.csv", stringsAsFactors = T)

ggplot(q, aes(x=fst, y=residual))+
  geom_smooth(method = "lm", color="grey75", fill="grey90")+
  geom_jitter()+
  xlab("Fst")+
  ylab("Over/Under Representation")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=14))

ggplot(q, aes(x=morph, y=residual))+
  geom_smooth(method = "lm", color="grey75", fill="grey90")+
  geom_jitter()+
  xlab("Morph Divergence")+
  ylab("Over/Under Representation")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=14))


