#Code by Matthew Kustra. For questions: mkustra@ucsc.edu
# R script that was used for making initial graphs/analysis for 
#"Microbes as manipulators of developmental life-history" 
#Matthew Kustra and Tyler Carrier

#Graphs specifically generated:
###Fig. 1B,C
### Fig S1B
###Fig. S2, S3
# 1. Setting up ------------------------------------------------
# * 1.a Loading up required libraries -------------------------------------
library(tidyverse)
library(viridis)
library(patchwork)
library(npreg)
# Theme -------------------------------------------------------------------
#* 1.d Setting Different plotting theme ----------------------
#Heatmaps
mytheme <-
  theme_classic() + theme(
    legend.position = "bottom",
    # this puts legend on the bottom
    # this makes the axis titles in bold,
    axis.line = element_line(
      color = "black", size =
        0
    ),
    # Makes the axis line black and  thicker
    text = element_text(size = 12, color = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    ),
    strip.background = element_rect(
      color = "black", fill = "white", size =
        2
    ),
    strip.text = element_text(color = "black")
  )
#function graphs
mytheme2 <-
  theme_classic() + theme(
    legend.position = "top",
    # this puts legend on the bottom
    axis.title = (element_text(face = "bold")),
    # this makes the axis titles in bold,
    axis.line = element_line(
      color = "black", size =
        1
    ),axis.text = element_text(color="black"),
    # Makes the axis line black and  thicker
    text = element_text(
      size = 18, face = "bold", color =
        "black"
    )
  )

# 2.Male Killing ------------------------------------------------------

# * 2.a Loading up/processing data -----------------------------------------------------
Data_all<-read_csv("MK_gridV.csv")%>%
  rename(MK=mk)


#Calculate stable egg size if they are near
#This accounts for situations where stability can be approximated even if not reached based on direction of evolution left or right of where the calculation was performed.
Data_all2<-Data_all%>%
  group_by(cost,MK,b)%>%
  summarise(minEgg=min(Egg_Size2[DIRF=="Smaller"]),maxEgg=max(Egg_Size2[DIRF=="Larger"]),RESWF=mean(RESWF))%>%
  ungroup()%>%
  mutate(Egg_Size2=(minEgg+maxEgg)/2,DIRF="Stable")%>%
  select(-minEgg,-maxEgg)%>%
  mutate(Stat="Calculated")


#Add back in stable egg size that were found 
Data_all3<-Data_all%>%
  select(cost,MK,b,Egg_Size2,DIRF,RESWF)%>%
  filter(DIRF=="Stable")%>%
  mutate(Stat="OG")%>%
  rbind(Data_all2)

#* 2.b  Figure S2----------------------------------------------------------------------
(b300<-Data_all2%>%
   filter(DIRF=="Stable")%>%
   filter(b%in%c(300))%>%
   ggplot(aes(x=MK,y=-cost,fill=Egg_Size2))+
   geom_tile()+
   facet_wrap(~b,ncol=5)+
   scale_fill_viridis(
     option = "magma",
     na.value = "grey",
     direction = -1
   )+mytheme+
   labs(x="Male killing rate",y="Benefit")+ scale_x_continuous(expand = c(0.0, 0.0)) +
   scale_y_continuous(
     expand = c(0.0, 0.0)))

(b900<-Data_all3%>%
    filter(DIRF=="Stable")%>%
    filter(b%in%c(900))%>%
    ggplot(aes(x=MK,y=-cost,fill=Egg_Size2))+
    geom_tile()+
    facet_wrap(~b,ncol=5)+
    scale_fill_viridis(
      option = "magma",
      na.value = "grey",
      direction = -1
    )+mytheme+
    labs(x="Male killing rate",y="Benefit")+
    scale_size_manual(values=c(1,0))+ scale_x_continuous(expand = c(0.0, 0.0)) +
    scale_y_continuous(
      expand = c(0.0, 0.0)))

(b1500<-Data_all3%>%
    filter(DIRF=="Stable")%>%
    filter(b%in%c(1500))%>%
    ggplot(aes(x=MK,y=-cost,fill=Egg_Size2))+
    geom_tile()+
    facet_wrap(~b,ncol=5)+
    scale_x_continuous(expand = c(0.0, 0.0)) +
    scale_y_continuous(
      expand = c(0.0, 0.0))+

    scale_fill_viridis(
      option = "magma",
      na.value = "grey",
      direction = -1
    )+mytheme+
    labs(x="Male killing rate",y="Benefit"))

(b3000<-Data_all3%>%
    filter(DIRF=="Stable")%>%
    filter(b%in%c(3000))%>%
    ggplot(aes(x=MK,y=-cost,fill=Egg_Size2))+
    geom_tile()+
    facet_wrap(~b,ncol=5)+
    scale_x_continuous(expand = c(0.0, 0.0)) +
    scale_y_continuous(
      expand = c(0.0, 0.0))+
    scale_fill_viridis(
      option = "magma",
      na.value = "grey",
      direction = -1
    )+mytheme+
    labs(x="Male killing rate",y="Benefit"))


#Rough version of Figure S2. Manually cleaned up
(Fig.S2<-((b300+b900)/(b1500+b3000))&labs(fill="Egg Size um")&theme(legend.position = "right",aspect.ratio = 1))

ggsave("FigS2_V.SVG",height = (2*183),width=(2*183),units="mm")


# * 2.c Male killing part of Fig. 1B  ---------------------------------------------------------------
Data_all%>%
  filter(Egg_Size2<=500,(MK==0.84&cost==-1.5)|(MK==0&cost==0)|(MK==0.84&cost==0),DIRF!="Either",b%in%c(1000))%>%
  mutate(cat=ifelse((MK==0.84&cost==-1.5),"EG+MK",ifelse((MK==0&cost==0),"No Inf","MK Only")))%>%
  mutate(cat=factor(cat,levels=rev(c("No Inf","MK Only","EG+MK"))))%>%
  ggplot(aes(x=Egg_Size2,y=cat))+
  geom_tile(aes(fill=DIRF))+
  facet_wrap(~b,ncol=5)+
  mytheme+
  scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_discrete(
    expand = c(0.0, 0.0))+
  scale_fill_manual(values=c("#364678","#8AAFD0","white"),name="Evolutionary trajectory")+
  xlab("Egg diameter (um)")+
  ylab("Infection type")+theme(legend.background = element_rect(fill="grey"),legend.position = "right")

ggsave("Fig1B_MK_V.svg",height = 183/1.5,width=1.5*183,units="mm")

# * 2.d Calculations for relative host survival -----------------------------------------------------------
#Filter out for an example enhanced growth and male killing
DataTy<-Data_all3%>%
    filter((MK==0.84&cost==-1.5)|(MK==0&cost==0)|(MK==0.84&cost==0)|(MK==0&cost==-1.5))%>%
    mutate(cat=ifelse((MK==0.84&cost==-1.5),"EG+MK",ifelse((MK==0&cost==0),"No Inf",ifelse((MK==0.84&cost==0),"MK Only","EG Only"))))
  

#just select get EG+Mk and infected
Dataty2<-DataTy%>%filter(cat!="MK Only",cat!="EG Only") %>%
  drop_na()%>%
  filter(Egg_Size2!=-Inf)

#spline fit uninfected data egg size predicting B parameter
UMod<-ss(Dataty2$Egg_Size2[Dataty2$cat=="No Inf"],Dataty2$b[Dataty2$cat=="No Inf"],nknots=10)

#spline fit infected data
InfMod<-ss(Dataty2$Egg_Size2[Dataty2$cat=="EG+MK"],Dataty2$b[Dataty2$cat=="EG+MK"],nknots=10)

#calculate predictions from spline fit
INF<-predict(InfMod,seq(50,725))%>%
  rename(PredInf=y,InfSE=se,Egg=x)

#calculate predictions from spline fit
U<-predict(UMod,seq(50,725))%>%
  rename(PredU=y,uSE=se)%>%
  select(-x)
#bind predictions then subtract differencs
All<-cbind(INF,U)%>%
  mutate(Diff=PredU-PredInf)

#calculate difference in survival of infected relative to uninfected
All$SurvDif<-exp(-All$PredInf/All$Egg)/exp(-All$PredU/All$Egg)


# * 2.e Male Killing Figure S1B ------------------------------------------------------
#make plot of stable egg size for given b parameter
(FigS1B<-ggplot(Dataty2,aes(x=b,y=Egg_Size2,color=cat,shape=cat))+
  geom_line(size=1)+
  mytheme2+
  xlab("B parameter")+
  ylab("Stable Egg size (um)")+
  scale_color_manual(values = c("#33487D","#7D8F99"),name="")+
  scale_shape_manual(name="",values=c(16,17))+theme(legend.position = c(0.2,0.88))+scale_y_continuous(limits=c(0,750),breaks=c(0,250,500,750)))


  

#survival difference at a given egg size based on b parameter needed to achieve that.

# * 2.f Figure 1C Male killing ---------------------------------------------------------------
(Fig1cMK<-ggplot(All,aes(x=Egg,y=SurvDif))+
  geom_line(size=1,color="#33487D") +mytheme2+
  ylab("Relative host survival")+
  xlab("Egg size (um)")+
  scale_x_continuous(limits=c(0,750),breaks=c(0,250,500,750))+
   geom_vline(aes(xintercept=100)))
  #scale_y_continuous(breaks=c(0,5,10,15),limits=c(0,15)))
#write_csv(All,"Figure1CMaleKillingVolume.csv")

#ggsave("Fig1CMK.svg",height = 183/1.5,width=1.5*183,units="mm")
#Manually integrated with Feminization to make Fig.1C later on.


# Feminization grid results -----------------------------------------------
# 3. Feminization ------------------------------------------------------
# * 3.a Loading up/processing data ----------------------------------------


Data_allF<-read_csv("Fem_gridV.csv")

#Calculate stable egg size if they are near
Data_all2F<-Data_allF%>%
  group_by(cost,Feminization,b)%>%
  summarise(minEgg=min(Egg_Size2[DIRF=="Smaller"],na.rm=T),maxEgg=max(Egg_Size2[DIRF=="Larger"],na.rm=T),RESWF=mean(RESWF))%>%
  ungroup()%>%
  mutate(Egg_Size2=(minEgg+maxEgg)/2,DIRF="Stable")%>%
  select(-minEgg,-maxEgg)%>%
  mutate(Stat="Calculated")


#Add back in stable egg size that were found 
Data_all3F<-Data_allF%>%
  select(cost,Feminization,b,Egg_Size2,DIRF,RESWF)%>%
  filter(DIRF=="Stable")%>%
  mutate(Stat="OG")%>%
  rbind(Data_all2F)

# * 3.b Fig S3----------------------------------------------------------------------
(b300F<-Data_all2F%>%
   filter(DIRF=="Stable")%>%
   filter(b%in%c(300))%>%
   ggplot(aes(x=Feminization,y=-cost,fill=Egg_Size2))+
   geom_tile()+
   facet_wrap(~b,ncol=5)+
   scale_fill_viridis(
     option = "magma",
     na.value = "grey",
     direction = -1
   )+mytheme+
   labs(x="Feminization rate",y="Benefit")+ scale_x_continuous(expand = c(0.0, 0.0)) +
   scale_y_continuous(
     expand = c(0.0, 0.0)))

(b900F<-Data_all3F%>%
    filter(DIRF=="Stable")%>%
    filter(b%in%c(900))%>%
    ggplot(aes(x=Feminization,y=-cost,fill=Egg_Size2))+
    geom_tile()+
    facet_wrap(~b,ncol=5)+
    scale_fill_viridis(
      option = "magma",
      na.value = "grey",
      direction = -1
    )+mytheme+
    labs(x="Feminization rate",y="Benefit")+
    scale_size_manual(values=c(1,0))+ scale_x_continuous(expand = c(0.0, 0.0)) +
    scale_y_continuous(
      expand = c(0.0, 0.0)))

(b1500F<-Data_all3F%>%
    filter(DIRF=="Stable")%>%
    filter(b%in%c(1500))%>%
    ggplot(aes(x=Feminization,y=-cost,fill=Egg_Size2))+
    geom_tile()+
    facet_wrap(~b,ncol=5)+
    scale_fill_viridis(
      option = "magma",
      na.value = "grey",
      direction = -1
    )+mytheme+
    labs(x="Feminization rate",y="Benefit")+ scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(
    expand = c(0.0, 0.0)))

(b3000F<-Data_all3F%>%
    filter(DIRF=="Stable")%>%
    filter(b%in%c(3000))%>%
    ggplot(aes(x=Feminization,y=-cost,fill=Egg_Size2))+
    geom_tile()+
    facet_wrap(~b,ncol=5)+
    scale_fill_viridis(
      option = "magma",
      na.value = "grey",
      direction = -1
    )+mytheme+
    labs(x="Feminization rate",y="Benefit")+ scale_x_continuous(expand = c(0.0, 0.0)) +
    scale_y_continuous(
      expand = c(0.0, 0.0)))


(FigS3F<-(((b300F+b900F)/(b1500F+b3000F))&labs(fill="Egg Size um")&theme(legend.position = "right",aspect.ratio = 1)))


ggsave("Fig_S3F_v.svg",height = 183,width=183,units="mm")

# * 3.c Feminization Figure 1B ---------------------------------------------------------------
Data_allF%>%
  filter(Egg_Size2<=500,(Feminization==0.845&cost==-1.5)|(Feminization==0.5&cost==0)|(Feminization==0.845 &cost==0),DIRF!="Either",b%in%c(1000))%>%
  mutate(cat=ifelse((Feminization==0.845&cost==-1.5),"EG+F",ifelse((Feminization==0.5&cost==0),"No Inf","F Only")))%>%
  mutate(cat=factor(cat,levels=rev(c("No Inf","EG Only","F Only","EG+F"))))%>%
  ggplot(aes(x=Egg_Size2,y=cat))+
  geom_tile(aes(fill=DIRF))+
  facet_wrap(~b,ncol=5)+
  mytheme+
  scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_discrete(
    expand = c(0.0, 0.0))+
  scale_fill_manual(values=c("#364678","#8AAFD0","white"),name="Evolutionary trajectory")+
  xlab("Egg diameter (um)")+
  ylab("Infection type")+theme(legend.background = element_rect(fill="grey"),legend.position = "right")

#ggsave("FemFig1b_V.svg",height = 183/1.5,width=1.5*183,units="mm")

# * 3.d Calculations for relative host fitness -----------------------------------------------------------
#Filter out for an example enhanced growth and Feminization
DataTyF<-Data_all3F%>%
  filter((Feminization==0.845&cost==-1.5)|(Feminization==0.5&cost==0)|(Feminization==0.845&cost==0)|(Feminization==0.5&cost==-1.5),DIRF!="Either")%>%
  mutate(cat=ifelse((Feminization==0.845&cost==-1.5),"EG+F",ifelse((Feminization==0.5&cost==0),"No Inf",ifelse((Feminization==0.845&cost==0),"F Only","EG Only"))))%>%
  mutate(cat=factor(cat,levels=rev(c("No Inf","EG Only","F Only","EG+F"))))

#pm;u get EG+Mk and infected
Dataty2F<-DataTyF%>%filter(cat!="F Only",cat!="EG Only") %>%
  drop_na()%>%
  filter(Egg_Size2!=-Inf)
#spline fit uninfected data egg size predicting parameter
UModF<-ss(Dataty2F$Egg_Size2[Dataty2F$cat=="No Inf"],Dataty2F$b[Dataty2F$cat=="No Inf"],nknots=10)

#spline fit infected data
InfModF<-ss(Dataty2F$Egg_Size2[Dataty2F$cat=="EG+F"],Dataty2F$b[Dataty2F$cat=="EG+F"],nknots=10)

#calculate predictions from spline fit
INFF<-predict(InfModF,seq(50,725))%>%
  rename(PredInf=y,InfSE=se,Egg=x)

#calculate predictions from spline fit
UF<-predict(UModF,seq(50,725))%>%
  rename(PredU=y,uSE=se)%>%
  select(-x)
#bind predictions then subtract differencs
AllF<-cbind(INFF,UF)%>%
  mutate(Diff=PredU-PredInf)

#calculate difference in survival of infected relative to uninfected
AllF$SurvDif<-exp(-AllF$PredInf/AllF$Egg)/exp(-AllF$PredU/AllF$Egg)

# * 3.e Feminization Figure S1B -------------------------------------------
#make plot of stable egg size for given b parameter
(FigS1BFem<-ggplot(Dataty2F,aes(x=b,y=Egg_Size2,color=cat,shape=cat))+
    geom_line(size=1)+
    mytheme2+
    xlab("B parameter")+
    ylab("Stable Egg size (um)")+
    scale_color_manual(values = c("#7FB2D5","#7F909A"),name="")+
    scale_shape_manual(name="",values=c(16,17))+theme(legend.position = c(0.2,0.88))+scale_y_continuous(limits=c(0,750),breaks=c(0,250,500,750)))

#ggsave("FemFigS1b.svg",height = 183/1.5,width=1.5*183,units="mm")


# * 3.f Figure 1C Feminization --------------------------------------------
#survival difference at a given egg size based on b parameter needed to achieve that.
(Fig1cFem<-ggplot(AllF,aes(x=Egg,y=SurvDif))+
    geom_line(size=1,color="#7FB2D5") +mytheme2+
    ylab("Survival Infected/Uninfected")+
    xlab("Egg size (um)")+
    scale_x_continuous(limits=c(0,750),breaks=c(0,250,500,750)))


#ggsave("Fig1CFem.svg",height = 183/1.5,width=1.5*183,units="mm")


Fig1cmk<-All%>%
  select(Egg,SurvDif)%>%
  mutate(Cat="MK+EG")
Fig1cF<-AllF%>%
  select(Egg,SurvDif)%>%
  mutate(Cat="FEM+EG")

Fig1cD<-rbind(Fig1cmk,Fig1cF)

write_csv(Fig1cD,"Figure1C-V.csv")
(Fig1c<-ggplot(Fig1cD,aes(x=Egg,y=SurvDif,color=Cat,group=Cat))+
    geom_line(size=1) +mytheme2+
    ylab("Survival Infected/Uninfected")+
    xlab("Egg size (um)")+
    scale_x_continuous(limits=c(0,750),breaks=c(0,250,500,750)))

