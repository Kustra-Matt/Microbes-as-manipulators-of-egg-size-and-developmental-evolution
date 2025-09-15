#Code by xxx. For questions: xxxxx
# R script that was used for making initial graphs/analysis for 
#"Microbes as manipulators of developmental life-history" 
#xxxx

#Graphs specifically generated:
### Fig.2,3
### Fig S6,7,8,9,10,11

# 1. Setting up ------------------------------------------------
# * 1.a Loading up required libraries -------------------------------------
library(tidyverse)
library(ggprism)
library(data.table)
library(viridis)
library(patchwork)

# Theme -------------------------------------------------------------------
#* 1.b Setting Different plotting theme ----------------------
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
    text = element_text(size = 20, color = "black"),
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


# 2. Figure 2 ------------------------------------------------------
#Load up data
Fig2Data<-read_csv("Data/Fig2QG.csv")

#Make Figure 2 
(Fig2<-ggplot(Fig2Data,aes(x=Generation,y=EggSize,color=Cat2,linetype=GRCat))+
  geom_line(size=1)+
  scale_linetype_manual(values=c("dashed","solid"))+
  scale_color_manual(values=c("#83B4D7","#3F5384","#C9CACC"))+
  theme_prism()+
  scale_x_continuous(guide=guide_prism_minor(),
                     limits = c(0, 1200),
                     breaks = c(0,300, 600,900,1200)) +
  scale_y_continuous(guide=guide_prism_minor(),
                     limits = c(100, 400),
                     breaks = c(100,200, 300, 400)))

#ggsave("Figure2.svg",height = 183,width=(2*183),units="mm")

# 3. Figure S5 ----------------------------------------------------------------
#Load up data
FigS6Data<-read_csv("Data/FigS6LossQG.csv")

#Sex ratio for part A of Figure S5
(S6ASR<-Fig2Data%>%
  filter(Cat%in%c("MK + EG","Fem + EG"))%>%
  ggplot(aes(x=Generation,y=SexRatio,color=Cat2))+
  geom_line(size=1,linetype="dashed")+
  scale_color_manual(values=c("#83B4D7","#3F5384"))+
  theme_prism()+
  scale_x_continuous(guide=guide_prism_minor(),
                     limits = c(0, 4000),
                     breaks = c(0,1000, 2000,3000,4000)) +
  scale_y_continuous(guide=guide_prism_minor(),
                     limits = c(0, 0.5),
                     breaks = c(0,0.25, 0.5)))
  
#ggsave("FigureS6_SR.svg",height = 183,width=(2*183),units="mm")
#Loss of microbe time series
#first need to shift everything by 1200
(S6BLoss<-FigS6Data%>%
    mutate(Generation=Generation+1200)%>%
    ggplot(aes(x=Generation,y=EggSize,color=Cat))+
    geom_line(size=1)+
    scale_color_manual(values=c("#83B4D7","#3F5384"))+
    theme_prism()+
    scale_x_continuous(guide=guide_prism_minor(),
                       limits = c(1200, 2400),
                       breaks = c(1200,1800, 2400)) +
    scale_y_continuous(guide=guide_prism_minor(),
                       limits = c(100, 400),
                       breaks = c(100,200,300,400)))

#ggsave("FigureS6_Loss.svg",height = 183,width=(2*183),units="mm")

#Figure S6 was assembled in illustrator based on these pieces and Figure 2 pieces.

# 4. Figure 3 Sex ratio,  V, B -----------------------------------------------
Data_all<-fread("Data/MK_QGPop.csv")
#Filtering for transition only and enhanced growth only
#then finding the minimum male killing rate needed for transition
D2<-Data_all%>%
  filter(EggSize>300,GR==0.5)%>%
  group_by(B,V_b) %>% slice(which.min(MK))

#Just some summary statistics for values in this figure.
D2%>%
  mutate(`Proportion Female`=1-SexRatio)%>%
  filter(Generation==3999)%>%
  summary()

#Figure 3A
(a<-D2%>%
    mutate(`Proportion Female`=1-SexRatio)%>%
    filter(Generation==3999)%>%
    ggplot(aes(x=B,y=V_b))+
    geom_tile(aes(fill=`Proportion Female`))+
    scale_fill_distiller(
      palette = "Blues",
      na.value = "grey",
      direction = 1,
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",frame.linewidth=0.5
      ))+
    mytheme+
    labs(x="B parameter",y="V parameter")+ scale_x_continuous(expand = c(0.0, 0.0)) +
    scale_y_continuous(expand = c(0.0, 0.0))+
    theme(legend.position = "right",
          axis.text = element_text(color="black"))+
    theme(panel.background = element_rect(fill="grey50")))

#Figure 3B
(b<-Data_all%>%
    filter(Generation==3999,GR==0.5)%>%
    group_by(B,V_b)%>%
    summarise(Min=min(EggSize),Max=max(EggSize),Mean=mean(EggSize),Median=median(EggSize))%>%
    pivot_longer(cols = -c(B,V_b),values_to = "Egg size",names_to = "Metric")%>%
    filter(Metric!="Median")%>%
    mutate(Metric=factor(Metric,levels=c("Min","Mean","Max")))%>%
    filter(Metric=="Mean")%>%
    ggplot(aes(x=B,y=V_b))+
    geom_tile(aes(fill=`Egg size`))+
    scale_fill_viridis(option = "mako",na.value = "grey",direction = -1,
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",frame.linewidth=0.5)
    )+
    mytheme+
    labs(x="B parameter",y="V parameter")+ scale_x_continuous(expand = c(0.0, 0.0)) +
    scale_y_continuous(expand = c(0.0, 0.0))+
    theme(legend.position = "right",
          axis.text = element_text(color="black"))+
    theme(panel.background = element_rect(fill="grey50")))

#Putting it all together
(fig3<-((a+b)&theme(legend.position = "right"))+ plot_layout(guides = "collect")+
    plot_annotation(tag_levels = "A",tag_suffix = ")"))
#Legend fixed in illustrator afterwards
#ggsave("Figure3_MK.svg",height = 183,width=(2*183),units="mm")



# 5. Figure S7 ------------------------------------------------------------
#Load up Data
Data_S7<-fread("Data/FigS7GQ.csv")

#Set up labels for this graph
blabs<-c("B = 300","B = 600")
names(blabs)=c(300,600)
vlabs<-c("V = 250","V = 750")
names(vlabs)=c(250,750)

#Make the figure
(FigS7<-ggplot(Data_S7%>%filter(GR>=0,Generation>3000,V_b%in%c(250,750),B%in%c(300,600)),aes(x=MK,y=GR))+
  geom_tile(aes(fill=EggSize))+
  facet_grid(V_b~B,labeller = labeller(B=blabs,V_b=vlabs))+
  scale_fill_viridis(option = "mako",na.value = "grey",direction = -1,
                     guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",frame.linewidth=0.5)
  )+
  mytheme+
  labs(x="Male killing rate",y="Benefit",fill="Egg size (\U003BCM)")+ scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(expand = c(0.0, 0.0))+theme(legend.position = "right"))

#ggsave("FigS7_V.svg",height = 183,width=1.3*183,units="mm")


# 6. Figure S8 Feminization version of Figure 3 -----------------------------------------------
#Read in data
FData_all<-fread("Data/Fem_QGPop.csv")
#Filter out for transition and miniumum feminization needed for manipultion
#MK is actually feminization
FD2<-FData_all%>%
  filter(EggSize>300,GR==0.5)%>%
  group_by(B,V_b) %>% slice(which.min(MK))

#Fig S7A
Fa<-FD2%>%
  mutate(`Proportion Female`=1-SexRatio)%>%
  filter(Generation==3999)%>%
  ggplot(aes(x=B,y=V_b))+
  geom_tile(aes(fill=`Proportion Female`))+
  scale_fill_distiller(
    palette = "Blues",
    na.value = "grey",
    direction = 1,
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",frame.linewidth=0.5
    ))+
  mytheme+
  labs(x="B parameter",y="V parameter")+ scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(expand = c(0.0, 0.0))+
  theme(legend.position = "right",
        axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(fill="grey50"))

#Fig S7B
Fb<-FData_all%>%
  filter(Generation==3999,GR==0.5)%>%
  group_by(B,V_b)%>%
  summarise(Min=min(EggSize),Max=max(EggSize),Mean=mean(EggSize),Median=median(EggSize))%>%
  pivot_longer(cols = -c(B,V_b),values_to = "Egg size",names_to = "Metric")%>%
  filter(Metric!="Median")%>%
  mutate(Metric=factor(Metric,levels=c("Min","Mean","Max")))%>%
  filter(Metric=="Mean")%>%
  ggplot(aes(x=B,y=V_b))+
  geom_tile(aes(fill=`Egg size`))+
  scale_fill_viridis(option = "mako",na.value = "grey",direction = -1,
                     guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",frame.linewidth=0.5)
  )+
  mytheme+
  labs(x="B parameter",y="V parameter")+ scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(expand = c(0.0, 0.0))+
  theme(legend.position = "right",
        axis.text = element_text(color="black"))+
  theme(panel.background = element_rect(fill="grey50"))

#Putting it all together
(FfigS8<-((Fa+Fb)&theme(legend.position = "right"))+ plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "A",tag_suffix = ")"))

#ggsave("FigureS8_Fem.svg",height = 183,width=(2*183),units="mm")


# 7. Percentage transition and summary statistics -------------------------------------------------------
#Of all combinations, what % stays plank, what transitions, and what stays lecithotroph

#  * 7.a Male killing -----------------------------------------------------
#categorizing things as lecitotrophy if greater than 300 um.
All<-Data_all%>%filter(Generation>3000)%>%
  filter(!(GR==0&MK==0))%>%
  #filter(GR==0)%>%
  mutate(Transition=ifelse(EggSize>=300,"Lecithotrophy","Planktotrophy"))

#Calculating what default is without infection at different B parameters
NoInf<-Data_all%>%filter(GR==0,Generation>3000,MK==0)%>%
  mutate(NoInf=ifelse(EggSize>=300,"Lecithotrophy","Planktotrophy"))%>%
  group_by(B)%>%
  summarise(NoInf=max(NoInf),EggSizeUn=mean(EggSize))

#Joining everything together
All_inf<-left_join(All,NoInf,by=c("B"))

#Calculate percentages of transitions for combinations of GR and Male killing
All_inf%>%
  mutate(LHT=ifelse(Transition==NoInf,NoInf,"P to L"))%>%
  mutate(MKstat=ifelse(MK>0,"MK","NoMK"))%>%
  mutate(Gstat=ifelse(GR>0,"GR","NoGR"))%>%
  group_by(Gstat,MKstat,LHT)%>%
  summarize(N=length(LHT))%>%
  group_by(MKstat,Gstat)%>%
  mutate(Total=sum(N))%>%
  mutate(Percent=(N/Total)*100)%>%
 # filter(LHT=="P to L")%>%
  select(-c(N,Total))

#Summary statistics
All_inf%>%
  mutate(LHT=ifelse(Transition==NoInf,NoInf,"P to L"))%>%
  filter(LHT=="P to L")%>%
  summarise(Min=min(EggSize),Max=max(EggSize),Mean=mean(EggSize),Median=median(EggSize),SD=sd(EggSize))

# * 7.b Feminization ------------------------------------------------------
FAll<-FData_all%>%filter(Generation>3000)%>%
  filter(!(GR==0&MK==0.5))%>%
  mutate(Transition=ifelse(EggSize>=300,"Lecithotrophy","Planktotrophy"))

FNoInf<-FData_all%>%filter(GR<0.5,MK<0.501)%>%
  mutate(NoInf=ifelse(EggSize>=300,"Lecithotrophy","Planktotrophy"))%>%
  group_by(B)%>%
  summarise(NoInf=max(NoInf),EggSizeUn=mean(EggSize))

#Merge dataframes
FAll_inf<-left_join(FAll,FNoInf,by=c("B"))

#Calculate percentages of transitions for combinations of GR and Feminization
FAll_inf%>%
  mutate(LHT=ifelse(Transition==NoInf,NoInf,"P to L"))%>%
  mutate(MKstat=ifelse(MK>0.5,"MK","NoMK"))%>%
  mutate(Gstat=ifelse(GR>0,"GR","NoGR"))%>%
  group_by(Gstat,MKstat,LHT)%>%
  summarize(N=length(LHT))%>%
  group_by(MKstat,Gstat)%>%
  mutate(Total=sum(N))%>%
  mutate(Percent=(N/Total)*100)%>%
  #filter(LHT=="P to L")%>%
  select(-c(N,Total))


#Summary statistics on egg size
FAll_inf%>%
  mutate(LHT=ifelse(Transition==NoInf,NoInf,"P to L"))%>%
  filter(LHT=="P to L")%>%
  summarise(Min=min(EggSize),Max=max(EggSize),Mean=mean(EggSize),Median=median(EggSize),SD=sd(EggSize))

# 8. Fig S9 -----------------------------------------------
(Fem_E<-FAll_inf%>%
   filter(V_b%in%c(300,400,500),B%in%c(250),GR==0.5)%>%
   ggplot(aes(x=MK,y=EggSize,color=factor(V_b),group=factor(V_b)))+
   geom_line(linewidth=1.25)+
   theme_prism()+
   geom_hline(aes(yintercept=138,color="No microbe"),linewidth=1.25)+
   xlab("Feminization rate")+
   ylab("Egg size (\U003BCM)")+
   scale_color_manual(values=c("#562367","#515797","#338A9F","#7D8F99"))+
   scale_y_continuous(guide=guide_prism_minor(),limits=c(100,650),breaks=c(100,375,650),expand = c(0, 0))+
   scale_x_continuous(guide=guide_prism_minor(),limits=c(0.5,1.0),breaks=c(0.5,0.75,1.0),expand = c(0.0, 0.0)))


(MK_E<-All_inf%>%
    filter(V_b%in%c(300,400,500),B%in%c(600),GR==0.5)%>%
    ggplot(aes(x=MK,y=EggSize,color=factor(V_b),group=factor(V_b)))+
    geom_line(linewidth=1.25)+
    theme_prism()+
    geom_hline(aes(yintercept=138,color="No microbe"),linewidth=1.25)+
    xlab("Male killing rate")+
    ylab("Egg size (\U003BCM)")+
    scale_color_manual(values=c("#562367","#515797","#338A9F","#7D8F99"))+
    scale_y_continuous(guide=guide_prism_minor(),limits=c(100,650),breaks=c(100,375,650),expand = c(0, 0))+
    scale_x_continuous(guide=guide_prism_minor(),limits=c(0,1.0),breaks=c(0,0.5,1.0),expand = c(0.0, 0.0)))


((Fem_E+ggtitle("Feminization, enhanced growth")+MK_E+ggtitle("Male killing, enhanced growth"))&theme(legend.position = "right",plot.title = element_text(hjust = 0.5,face = "italic",size=14)))+ plot_layout(guides = "collect")

#ggsave("FigS9.svg",height = 183,width=(2*183),units="mm")

# 9. Fig S10 Manipulation vs V -------------------------------------------------------
#Mk vs V
(S10B<-Data_all%>%
   filter(Generation==3999,GR==0.5)%>%
   group_by(MK,V_b)%>%
   summarise(Min=min(EggSize),Max=max(EggSize),Mean=mean(EggSize),Median=median(EggSize))%>%
   pivot_longer(cols = -c(MK,V_b),values_to = "Egg size",names_to = "Metric")%>%
   filter(Metric=="Max")%>%
   ggplot(aes(x=MK,y=V_b))+
   geom_tile(aes(fill=`Egg size`))+
   scale_fill_viridis(option = "mako",na.value = "grey",direction = -1,
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",frame.linewidth=0.5)
   )+
   mytheme+
   labs(x="Male killing rate",y="V parameter")+ scale_x_continuous(expand = c(0.0, 0.0)) +
   scale_y_continuous(expand = c(0.0, 0.0))+
   theme(legend.position = "right",
         axis.text = element_text(color="black"))+
   theme(panel.background = element_rect(fill="grey50")))
#Fem vs V
(S10A<-FData_all%>%
    filter(Generation==3999,GR==0.5)%>%
    group_by(MK,V_b)%>%
    summarise(Min=min(EggSize),Max=max(EggSize),Mean=mean(EggSize),Median=median(EggSize))%>%
    pivot_longer(cols = -c(MK,V_b),values_to = "Egg size",names_to = "Metric")%>%
    filter(Metric=="Max")%>%
    ggplot(aes(x=MK,y=V_b))+
    geom_tile(aes(fill=`Egg size`))+
    scale_fill_viridis(option = "mako",na.value = "grey",direction = -1,
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",frame.linewidth=0.5)
    )+
    mytheme+
    labs(x="Feminization rate",y="V parameter")+ scale_x_continuous(expand = c(0.0, 0.0)) +
    scale_y_continuous(expand = c(0.0, 0.0))+
    theme(legend.position = "right",
          axis.text = element_text(color="black"))+
    theme(panel.background = element_rect(fill="grey50")))


S10A+S10B+ plot_layout(guides = "collect")

#ggsave("FigS10.svg",height = 90,width=(3*90),units="mm")
# 10. Fig S11 Disitributions. ------------------------------------------------------------
# * 10.A S11A --------------------------------------------------------------
#Distribution of egg sizes for male killing that resulted in a life-history transition

(MKDis<-All_inf%>%
  mutate(LHT=ifelse(Transition==NoInf,NoInf,"P to L"))%>%
  filter(LHT=="P to L")%>%
  ggplot(aes(EggSize))+
  geom_histogram(binwidth = 20,fill="#EAEAE9",color="black")+
  mytheme2+
  xlab("Egg size (\U003BCM)")+
  ylab("Number of observations")+
  scale_x_continuous(limits=c(270,930),breaks=c(300,450,600,750,900),expand=c(0,0.044))+
   scale_y_continuous(limits=c(0,14000),breaks=c(0,3500,7000,10500,14000)))

#Distribution of egg sizes for feminizations that resulted in a life-history transition 
(FemDis<-FAll_inf%>%
  mutate(LHT=ifelse(Transition==NoInf,NoInf,"P to L"))%>%
  filter(LHT=="P to L")%>%
  ggplot(aes(EggSize))+
  geom_histogram(binwidth = 20,fill="#EAEAE9",color="black")+
  mytheme2+
  xlab("Egg size (\U003BCM)")+
  ylab("Number of observations")+
    scale_x_continuous(limits=c(270,930),breaks=c(300,450,600,750,900),expand=c(0,0.044))+
    scale_y_continuous(limits=c(0,14000),breaks=c(0,3500,7000,10500,14000)))


(Figs11A<-((FemDis+ggtitle("Feminization, enhanced growth")+MKDis+ggtitle("Male killing, enhanced growth"))&theme(legend.position = "right",plot.title = element_text(hjust = 0.5,face = "italic",size=14))))


# * 10.B S11B Relative egg size ---------------------------------------
(MKDisRel<-All_inf%>%
    mutate(EggVolume=(pi*EggSize^3)/6)%>%
    mutate(EggVolumeUn=(pi*EggSizeUn^3)/6)%>%
    mutate(RelEgg=EggVolume/EggVolumeUn)%>%
   mutate(LHT=ifelse(Transition==NoInf,NoInf,"P to L"))%>%
   filter(LHT=="P to L")%>%
    ggplot(aes(RelEgg))+
    #geom_boxplot()+
    geom_histogram(fill="#EAEAE9",color="black")+
    mytheme2+
    xlab("Relative egg volume")+
    ylab("Observations (x1000)")+
   # scale_y_continuous(limits=c(0,50000),breaks=c(0,25000,50000),labels=c(0.0,25.0,50),expand=c(0,0))+
    scale_x_log10(limits=c(1,120),breaks=c(1,10,100), oob = scales::oob_keep,expand=c(0,0.044)))

(FDisRel<-FAll_inf%>%
    mutate(EggVolume=(pi*EggSize^3)/6)%>%
    mutate(EggVolumeUn=(pi*EggSizeUn^3)/6)%>%
    mutate(RelEgg=EggVolume/EggVolumeUn)%>%
    mutate(LHT=ifelse(Transition==NoInf,NoInf,"P to L"))%>%
    filter(LHT=="P to L")%>%
    ggplot(aes(RelEgg))+
    #geom_boxplot()+
    geom_histogram(fill="#EAEAE9",color="black")+
    mytheme2+
    xlab("Relative egg volume")+
    ylab("Observations (x1,000)")+
    #scale_y_continuous(limits=c(0,50000),breaks=c(0,25000,50000),labels=c(0.0,25.0,50),expand=c(0,0))+
    scale_x_log10(limits=c(1,120),breaks=c(1,10,100), oob = scales::oob_keep,expand=c(0,0.044)))

#Putting it all together
(FigS11B<-((FDisRel+ggtitle("Feminization, enhanced growth")+MKDisRel+ggtitle("Male killing, enhanced growth"))&theme(legend.position = "right",plot.title = element_text(hjust = 0.5,face = "italic",size=14)))+ plot_layout(guides = "collect"))

(Figs11A/FigS11B)
#ggsave("FigS11.svg",height = 183,width=(2*183),units="mm")


All_inf%>%
  mutate(EggVolume=(pi*EggSize^3)/6)%>%
  mutate(EggVolumeUn=(pi*EggSizeUn^3)/6)%>%
  mutate(RelEgg=EggVolume/EggVolumeUn)%>%
  mutate(LHT=ifelse(Transition==NoInf,NoInf,"P to L"))%>%
  filter(LHT=="P to L")%>%
  summarise(mean=mean(RelEgg),max=max(RelEgg))

FAll_inf%>%
  mutate(EggVolume=(pi*EggSize^3)/6)%>%
  mutate(EggVolumeUn=(pi*EggSizeUn^3)/6)%>%
  mutate(RelEgg=EggVolume/EggVolumeUn)%>%
  mutate(LHT=ifelse(Transition==NoInf,NoInf,"P to L"))%>%
  filter(LHT=="P to L")%>%
  summarise(mean=mean(RelEgg),max=max(RelEgg))

