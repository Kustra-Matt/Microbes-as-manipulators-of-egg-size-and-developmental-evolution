#Code by xxxx. For questions: xxxx
# R script that was used for making initial graphs/analysis for 
#"Microbes as manipulators of developmental life-history" 
#xxxx

#Graphs specifically generated:
###Fig. S1A
### Fig S5
# 1. Loading up libraries -------------------------------------------------
library(tidyverse)
library(viridis)
library(ggprism)
library(patchwork)
# 2. Setting up themes for plotting --------------------------------------
mytheme <-
  theme_prism() + theme(
    legend.position = "top",
    # this puts legend on the bottom
    axis.title = (element_text(face = "bold")),
    # this makes the axis titles in bold,
    axis.line = element_line(
      color = "black", size =
        0
    ),
    # Makes the axis line black and  thicker
    text = element_text(
      size = 18, face = "bold", color =
        "black"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    )
  ) # makes all the text larger and bold
theme_set(mytheme)


# 3. Equations ------------------------------------------------------------

#Larval Mortality Equation 6
deathl <- function(Mx, b) {
  return((exp(-b / Mx)))
}


#V parameter
Eggfunc<-function(egg,max,v){
  return((max*egg)/(v+egg))
}

#B influence on V parameter equation 15
Eggfunc2 <- function(egg, b,v) {
  B<-(3000*egg)/(v+egg)
  return((exp(-(b+B) / egg)))
}



# 4. Fig S1A Graph ------------------------------------------------------
## First make fake data frame to populate
B_parm<-data.frame("Egg Size"=seq(40,3000,1))

#Populate with calculations
for (b in seq(600,3000,600)){
  B_parm[[paste(b)]]<-deathl(B_parm$Egg.Size,b)
}

#Make Figure S1 showing B parm affect on survival across egg sizes
(S1A<-B_parm%>%
  pivot_longer(-Egg.Size,names_to="Bparm",values_to="Survival")%>%
  ggplot(aes(x = Egg.Size,y=Survival,color=Bparm)) +
  xlab("Egg size (um)") +
  ylab("Larval Survival") +
  theme(text = element_text(size = 15)) +
  labs(color = "Functions") +
  geom_line()+
  scale_y_continuous(limits=c(0,0.9),breaks=c(0,0.3,0.6,0.9))+
  scale_color_viridis(discrete=T,direction = -1))

#colors adjusted in touch up.
#ggsave("FigS1A.SVG",height = 183,width=183,units="mm")


# 5. Fig S5 ---------------------------------------------------------------
#initailizing an empty data framing
Data <- data.frame(Egg = 0)

#Panel A of S5 looking at b modification
(s5a<-ggplot(data = Data, mapping = aes(egg= Egg)) +
  xlab("Egg size (\U003BCM)") +
  ylab("B parameter modifier")+
  scale_x_continuous(guide=guide_prism_minor(),
    limits = c(0, 1500),
    breaks = c(0,500, 1000,1500)) +
  scale_y_continuous(guide=guide_prism_minor(),
    limits = c(0, 3000),
    breaks = c(0,1000, 2000, 3000)) +
  stat_function(
    fun = Eggfunc,
    geom="line",
    args = (list(
      max =3000,v=250)),mapping = aes(color = "v=250"))+
  stat_function(
    fun = Eggfunc,
    geom="line",
    args = (list(
      max =3000,v=500)),mapping = aes(color = "v=500"))+
  stat_function(
    fun = Eggfunc,
    geom="line",
    args = (list(
      max =3000,v=750)),mapping = aes(color = "v=750"))+
    scale_color_manual(values=c("#521E62","#515797","#338A9F")))

#Panel B of S4 looking at survival after modification.
(s5b<-ggplot(data = Data, mapping = aes(egg= Egg)) +
  xlab("Egg size (\U003BCM)") +
  ylab("Larval survival (%)")+
  scale_x_continuous(guide=guide_prism_minor(),
    limits = c(0, 1500),
    breaks = c(0,500, 1000,1500)) +
  scale_y_continuous(guide=guide_prism_minor(),
                     limits = c(0, 0.3),
                     breaks = c(0,0.1, 0.2,0.3))+
  stat_function(
    fun = Eggfunc2,
    geom="line",
    args = (list(
      b =250,v=250)),mapping = aes(color = "v=250"))+
  stat_function(
    fun = Eggfunc2,
    geom="line",
    args = (list(
      b =250,v=500)),mapping = aes(color = "v=500"))+
  stat_function(
    fun = Eggfunc2,
    geom="line",
    args = (list(
      b =250,v=750)),mapping = aes(color = "v=750"))+
    scale_color_manual(values=c("#521E62","#515797","#338A9F")))


((s5a+s5b)&theme(legend.position = "right"))+plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "A",tag_suffix = ")")

#ggsave("FigS5.svg",height = 183,width=(2*183),units="mm")

