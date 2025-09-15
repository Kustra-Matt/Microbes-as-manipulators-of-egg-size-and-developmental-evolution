#Code by xxxx. For questions: xxxx
# R script that was used for making initial graphs/analysis for 
#"Microbes as manipulators of developmental life-history" 
#xxxxx

#Analysis and graph for Figure S2

# 1.  Load up libraries ---------------------------------------------------
library(tidyverse)
library(npreg)
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
# 2.  Make function to test equality for numbers (avoids rounding  --------
is_equal_tol <- function(x, y, tol = .Machine$double.eps) {
  abs(x - y) < tol
}

# 3.Male killing first ------------------------------------------------------

#make a vector of all male killing
mks<-seq(0.03,0.99,0.03)

#load up data
Data_all<-read_csv("Data/MK_grid.csv")%>%
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


#make empty data frame to fill
results<-data.frame("MK"=mks,peak=0,RelSurv=0)

#For loop that calculates the second relative survival peak (both egg size of the peak and what that peak is)
for(mk in mks){
  DataTy<-Data_all3%>%
    filter((is_equal_tol(MK,mk)&cost==-1.5)|(MK==0&cost==0)|(is_equal_tol(MK,mk)&cost==0)|(MK==0&cost==-1.5))%>%
    mutate(cat=ifelse((is_equal_tol(MK,mk)&cost==-1.5),"EG+MK",ifelse((MK==0&cost==0),"No Inf",ifelse((is_equal_tol(MK,mk)&cost==0),"MK Only","EG Only"))))
  
  
  #just select get EG+Mk and infected
  Dataty2<-DataTy%>%filter(cat!="MK Only",cat!="EG Only") %>%
    drop_na()%>%
    filter(Egg_Size2!=-Inf)
  
  #spline fit uninfected data egg size predicting B parameter
  UMod<-ss(Dataty2$Egg_Size2[Dataty2$cat=="No Inf"],Dataty2$b[Dataty2$cat=="No Inf"],nknots=10)
  
  #spline fit infected data
  InfMod<-ss(Dataty2$Egg_Size2[Dataty2$cat=="EG+MK"],Dataty2$b[Dataty2$cat=="EG+MK"],nknots=10)
  
  #calculate predictions from spline fit
  INF<-predict(InfMod,seq(100,725))%>%
    rename(PredInf=y,InfSE=se,Egg=x)
  
  #calculate predictions from spline fit
  U<-predict(UMod,seq(100,725))%>%
    rename(PredU=y,uSE=se)%>%
    select(-x)
  
  All<-cbind(INF,U)%>%
    mutate(Diff=PredU-PredInf)
  
  #calculate difference in survival of infected relative to uninfected
  All$SurvDif<-exp(-All$PredInf/All$Egg)/exp(-All$PredU/All$Egg)
  
  
  results[is_equal_tol(results$MK,mk),2:3]<-All%>%
    filter(Egg>100)%>%
    filter(SurvDif==max(SurvDif))%>%
    select(Egg,SurvDif)
}



# 4. Feminization ------------------------------------------------------

#Make a vector of all feminziation rates to go through
fems<-seq(0.5,0.995,0.015)

#Read in feminization data
Data_allF<-read_csv("Data/Fem_grid.csv")

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

#make empty data frame to fill for results
resultsF<-data.frame("Feminization"=fems,peak=0,RelSurv=0)


#For loop that calculates the second relative survival peak (both egg size of the peak and what that peak is)
for(f in fems){
  DataTyF<-Data_all3F%>%
    filter((is_equal_tol(Feminization,f)&cost==-1.5)|(Feminization==0.5&cost==0)|(is_equal_tol(Feminization,f)&cost==0)|(Feminization==0.5&cost==-1.5),DIRF!="Either")%>%
    mutate(cat=ifelse((is_equal_tol(Feminization,f)&cost==-1.5),"EG+F",ifelse((Feminization==0.5&cost==0),"No Inf",ifelse((is_equal_tol(Feminization,f)&cost==0),"F Only","EG Only"))))%>%
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
  
  
  resultsF[is_equal_tol(resultsF$Feminization,f),2:3]<-AllF%>%
    filter(Egg>100)%>%
    filter(SurvDif==max(SurvDif))%>%
    select(Egg,SurvDif)
}


# 5.  Plotting Figure S2 -------------------------------------------------------

#please note that these graphs were put together in illustrator to make the final product of Figure S2.

#Plotting feminzaiton peak 
(FemPeak<-ggplot(resultsF,aes(x= Feminization,y=peak))+geom_line(color="#8EB1D0",size=1.5)+
   mytheme+
   labs(x="Feminization rate",y="Egg size (um)")
)

#Plotting relative survival at the peak
(FemSurv<-ggplot(resultsF,aes(x= Feminization,y= RelSurv))+geom_line(color="#8EB1D0",size=1.5,linetype="dashed")+
    mytheme+
    labs(x="Feminization rate",y="Relative survival")
)


#Plotting male killing peak
(MKPeak<-ggplot(results,aes(x= MK,y=peak))+geom_line(color="#45517D",size=1.5)+
    mytheme+
    labs(x="Male killing rate",y="Egg size (um)")
)

#Plotting male killing survival at peak
(MKSurv<-ggplot(results,aes(x= MK,y= RelSurv))+geom_line(color="#45517D",size=1.5,linetype="dashed")+
    mytheme+
    labs(x="Male killing rate rate",y="Relative survival")
)

