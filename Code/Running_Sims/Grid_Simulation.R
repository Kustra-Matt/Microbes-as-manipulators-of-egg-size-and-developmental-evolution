#Code by xxx
# R script that was used for grid based analysis for:
#"Microbes as manipulators of developmental life-history" 
#xxxx
# 1. Setting up ------------------------------------------------
# * 1.a Loading up required libraries -------------------------------------
library(tidyverse)
library(doParallel)
#library(feather)
#library(viridis)
# * 1.b Clear working space -------------------------------------
rm(list = ls())

# 2. Functions ------------------------------------------------

# * 2.a Number of zygote function -----------------------------------------
# Males is male adult density (Males per m^2)
# sigma is egg size (cross sectional area of the egg; mm^2)
# v is sperm speed (mm/s)
# Fe is fertilization efficiency
# tb is time for polyspermy block (s)
# sr is sex ratio (proportion of population that is male)
# S is sperm density per unit male density (sperm/uL/individual/m^2)
# E is egg density per unit female density (egg/uL/individual/m^2)
# tau is half life of sperm (s)
# c is cost due to infection.
# nif is density of infected females (females per m^2)
# Fertilization dynamics function
# psi is settlment constant (uL/m^2)
number_zygotes <- function(Males, sigma, v, Fe, tb, sr, S, E, tau, nif, c, psi = 1) {
  # beta0=collision rate constant; estimated by egg size (sigma) and velocity(v)
  S0 <- Males * S
  Females <- (Males / sr) - Males
  # E0<-S0*E*(Females/Males)#Egg density should be proportional to #females
  E0 <- E * (Females - nif) + E * (1 - c) * nif
  beta0 <- sigma * v
  tau <- tau
  # x=Average number of potential fertilizaing sperm
  x <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tau)))
  # b=mean number of extr fertilizing sperm that will contact an egg in a time period tb in a population of eggs already contacted
  b <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tb)))
  prop_mono <- 1 - exp(-1 * x) - ((1 - exp(-1 * x) - x * exp(-1 * x)) *
                                    (1 - exp(-1 * b)))
  return(prop_mono * E0 * psi)
}

#number_zygotes(0.9,21,100,10,10,0.5,10,60,100,0.2,0.5)

# * 2.b Population recursion equations ------------------------------------
# Population dynamic functions for feminization
# Nz is density of settled zygotes (individuals/m^2)
# Nif is density of in infected female (individuals/m^2)
# Nuf is density of uninfected females (individuals/m^2)
# Nim is density of infected males (individuals/m^2)
# Num is density of uninfected males (individuals/m^2)
# K is adult carrying capacity (individuals/m^2)
# d is mortality rate of changee
# Ma is base adult mortality
# ml is larval mortality
# R is feminzation rate
# Feminization models
# (nuf/Females) + (1-c)*(nuf/Females)
# infected females
Dens_F_i_t1f <- function(Nz, Nif, Nuf, Nim, Num, K, d, Ma, Ml, R, c) {
  Nt <- Nif + Nuf + Nim + Num
  Nzif <- Nz * (R * (((1 - c) * Nif) / (Nif * (1 - c) + Nuf)))
  return((Nzif * (1 - Ml) + Nif) * (1 - Ma / (1 + exp(-d * (
    Nt + Nz * (1 - Ml) - K / 2
  )))))
}
# uninfected females
Dens_F_u_t1f <- function(Nz, Nif, Nuf, Nim, Num, K, d, Ma, Ml, c) {
  Nt <- Nif + Nuf + Nim + Num
  Nzuf <- Nz * (Nuf / (Nif * (1 - c) + Nuf)) * 0.5
  return((Nzuf * (1 - Ml) + Nuf) * (1 - Ma / (1 + exp(-d * (
    Nt + Nz * (1 - Ml) - K / 2
  )))))
}
# infected males
Dens_M_i_t1f <- function(Nz, Nif, Nuf, Nim, Num, K, d, Ma, Ml, R, c, Mk =
                           0) {
  Nt <- Nif + Nuf + Nim + Num
  Nzim <- Nz * ((1 - R) * ((Nif * (1 - c)) / (Nif * (1 - c) + Nuf))) *
    (1 - Mk)
  return((Nzim * (1 - Ml) + Nim) * (1 - Ma / (1 + exp(-d * (
    Nt + Nz * (1 - Ml) - K / 2
  )))))
}
# uninfected males
Dens_M_u_t1f <- function(Nz, Nif, Nuf, Nim, Num, K, d, Ma, Ml, c) {
  Nt <- Nif + Nuf + Nim + Num
  Nzmf <- Nz * (Nuf / (Nif * (1 - c) + Nuf)) * 0.5
  return((Nzmf * (1 - Ml) + Num) * (1 - Ma / (1 + exp(-d * (
    Nt + Nz * (1 - Ml) - K / 2
  )))))
}

# * 2.c Vectorized equal function -----------------------------------------
# vectorized equal
# tol is how close doubles can be together
# taken from https://stackoverflow.com/questions/35097815/vectorized-equality-testing
is_equal_tol <- function(x, y, tol = .Machine$double.eps) {
  abs(x - y) < tol
}

# 3. Run the simulation for Figure 4A (feminization)-----------------------------------------
# simulate function
# initialize things
# register number of cores for parallel computing
registerDoParallel(23)
# initalize vector of different egg sizes
eggs <- seq(40, 1500, 10)
# convert egg sizes to cross sectional area
eggs <- ((pi * (eggs / 1000)^2) / 4)

eggs2<-(((eggs* 4)/pi)^0.5)*1000

#Different B parameter values
bss<-seq(100,3000,100)

#different feminization values
fss<-seq(0.5,0.995,0.015)
#different male killing values
mss<-seq(0,0.99,0.03)
#Different enhanced growth rates
#coded as negative because originally was considered a cost.
cs<-seq(-3,0,0.1)

# 4. Run the grid based simulations  --------
msims <-
  foreach(e = eggs, .combine = rbind) %:%
  #Replace f=fss with m = mss for male killing values a
  #comment out the whole line for uninfected runs
  foreach(m = mss, .combine = rbind) %:%
  foreach(c = cs, .combine = rbind) %:%
  foreach(
    b_ss = bss,
    .combine = rbind
  ) %dopar% {
    # For mK make r<-0.5 and mk<-m
    # for uninfected runs make r <-0.5 and mk <-0
    r <- 0.5
    mk<-m
    # set flag to false, this lets us know whether population sizes of equilibriated
    flags <- c(FALSE, FALSE, FALSE, FALSE)
    # keep track of any numerical errors
    er <- "None"
    # intilize starting denisties
    nt0 <- c(.1, .1,0,0)
    v<-0.14
    # calculate number of eggs based on egg size. Number taken from same overal reproductive value as H. tub.
    #egg in mm of resident
    e2<-((((e* 4)/pi)^0.5)*1000)
    numegg <- 0.112 / ((e2/2000)^3*(4/3)*pi)
    # Egg size of lower mutant
    eL<-e2-10
    #numegg of lower mutant
    numeggL<-0.112 /((eL/2000)^3*(4/3)*pi)
    #egg size of higher mutant
    eH<-e2+10
    #numegg of higher mutant
    numeggH<-0.112 /((eH/2000)^3*(4/3)*pi)
    #larval mortality rate of resident
    ml <-1-exp(-b_ss/e2)
    #larval mortality or higher mutant
    mlH<-exp(-b_ss/eH)
    #larval mortality of lower mutant
    mlL<-exp(-b_ss/eL)
    sl<-exp(-b_ss/e2)
    
    #sperm resident
    s<-98/v
    vlH<-v+0.01
    spH<-98/vlH
    vlL<-v-0.01
    spL<-98/vlL
    # set generation time to time
    time <- 1
    # Keep going if flags have not stabilized
    while (!all(flags)) {
      nz <-
        number_zygotes(
          Males = nt0[2] + nt0[4],
          sigma = e,
          v = v,
          Fe = 0.09444,
          tb = 1,
          sr = nt0[2] + nt0[4] / sum(nt0),
          S = s,
          E = numegg,
          tau = 5400,
          nif = nt0[1],
          c = c
        )
      nt1 <-
        c(
          Dens_F_i_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            R = r,
            c = c
          ),
          Dens_M_i_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            R = r,
            c = c,
            Mk = mk
          ),
          Dens_F_u_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            c = c
          ),
          Dens_M_u_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            c = c
          )
        )
      # test if there is any change
      flags <- is_equal_tol(nt0, nt1, tol = .Machine$double.eps)
      if (any(is.na(flags))) {
        er <- "NAs"
        break
      }
      # some simulations result in cycles
      if (time > 100000) {
        er <- "Exceeded Time"
        break
      }
      nt0 <- nt1
      time <- time + 1
    }
    #caluclate current density
    dens <- sum(nt0)
    #calculate current sex ratio
    sr <- (nt0[1] + nt0[3]) / dens
    #calculate fitness of current egg size
    resW<-number_zygotes(
      Males = nt0[2] + nt0[4],
      sigma = e,
      v = v,
      Fe = 0.09444,
      tb = 1,
      sr = nt0[2] + nt0[4] / sum(nt0),
      S = s,
      E = numegg,
      tau = 5400,
      nif = nt0[1],
      c = c
    )*sl
    #calculate fitness of greater egg size
    hW<-number_zygotes(
      Males = nt0[2] + nt0[4],
      sigma = ((pi * (eH / 1000)^2) / 4),
      v = v,
      Fe = 0.09444,
      tb = 1,
      sr = nt0[2] + nt0[4] / sum(nt0),
      S =s,
      E = numegg,
      tau = 5400,
      nif = nt0[1],
      c = c
    )/numegg*numeggH*mlH
    #calculate fitness of smaller mutant egg size
    lW<-number_zygotes(
      Males = nt0[2] + nt0[4],
      sigma = ((pi * (eL / 1000)^2) / 4),
      v = v,
      Fe = 0.09444,
      tb = 1,
      sr = nt0[2] + nt0[4] / sum(nt0),
      S = s,
      E = numegg,
      tau = 5400,
      nif = nt0[1],
      c = c
    )/numegg*numeggL*mlL
    #See what direciton evolution goes based on fitness of mutants vs resident.
    if (any(is.na(c(lW,resW,hW)))) {
      er <- "NAs"
      dir<-NA
    }else if(lW >resW &hW > resW){
      dir<-"Either"
    }else if(lW>resW){
      dir<-"Smaller"
    }else if(hW>resW){
      dir<-"Larger"
}else{
     dir<-"Stable"
}
    
    result <-
      data.frame(
        DIRF=dir,
        RESWF=resW,
        LWF=lW,
        HWF=hW,
        Egg_Size = e,
        Egg_Size2=e2,
        Larval_m = ml,
        Egg_number = numegg,
        Feminization = r,
        Time = time,
        Density = dens,
        Sex_ratio = sr,
        Error = er,
        nif0 = nt0[1],
        nim0 = nt0[2],
        cost = c,
        mk = 0,
        b =b_ss
      )
    return(result)
  }
#Save results
write.csv(msims,"GridResults.csv")

