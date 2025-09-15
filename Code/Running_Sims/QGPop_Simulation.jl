#Code by xxx. For questions: xxxx
# Julia script that was used for male killing quantitative genetic population model for:
#"Microbes as manipulators of developmental life-history" 
#xxxxx

#Load up packages

using Distributed
#add the number of procesess i.e. cores being used
addprocs(23)

#Make function
@everywhere using Random, Distributions, StatsBase, GLM, DataFrames,CSV,SharedArrays, LogExpFunctions


#Styan equation
@everywhere function Fitness(E0,S0,sigma, v, tb,tau,B,gr,psi,Fe,v_b)
    # beta0=collision rate constant; estimated by egg size (sigma) and velocity(v)
    sigma2=((pi * (sigma / 1000)^2) / 4)
    vol=(sigma/2000)^3*(4/3)*pi
    # E0<-S0*E*(Females/Males)#Egg density should be proportional to #females
    beta0 = sigma2 * v
    # x=Average number of potential fertilizaing sperm
    x = Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tau)))
    # b=mean number of extra fertilizing sperm that will contact an egg in a time period tb in a population of eggs already contacted
    b = Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tb)))
    prop_mono = 1 - exp(-1 * x) - ((1 - exp(-1 * x) - x * exp(-1 * x)) * (1 - exp(-1 * b)))
    #If no enhanced growth just normal B equation (Eq.6 in text)
    if gr <=0
        Surv=exp(-(B/sigma))
    # with enhanced growth Use Eq. 14 and 15 to modify and change the effective B parameter
    else
        b_p = (3000*sigma)/(v_b+sigma)
        Surv=exp(-((B+b_p)/sigma))
    end
    Feggs=prop_mono*(0.112/vol)*Surv*psi*(1+gr)
    return(Feggs)
  end

 #Custom function to find nearest egg size bin.
@everywhere  function searchsortednearest(a,x)
    idx = searchsortedfirst(a,x)
    if (idx==1); return idx; end
    if (idx>length(a)); return length(a); end
    if (a[idx]==x); return idx; end
    if (abs(a[idx]-x) < abs(a[idx-1]-x))
       return idx
    else
       return idx-1
    end
 end

 ###### For Male killing

@everywhere function Simulation(Sp,v, tb,tau,B,gr,psi,Fe,h2,fine,mu0,sd0,dens0,K,mort,mk,v_b,generation)
#Time, mean, density, sex ratio, mk, gr, b
  dfall=zeros(generation,8)
  #Starting Female density
    fdens = dens0*0.5
    #Starting male density
    mdens= dens0*0.5
    #Starting sex ratio
    sr = mdens/(fdens+ mdens)
    #Possible range of egg sizes
    eggs=collect(10:fine:3000)
    #total population matrix. Columns are eggs, density, and fitness
    pop=[eggs zeros(length(eggs)) zeros(length(eggs))]
    #Starting population. 
    start=collect(mu0-(10*sd0):fine:mu0+(10*sd0))
    #declare distribution
    Dis=Normal(mu0,sd0)
    #calculate probability distribution
    y=pdf.(Dis,start)
    #convert probability distribution to actual densities. 
    ystand=(y/sum(y))*(fdens)
     #Find closest match to minimum in distribution
    minI=findall(pop[:,1].==minimum(start))[1]
    #Find closest match to maximum in distribution
    maxI=findall(pop[:,1].==maximum(start))[1]
    #add the starting population to the population vector. 
    pop[minI:maxI,2]=ystand
    #record starting condition.
    dfall[1,:]=[0,mu0,dens0,sr,mk,gr,B,v_b]
    #Running the simulation for generations
    for i in 2:1:generation
      #Just focus on real population or assuming anything with p<10^-30 can be ignored.
      relpop=pop[pop[:,2].>10^-30,:]
      #focus on real population but log transform
      #1.13/sigma2 is the number eggs
      E0= sum(exp.(log.(0.112 ./(((relpop[:,1]./2000).^3).*(4/3).*pi)) +log.(relpop[:,2])))*(1+gr)
      #calcualte sperm density for fertilization function
      S0=Sp*mdens
      #Calculate fitness of offspring
      fit=Fitness.(E0,S0,relpop[:,1],v,tb,tau,B,gr,psi,Fe,v_b)
      #If fitness was equal or less than zero, assign it a really small number to help with log transformations
      fit[fit .<= 0.0].=10^-30
      #Log transform data due to small numbers.
      relpop[:,3]=log.(fit)
      #calculate mean fitness as sum of proportion of indiviuals for each egg group multiplied by their Fitness
      #Note the log transformation due to undeflow issues.
      meanfitness=sum(exp.(relpop[:,3]+log.(relpop[:,2]./sum(relpop[:,2]))))
      #end if the population crashes
      if meanfitness ==0
        break
      end
      #Calculation selection differential (Eq.16)
      S=(1 ./meanfitness)*cov(exp.(relpop[:,3]),relpop[:,1])
      #Calculate response to selection for next (Eq.17)
      R=h2*S
      #New offspring
      #need to weight by density
      prevmu=sum(exp.(log.(relpop[:,1])+log.(relpop[:,2]./sum(relpop[:,2]))))
      #need to get closest match to offspring mean
      offmu=relpop[searchsortednearest(relpop[:,1],prevmu+R),1]
      #now create new egg density slot
      offMax=minimum([offmu+(10*sd0),3000])
      offMin=maximum([10,offmu-(10*sd0)])
      offeggs=collect(offMin:fine:offMax)
      #now probability distribution
      Dis2=Normal(offmu,sd0)
      off=pdf.(Dis2,offeggs)
      #total production weighted by density
      produc=sum(exp.(relpop[:,3]+log.(relpop[:,2])))
      #weight by density 
      offstand=(off ./sum(off)) .* produc
      #update densities from mortality
      pop[:,2]=pop[:,2].* mort
      #add new densities
      minN=searchsortednearest(pop[:,1],offMin)
      maxN=searchsortednearest(pop[:,1],offMax)
      pop[minN:maxN,2]=((offstand*0.5)+pop[minN:maxN,2])
      mdens=(produc*(1-mk)*0.5)+ (mdens*mort)
      fdens=sum(pop[:,2])
      #calculate new sex ratio
      sr=mdens/(mdens+fdens)
      if (mdens+fdens)>K 
          mdens=K*sr
          fdens=K*(1-sr)
      end
      #convert anything less than zero (underflow) to a very small number.
      pop[pop .<= 0.0].=10^-30
      #Calculate population mean with log functions to account for underflow
      newpopmean=sum(exp.(log.(pop[:,1])+log.(pop[:,2]./sum(pop[:,2]))))
      pop[:,2]=(pop[:,2]./sum(pop[:,2])).*fdens
      #Time, mean, density, sex ratio, mk, gr, b
      dfall[i,:]=[i-1,newpopmean,(mdens+fdens),sr,mk,gr,B,v_b]
  end
  return(dfall)
end

#Sp,sigma, v, tb,tau,B,gr,psi,Fe,h2,fine,mu0,sd0,K,mort,mk,generation
#Sperm per male
@everywhere Spp = 700
#sperm velocity
@everywhere vp= 0.14
#Time to block polyspermy
@everywhere tbp = 1
#Sperm half life
@everywhere taup = 5400
#Fertilization efficiency
@everywhere Fep=0.09444
#Adult mortality
@everywhere mortp=0.9
#Starting density
@everywhere dens0p=1
#Heritability
@everywhere h2p=0.05
#Fineness of bins for egg size distribution
@everywhere finep=0.2
#Settlement constant
@everywhere psip=1000
#Starting standard deviation
@everywhere sd0p=9
#starting average
@everywhere mu0p=180
#Number of generations
@everywhere generationp=4000
#maximum population size.
@everywhere Kp=2.6


# list of male killing rates
#mks=0.84:0.02:0.98
mks = 0:0.06:0.96
#list of growth rates
@everywhere gs= 0:0.5:2
#list of b parameter values
@everywhere bs=300:300:1200
#List of v parameter values
@everywhere vs=100:50:1000

#Collect vector so I can parralelly run combinations
@everywhere things = vec(collect(Iterators.product(gs,bs,vs)))
#change working directory to where I want to save results.
@everywhere cd("xxxx")

#loop throuhg male killing values 
for m in mks
    #parallelel compute other combinations
@everywhere data=SharedArray{Float64}(length(things),8)
@sync @distributed for i in 1:length(things)
    g,b,v=things[i]
    temp=Simulation(Spp,vp,tbp,taup,b,g,psip,Fep,h2p,finep,mu0p,sd0p,dens0p,Kp,mortp,m,v,generationp)
    extract=temp[vec(.!any(isnan.(temp), dims=2)) .& vec(.!all(iszero.(temp), dims=2)), :]
    data[i,:]=extract[end,:]
end
data2=DataFrame(data,["Generation","EggSize","Density","SexRatio","MK","GR","B","V_b"])
CSV.write(string(m,"_IPM_24-5-30.csv"),  string.(data2))
end






#########For Feminization
#######for more indepth code commmnting see above.
@everywhere function SimulationFem(Sp,v, tb,tau,B,gr,psi,Fe,h2,fine,mu0,sd0,dens0,K,mort,Rd,v_b,generation)
  #Time, mean, density, sex ratio, mk, gr, b
    dfall=zeros(generation,8)
      fdens = dens0*0.5
      mdens= dens0*0.5
      sr = mdens/(fdens+ mdens)
      eggs=collect(10:fine:3000)
      #eggs, density, fitness
      pop=[eggs zeros(length(eggs)) zeros(length(eggs))]
      start=collect(mu0-(10*sd0):fine:mu0+(10*sd0))
      #declare distribution
      Dis=Normal(mu0,sd0)
      y=pdf.(Dis,start)
      ystand=(y/sum(y))*(fdens)
      minI=findall(pop[:,1].==minimum(start))[1]
      maxI=findall(pop[:,1].==maximum(start))[1]
      pop[minI:maxI,2]=ystand
      dfall[1,:]=[0,mu0,dens0,sr,Rd,gr,B,v_b]
    #plot(start,ystand)
    #calculating fitness of each bin
  
      for i in 2:1:generation
        relpop=pop[pop[:,2].>10^-30,:]
        #focus on real population but log transform
        #E0,S0,sigma, v, tb,tau,B,gr,psi,Fe
        #1.13/sigma2
        E0= sum(exp.(log.(0.112 ./(((relpop[:,1]./2000).^3).*(4/3).*pi)) +log.(relpop[:,2])))*(1+gr)
        S0=Sp*mdens
        fit=Fitness.(E0,S0,relpop[:,1],v,tb,tau,B,gr,psi,Fe,v_b)
        fit[fit .<= 0.0].=10^-30
        relpop[:,3]=log.(fit)
        #change in mean
        meanfitness=sum(exp.(relpop[:,3]+log.(relpop[:,2]./sum(relpop[:,2]))))
        if meanfitness ==0
          break
        end
        S=(1 ./meanfitness)*cov(exp.(relpop[:,3]),relpop[:,1])
        R=h2*S
        #New offspring
        #need to weight by density
        prevmu=sum(exp.(log.(relpop[:,1])+log.(relpop[:,2]./sum(relpop[:,2]))))
        #need to get closest match to offspring mean
        offmu=relpop[searchsortednearest(relpop[:,1],prevmu+R),1]
        #now create new egg density slot
        offMax=minimum([offmu+(10*sd0),3000])
        offMin=maximum([10,offmu-(10*sd0)])
        offeggs=collect(offMin:fine:offMax)
        #now probability distribution
        Dis2=Normal(offmu,sd0)
        off=pdf.(Dis2,offeggs)
        #total production weighted by density
        produc=sum(exp.(relpop[:,3]+log.(relpop[:,2])))
        #weight by density 
        offstand=(off ./sum(off)) .* produc
        #update densities from mortality
        pop[:,2]=pop[:,2].* mort
        #add new densities
        minN=searchsortednearest(pop[:,1],offMin)
        maxN=searchsortednearest(pop[:,1],offMax)
        pop[minN:maxN,2]=((offstand .* Rd)+pop[minN:maxN,2])
        #Feminzation diverges from male killing here when calculating total densities in the next generation.
        mdens=(produc*(1-Rd))+ (mdens*mort)
        fdens=sum(pop[:,2])
        #sex ratio
        sr=mdens/(mdens+fdens)
        if (mdens+fdens)>K 
            mdens=K*sr
            fdens=K*(1-sr)
        end
        #if(any(pop.<0))
        #  dfall[i,:]=[i-1,-10,-10,sr,R,gr,B]
        #  break
        #end
        pop[pop .<= 0.0].=10^-30
        #print(pop)
        newpopmean=sum(exp.(log.(pop[:,1])+log.(pop[:,2]./sum(pop[:,2]))))
        pop[:,2]=(pop[:,2]./sum(pop[:,2])).*fdens
        #Time, mean, density, sex ratio, mk, gr, b
        dfall[i,:]=[i-1,newpopmean,(mdens+fdens),sr,Rd,gr,B,v_b]
    end
    return(dfall)
  end
#Sp,sigma, v, tb,tau,B,gr,psi,Fe,h2,fine,mu0,sd0,K,mort,mk,generation
@everywhere Spp = 700
@everywhere vp= 0.14
@everywhere tbp = 1
@everywhere taup = 5400
@everywhere Fep=0.09444
@everywhere mortp=0.9
@everywhere dens0p=1
@everywhere h2p=0.05
@everywhere finep=0.2
@everywhere psip=1000
@everywhere sd0p=9
@everywhere mu0p=180
@everywhere generationp=4000
@everywhere Kp=2.6


# list of things
fs=0.81:0.01:0.9
@everywhere gs= 0:0.5:0.5
@everywhere bs=250:10:750
@everywhere vs=250:10:750


@everywhere things = vec(collect(Iterators.product(gs,bs,vs)))
@everywhere cd("xxxx")

for f in fs
  @everywhere data=SharedArray{Float64}(length(things),8)
  @sync @distributed for i in 1:length(things)
      g,b,v=things[i]
      temp=SimulationFem(Spp,vp,tbp,taup,b,g,psip,Fep,h2p,finep,mu0p,sd0p,dens0p,Kp,mortp,f,v,generationp)
      extract=temp[vec(.!any(isnan.(temp), dims=2)) .& vec(.!all(iszero.(temp), dims=2)), :]
      data[i,:]=extract[end,:]
  end
  data2=DataFrame(data,["Generation","EggSize","Density","SexRatio","MK","GR","B","V_b"])
  CSV.write(string(f,"_IPM_FEM_24-5-30.csv"),  string.(data2))
end