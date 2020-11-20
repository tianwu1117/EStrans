#*********************************************************************************
#The aim of this script is to implement So et al 2011, Pawitan et al 2009, Gillet et al 2018, and 
#Wu & Sham 2020 to transform the genetic effect on the liability scale from log odds ratio, 
#disease prevalence(K), and minor allele frequency (f).
#*********************************************************************************

##So et al 2011 ####
So.trans <- function(f,lnOR,K){
  #genotype frequency
  Paa = (1-f)^2
  PAa = 2*f*(1-f)
  PAA = f^2
  
  #ORs
  OR1 <- exp(lnOR)
  OR2 <- OR1^2
  
  #to calculate faa, RR1 and RR2 by iteration
  faa <- 0.5
  
  repeat{  
    RR1 <- OR1/(1+faa[length(faa)]*(OR1-1))
    RR2 <- OR2/(1+faa[length(faa)]*(OR2-1))
    faa <- rbind(faa, K/(Paa+PAa*RR1+PAA*RR2))
    if(abs(faa[length(faa)]-faa[length(faa)-1])<10^-4){
      break
    }
  }
  faa <- faa[length(faa)]
  
  muaa=0
  fAa= RR1*faa
  fAA= RR2*faa 
  T = qnorm(1-faa) 
  muAa = T-qnorm(1-fAa)
  muAA = T-qnorm(1-fAA)
  mean.all= PAa*muAa+ PAA*muAA
  V.star= Paa*(muaa-mean.all)^2 + PAa*(muAa-mean.all)^2+ PAA*(muAA-mean.all)^2
  hsqd =  V.star/(1+V.star) 
  return(ifelse(lnOR > 0, sqrt(hsqd),-sqrt(hsqd)))
}

##Pawitan et al 2009 ####
Pawi.trans <- function(f,lnOR){
  V.star.pawi <- 2*f*(1-f)*lnOR^2
  Vg.pawi <- V.star.pawi/(V.star.pawi+pi^2/3)
  return(ifelse(lnOR > 0, sqrt(Vg.pawi),-sqrt(Vg.pawi)))
}

##Gillet et al 2018 ####
Gillet.trans <- function(f,lnOR,K){
  lnOR.s <- sqrt(2*f*(1-f))*lnOR
  beta.hat.trans <- qnorm(plogis(log(K/(1-K))+lnOR.s))-qnorm(K)
  return(beta.hat.trans)
}

##Wu&Sham, 2020 ####
LinApp.trans <- function(f,lnOR,K){
  return(K*(1-K)/(dnorm(-qnorm(K)))*lnOR*sqrt(2*f*(1-f)))
}


##Example
f <- 0.05
K <- 0.01
lnOR <- log(1.3)

lnOR.s <- sqrt(2*f*(1-f))*lnOR

So.beta <- So.trans(f,lnOR,K)
Pawitan.beta <- Pawi.trans(f,lnOR)
Gillet.beta <- Gillet.trans(f,lnOR,K)
LinApp.beta <- LinApp.trans(f,lnOR,K)

print(c(So.beta,Pawitan.beta,Gillet.beta,LinApp.beta))