# Simulating Capture-Recapture Data arising from Robust designs with potential for temporary emigration.
rm(list = ls())
library(RMark)

# Ensure the location of MARK is defined, i.e.
# MarkPath <- "D:/Program Files (x86)/MARK"
setwd("D:/james/Documents/STOR-i/R for STOR-i")

# Simulating robust data

Robust.Sim <- function(N0, Primary, Secondary,
                       Survival.Rates,Birth.Rates,
                       Capture.Rates,
                       Gamma.p = 1, Gamma.pp = 0) {
  # N0 is the initial population size, positive integer,
  # Primary are the number of primary sampling occasions, positive integer,
  # Secondary are the number of secondary sampling occasions, "",
  
  # Survival.Rates and Birth.Rates are the probabilities of each individual's survival or a new individual's birth.
  # These may be either constant, or a vector of length (P-1) specifiying time dependent rates. 
  
  # Gamma.p defines the probability of remaining in the unobserved population between primary occasions.
  # Gamma.pp defines the probability of moving from the observed to the unobserved population between primary occasions
  
  Total.Occasions = sum(Primary*Secondary)
  Gaps = Primary-1
  
  if(length(Capture.Rates)!=Total.Occasions){
    Capture.Rates = rep(Capture.Rates, length.out = Total.Occasions)
  }
  # Prepares constant capture rates for use.
  
  if(length(Survival.Rates)!=Gaps){
    Survival.Rates = rep(Survival.Rates,length.out = Gaps)
  }
  # Prepares constant survival rates for use.
  
  if(length(Birth.Rates)!=Gaps){
    Birth.Rates = rep(Birth.Rates,length.out = Gaps)
  }
  # Prepares constant birth rates for use.
  
  if(length(Gamma.p)!=Gaps){
    Gamma.p = rep(Gamma.p,length.out = Gaps)
  }
  # Prepares constant gammma.p for use
  
  if(length(Gamma.pp)!=Gaps){
    Gamma.pp = rep(Gamma.pp,length.out = Gaps)
  }
  # Prepares constant Gamma.pp for use
  
  
  Capture.Hist <- matrix(0,nrow = N0, ncol = Total.Occasions)
  # Sum accounts for varying secondary occasions admissible as a vector
  Observable.n = N0
  Unobservable.n = 0
  Unobservable.Indices = c()
  Superpop.n = N0
  Observable.Indices <- 1:N0
  Used.Indices <- 1:N0
  Occasion <- 1
  # We assume the secondary population is closed.
  for(i in 1:Gaps){
    # Sampling within closed secondary periods. This block will change depending on heterogeneity, time and behavioural effects
    for(j in 1:Secondary){
      n.Samp = rbinom(1,Observable.n, Capture.Rates[Occasion])
      Ind.Samp = sample(x = Observable.Indices, size = n.Samp, replace = FALSE)
      # Samples from the observable population according to capture probabilities.
      
      if(length(Observable.Indices) == 1 & n.Samp == 1){
        Ind.Samp = Observable.Indices
      }
      # To stop the undesirable feature of sample() in the case of sampling a length 1 vector
      
      Capture.Hist[Ind.Samp, Occasion] = 1
      Occasion = Occasion + 1
      # Records those seen in the first primary sampling occaison
    }
    # Births. This block will change depending if births are different for the unobserved and observed groups.
    Obs.Births = rpois(1,Observable.n*Birth.Rates[i])
    Unobs.Births = rpois(1,Unobservable.n*Birth.Rates[i])
    if(Obs.Births > 0){
      # Need to change birth.indices
      Obs.Birth.Indices = (Superpop.n+1):(Superpop.n+Obs.Births)
      New.Rows = matrix(0, nrow = Obs.Births, ncol = Total.Occasions)
      Capture.Hist = rbind(Capture.Hist, New.Rows)
      # Creates capture histories for the new observed individuals
      
      Used.Indices = c(Used.Indices, Obs.Birth.Indices)
      Superpop.n = Superpop.n + Obs.Births
      Observable.n = Observable.n + Obs.Births
      Observable.Indices = c(Observable.Indices, Obs.Birth.Indices)
      # Updates the indices used within the observable and super populations
      
    }
    if(Unobs.Births > 0){
      # Need to change birth.indices
      Unobs.Birth.Indices = (Superpop.n+1):(Superpop.n+Unobs.Births)
      New.Rows = matrix(0, nrow = Unobs.Births, ncol = Total.Occasions)
      Capture.Hist = rbind(Capture.Hist, New.Rows)
      # Creates capture histories for the new unobservable individuals
      
      Used.Indices = c(Used.Indices, Unobs.Birth.Indices)
      Superpop.n = Superpop.n + Unobs.Births
      Unobservable.n = Unobservable.n + Unobs.Births
      Unobservable.Indices = c(Unobservable.Indices, Unobs.Birth.Indices)
      # Updates the indices used within the unobservable and super populations
    }
    
    # Deaths.
    n.Obs.Deaths = rbinom(1,Observable.n, 1-Survival.Rates[i])
    Obs.Death.Indices = sample(x = Observable.Indices, size = n.Obs.Deaths, replace = FALSE)
    # Randomly samples deaths from the observable population
    
    if(length(Observable.Indices) == 1 & n.Obs.Deaths == 1){
      Obs.Death.Indices = Observable.Indices
    }
    # To stop the undesirable feature of sample() in the case of sampling a length 1 vector
    
    Observable.Indices = Observable.Indices[! Observable.Indices %in% Obs.Death.Indices]
    Observable.n = Observable.n - n.Obs.Deaths
    # Updates the living individuals.
    
    n.Unobs.Deaths = rbinom(1,Unobservable.n, 1-Survival.Rates[i])
    Unobs.Death.Indices = sample(x = Unobservable.Indices, size = n.Unobs.Deaths, replace = FALSE)
    # Randomly samples deaths from the unobservable population
    
    if(length(Unobservable.Indices) == 1 & n.Unobs.Deaths == 1){
      Unobs.Death.Indices = Unobservable.Indices
    }
    # To stop the undesirable feature of sample() in the case of sampling a length 1 vector
    
    Unobservable.Indices = Unobservable.Indices[! Unobservable.Indices %in% Unobs.Death.Indices]
    Unobservable.n = Unobservable.n - n.Unobs.Deaths
    # Updates the unobservable population for only living individuals.
    
    # Temporary Emigration stage
    # We account for 3 possible cases: random, Markovian, and none
    # Gamma.p is (U,U)
    # Gamma.pp is (O,U)
    if(!(Gamma.pp[i] == 0 & Gamma.p[i] ==  1)){
      # i.e if we have random or Markovian temporary emigration
      
      OU.n = rbinom(1,Observable.n,Gamma.pp[i])
      OU.Ind = sample(Observable.Indices, size = OU.n, replace = FALSE)
      # Samples individuals moving from observed to unobserved.
    
      if(length(Observable.Indices) == 1 & OU.n == 1){
        OU.Ind = Observable.Indices
      }
      # To deal with the behaviour of the sample function.
      
      OO.Ind = Observable.Indices[! Observable.Indices %in% OU.Ind]
      # Unsampled individuals remain observable
      
      UU.n = rbinom(1,Unobservable.n,Gamma.p[i])
      UU.Ind = sample(Unobservable.Indices, size = UU.n, replace = FALSE)
      # Samples those individuals who will remain unobserved
      
      if(length(Unobservable.Indices) == 1 & UU.n == 1){
        UU.Ind = Unobservable.Indices
      }
      # Dealing with the behaviour of the sample function.
      
      UO.Ind = Unobservable.Indices[! Unobservable.Indices %in% UU.Ind]
      # Remaining individuals become observable.
      
      Observable.Indices = c(OO.Ind, UO.Ind)
      Observable.n = length(Observable.Indices)
      Unobservable.Indices = c(OU.Ind, UU.Ind)
      Unobservable.n = length(Unobservable.Indices)
      # The observable and unobservable populations are hence updated.
    }
  }
  #  Final sampling occasion with no population update afterwards
  for(j in 1:Secondary){
    n.Samp = rbinom(1,Observable.n, Capture.Rates[Occasion])
    Ind.Samp = sample(x = Observable.Indices, size = n.Samp, replace = FALSE)
    
    if(length(Observable.Indices) == 1 & n.Samp == 1){
      ind.samp = Observable.Indices
    }
    # To stop the undesirable feature of sample() in the case of sampling a length 1 vector
    
    Capture.Hist[Ind.Samp, Occasion] = 1
    Occasion = Occasion + 1
  }
  Zeroes <- rep(0,times = sum(Primary*Secondary))
  Non.Zero.ch <- Capture.Hist[!apply(Capture.Hist==Zeroes,1,all),]
  Output <- as.data.frame(apply(Non.Zero.ch,1,function(row) paste(row,collapse = "")))
  colnames(Output) <- "ch"
  # Removes unobserved individuaks and prepares the capture history for use by RMark
  
  return(Output)
}

Test = Robust.Sim(N0 = 100,Primary = 4,Secondary = 3,
                  Survival.Rates = 0.95,Capture.Rates = 0.6,
                  Birth.Rates = 0.1, Gamma.p = 0.2,Gamma.pp = 0.2)
# Working as intended for no behavioral dependence.Birth and death equal for obs/unobs. Can newborns instantly migrate?- could  change the order to fix.
