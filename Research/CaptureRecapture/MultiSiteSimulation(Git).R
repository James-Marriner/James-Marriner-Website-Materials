# Simulating Capture-Recapture Data for Multi-site Models

rm(list = ls())
library(RMark)

# Locate MARK as follows, i.e.
# MarkPath <- "D:/Program Files (x86)/MARK"

setwd("D:/james/Documents/STOR-i/R for STOR-i")

Multi.Site.Simulate <- function(x0,Occasions,Locations,
                                Survival.Rates,Capture.Rates,
                                Transition.Rate.Matrix,
                                Start.n = NA) {
  # x0 - Initial Population Size, positive integer.
  # Occasions - Number of sampling occasions, positive integer.
  # Locations - Number of sampling locations, upto 10 supported. Positive Integer.
  
  # Survival Rates should either be a vector containing constant rates for each location, or,
  # a matrix of dimension L x (O-1) containing survival rates for each location and occasion 
  
  # Capture rates should either be a vector containing constant rates for each location, or,
  # a matrix of dimension L x O containing capture probabilities for each location at each occasion/
  
  # Transition.Rate.Matrix should be a matrix of dimension L x L with entry (i,j) corresponding to the probability
  # of moving from state i to state j.
  # Alternatively, if not constant in time, an array of dimension L x L x (O-1) where
  # entry [i,j,k] corresponds to the probability of moving from j to k in the i'th period. These rows must sum to 1.
  
  # Start.n may be a vector of length L if specifying the initial number of individuals in each Location.
  # If left as NA, the distribution of individuals amongst the population will be multinomial.
  
  # Current operation order is: Sample x0 with C(1) -> S(1), Psi(1) -> Sample x1 with C(2) -> ...
  #                            -> S(T-1), Psi(T-1) -> Sample x(T-1) with C(T) 
  # Where C corresponds to capture, S to Survival, and Psi to transition.
  
  Pop.Size <- rep(x0,length = Occasions+1)
  Capture.Hist <- matrix(0,nrow = x0,ncol = Occasions)
  Deaths <- matrix(0, nrow = Locations, ncol=(Occasions-1))
  # Records the population at each sampling occasion, the capture history and number of deaths at each location between each occasion.
  
  Living.Indices <- 1:x0
  Highest.Used.Ind <- x0
  
  Letters <- c("A","B","C","D","E","F","G","H","I","J")
  # Defines the location names.
  
  Dead.Archive <- vector(mode = "list",length = (Occasions-1))
  # Records the individuals that die by their index.
  
  Location.Indices <- vector(mode = "list",length = Locations)
  # Records the location of each individual.
  
  if(any(is.na(Start.n))){
    Start.Num <- rmultinom(1,size = x0, apply(Transition.Rate.Matrix,2,mean))
  }
  else if(sum(Start.n)!=x0){
    stop("Starting numbers do not sum to starting population size")
  }
  else{
    Start.Num <- t(t(Start.n))
  }
  # If no starting locations are specified, they are multinomially sampled.
  # If the initial distribution does not sum to the initial population, stops the function
  # These are returned as an L x 1 matrix.
  
  if(is.vector(Capture.Rates)){
    Capture.Rates <-  matrix(rep(Capture.Rates,times = Occasions),ncol = Locations, byrow = TRUE)
  }
  # Prepares constant capture rates for use.
  
  if(is.vector(Survival.Rates)){
    Survival.Rates <-  matrix(rep(Survival.Rates,times = (Occasions-1)),ncol = Locations, byrow = TRUE)
  }
  # Prepares constant survival rates for use.
  
  if(is.matrix(Transition.Rate.Matrix)){
    Transition.Rate.Matrix <- array(matrix(rep(Transition.Rate.Matrix,times = (Occasions-1)),ncol =Locations),
                               dim=c(Locations,Locations,(Occasions-1)))
  }
  # Prepares constant transition rates for use.
  
  All.For.Start <- Living.Indices
  for(Area in 1:Locations){
    Location.Temp <- sample(All.For.Start,size = Start.Num[Area,1])
    Location.Indices[[Area]] <- Location.Temp
    All.For.Start <- All.For.Start[! All.For.Start %in% Location.Temp]
  }
  # Randomly assigns individuals to each location based on the starting distribution
  
  
  n.Per.Loc <- rep(0, length = Locations)
  # Stores the number observed in each location
  
  Move.n <- matrix(0,nrow = Locations, ncol= Locations)
  # Records the number of individuals moving between locations
  
  New.Indices <- vector(mode = "list",length = Locations)
  # Will record the location of individuals afrer movement
  
  for(t in 1:(Occasions-1)){
    for(Area in 1:Locations) {
      # Sampling
      n.Per.Loc[Area] <- rbinom(1,length(Location.Indices[[Area]]),Capture.Rates[t,Area])
      Ind.Samp <- sample(Location.Indices[[Area]],n.Per.Loc[Area])
      Capture.Hist[Ind.Samp,t] <- Letters[Area]
      # Individuals are observed in ecah location according to their capture probabilities.
      
      Deaths[Area,t] <- rbinom(1,length(Location.Indices[[Area]]),1-Survival.Rates[t,Area])
      Dead.Ind <- sample(Location.Indices[[Area]],Deaths[Area,t])
      # Individuals die in each location according to their survival rates.
      
      
      if(length(Location.Indices[[Area]])==1 & Deaths[Area,t] == 1){
        Dead.Ind <- Location.Indices[[Area]]
      }
      # Corrects for the behaviour of the sample function if only one individual is present
      
      Dead.Archive[[t]] <- c(Dead.Archive[[t]] ,Dead.Ind)
      Living.Indices <- Living.Indices[! Living.Indices %in% Dead.Ind]
      Location.Indices[[Area]] <- Location.Indices[[Area]][! Location.Indices[[Area]] %in% Dead.Ind]
    }
    Pop.Size[t+1] <- Pop.Size[t]-sum(Deaths[,t])
    # Updates the total number of living individuals
    
    # Movement
    for(Area in 1:Locations){
      Move.n[Area,] <- t(rmultinom(1,size=length(Location.Indices[[Area]]),Transition.Rate.Matrix[Area,,t]))
      # Arranges the number to move in each location i.e.
      # rbind(c(1-to-1,1-to-2,1-to-3),
      #       c(2-to-1,2-to-2,2-to-3),
      #       c(3-to-1,3-to-2,3-to-3))
    }
    for(From in 1:Locations){
      for(To in 1:Locations){
        j.to.k.Ind <- sample(Location.Indices[[From]],as.numeric(Move.n[From,To]))
        
        if(length(Location.Indices[[From]])==1 && Move.n[From,To] == 1){
          j.to.k.Ind <- Location.Indices[[From]]
        }
        # Corrects for the behaviour of the sample function when only one individual may move.
        
        New.Indices[[To]] <- c(New.Indices[[To]],j.to.k.Ind)
        Location.Indices[[From]] <- Location.Indices[[From]][! Location.Indices[[From]] %in% j.to.k.Ind]
      }
      # Records where individuals are moving to and removes them from Location.Indices so they can only move to one location.
    }
    Location.Indices <- New.Indices
    New.Indices <- vector(mode = "list",length = Locations)
    # Updates the locations of the population.
  }
  
  # Final sampling occurence
  for(Area in 1:Locations) {
    n.Per.Loc[Area] <- rbinom(1,length(Location.Indices[[Area]]),Capture.Rates[Occasions,Area])
    Ind.Samp <- sample(Location.Indices[[Area]],n.Per.Loc[Area])
    Capture.Hist[Ind.Samp,Occasions] <- Letters[Area]
  }
  # Samples from each location according to their capture probabilities.

  Zeroes <- rep(0,times = Occasions)
  Non.Zero.ch <- Capture.Hist[!apply(Capture.Hist==Zeroes,1,all),]
  Output <- as.data.frame(apply(Non.Zero.ch,1,function(row) paste(row,collapse = "")))
  colnames(Output) <- "ch"
  
  # Prepares the capture history of all observed individuals for use with RMark
  
  return(Output)
}

multi.test1 <- Multi.Site.Simulate(x0 = 100,Occasions = 4,Locations = 3,
                                   Survival.Rates = c(0.95,0.95,0.95),
                                   Capture.Rates = c(0.55,0.55,0.55),
                                   Transition.Rate.Matrix = rbind(c(0.5,0.3,0.2),
                                                                  c(0.2,0.7,0.1),
                                                                  c(0.3,0.4,0.3)))
