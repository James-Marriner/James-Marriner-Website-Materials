# Simulating Capture-Recapture Data from a Closed Population.

rm(list = ls())
library(RMark)

# If fitting models to simulated data using RMark, designate the MARK directory. i.e.
# MarkPath <- "D:/Program Files (x86)/MARK"

# Also set working directory for MARK output.

#===========================#============================#===========================#============================#
# Simulates capture histories for a closed population compatible with the RMark Package.
# Allows for the possibility of heterogeneous capture probabilities across periods as well as recapture probabilities
# so as to incorporate trap shy effects.


Closed.Sampler <- function(N,Occasions,Cap.Probs,Recap.Probs = 0, Seed = 1){
  # N is the population size.
  # Occasions are the number of sampling/capture occasions
  # Cap.Probs are the probabilities of capture in each occasion. Expected to be a vector of length 1 or length Occasions
  # Recap.Probs are the probabilities of recapture in all subsequent capture occasions. Expected to be a vector of length 1 or length occasions - 1
  
  set.seed(Seed)
  
  Samples <- matrix(0,nrow = N, ncol = Occasions)
  # Stores the capture history of each individual
  
  if(length(Cap.Probs)==1){
    Cap.Probs <- rep(Cap.Probs, times = Occasions)
  }
  if(length(Recap.Probs)==1){
    Recap.Probs <- rep(Recap.Probs,times = (Occasions-1))
  }
  # Creates the necessary capture & recapture probability vectors if only integers given.
  
  Size.1 <- rbinom(1,N,Cap.Probs[1])
  Ind.1 <- sample(1:N,Size.1)
  Samples[Ind.1,1] <- 1
  Sampled.Indices <- Ind.1
  Unsampled.Ind <- (1:N)[! (1:N) %in% Sampled.Indices]
  
  # Draws the number of individuals observed in the first capture occasion then chooses the individuals randomly from the population.
  # Records the individuals observed and removes them from the Unsampled population for future occasions.
  
  for(i in 2:Occasions){
    # Samples using the recapture probabilities for individuals already encountered.
    Size.C <- rbinom(1,length(Sampled.Indices),Recap.Probs[i-1])
    Ind.C <- sample(Sampled.Indices, Size.C)
    Samples[Ind.C,i] <- 1
    
    # Samples from the unseen population using capture probabilities.
    Size.P <- rbinom(1,length(Unsampled.Ind),Cap.Probs[i])
    Ind.P <- sample(Unsampled.Ind,Size.P)
    Sampled.Indices <- c(Ind.P,Sampled.Indices)
    Samples[Ind.P,i] <- 1
    
    # Note: The order of operations ensures that individuals cannot be sampled twice in a single occasion.
  }
  Seen <- Samples[rowSums(Samples) > 0,]
  Output <- as.data.frame(apply(Seen,1,function(row) paste(row,collapse = "")))
  colnames(Output) <- "ch"
  
  # Creates capture histories compatible with RMark for encountered individuals.
  
  return(Output)
}

# Example: Homogeneous Trap shy effect
Trap.Shy <- closed.sampler(100,3,0.6,0.3)
Trap.Processed <- process.data(Trap.Shy,model = "Closed")
Trap.ddl <- make.design.data(Trap.Processed)

# Creates the set of models and first using RMark
Closed.Models <- function(pr,ddl){
  p.dot <- list(formula =~1, share = TRUE)
  # Allows for constant capture probabilities
  
  p.time <- list(formula =~time, share = TRUE)
  # Allows for time varying capture probabilities
  
  p.behav <- list(formula = ~1, share = FALSE)
  # Allows for trap behavioral effects
  
  p.behav.time <- list(formula = ~ time + c, share = FALSE) # Mbt needs work
  # Allows for trap behavioral and time effects
  
  closed.model.list = create.model.list("Closed")
  closed.results = mark.wrapper(closed.model.list,
                                data = pr, 
                                ddl = ddl,
                                output = FALSE)
  return(closed.results)
}

trap.res <- closed.models(trap.pr,trap.ddl)
trap.res

# Worst
M0 <- trap.res$p.dot
summary(M0)

# 2nd Worst
Mbt <- trap.res$p.time
summary(Mbt)

# 2nd Best - Overfitting/big errors I think. Only real difference is capture in period 3.
Mhuh <- trap.res$p.behav.time
summary(Mhuh)

# Best but very close
Mb <- trap.res$p.behav
summary(Mb)
# Pretty close to what I simulated

