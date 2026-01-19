# Simulating Capture-Recapture Data from an Open Population.

rm(list = ls())
library(RMark)

# If fitting models to simulated data using RMark, designate the MARK directory. i.e.
# MarkPath <- "D:/Program Files (x86)/MARK"

# Also set working directory for MARK output.

#===========================#============================#===========================#============================#

# Simulates capture histories for an open population compatible with the RMark Package.
# Allows for the possibility of time varying capture, survival, birth, death, immigration and emigration probabilities.


Pop.Simulate <- function(x0,Occasions, Seed = 1,
                         Birth.Rates = 0,Survival.Rates,
                         Emig.Rates = 0,Immig.Rates = 0,
                         Capture.Rates) {
  # x0 is the starting population size (integer)
  # Occasions is the number of observed occasions (integer)
  # Birth.Rates, Survival.Rates, Emig.Rates & Immig.Rates should be a single probability or a vector of length(Occasions-1)
  # Capture.Rats should be a constant probability or vector of length(Occasions)
  
  #Current operation order is: Sample x0 with Capture(1) -> Birth(1),Death(1),Emig(1),Immig(1) ->
  #                            Sample x1 with Capture(2) -> ... ->
  #                            Birth(T-1),Death(T-1),Emig(T-1),Immig(T-1) -> Sample x(T-1) with Capture(T)
  # Where T is the number of Occasions.
  
  set.seed(Seed)
  
  Pop.Size <- rep(x0,length = Occasions)
  Capture.Hist <- matrix(0,nrow = x0,ncol = Occasions)
  Births <- rep(0, length=(Occasions-1))
  Deaths <- rep(0, length=(Occasions-1))
  Emigrants <- rep(0, length = (Occasions-1))
  Immigrants <- rep(0,length = (Occasions-1))
  Living.Indices <- 1:x0
  Highest.Used.Ind <- x0
  
  # Every individual is assigned an index to track their status.
  
  Born.Archive <- list()
  Dead.Archive <- list()
  Emig.Archive <- list()
  Immig.Archive <- list()
  
  if(length(Capture.Rates)==1){
    Capture.Rates <- rep(Capture.Rates, times = Occasions)
  }
  
  if(length(Survival.Rates)==1){
    Survival.Rates <- rep(Survival.Rates,times = Occasions-1)
  }
  
  if(length(Birth.Rates)==1){
    Birth.Rates <- rep(Birth.Rates,times = Occasions-1)
  }
  
  if(length(Emig.Rates)==1){
    Emig.Rates <- rep(Emig.Rates,times = Occasions-1)
  }
  
  if(length(Immig.Rates)==1){
    Immig.Rates <- rep(Immig.Rates,times = Occasions-1)
  }
  # Creates rate vectors if specified as constant.
  
  for(i in 1:(Occasions-1)){
    # Sample
    n.Samp <- rbinom(1,length(Living.Indices),Capture.Rates[i])
    Ind.Samp <- sample(Living.Indices,n.Samp)
    Capture.Hist[Ind.Samp,i] <- 1
    # Observes individuals in the first occasion.
    
    # Births, deaths, and emigration numbers
    Births[i] <- rpois(1,Birth.Rates[i]*Pop.Size[i])
    
    # The number of births is sampled from a Poisson distribution according to the current birth rate and population size.
    if(Births[i]!=0){
      New.Indices = (Highest.Used.Ind+1):(Highest.Used.Ind+Births[i])
      Living.Indices <- c(Living.Indices,New.Indices)
      Born.Archive[[i]] <- New.Indices
      New.Rows <- matrix(0,nrow = Births[i],ncol = Occasions)
      Capture.Hist <- rbind(Capture.Hist,New.Rows)
      Highest.Used.Ind <- Highest.Used.Ind+Births[i]
      
      # If new individuals are born, they are assigned an index and the stored information is updated.
    }
    Deaths[i] <- rbinom(1,length(Living.Indices),1-Survival.Rates[i])
    Dead.Ind <- sample(Living.Indices,Deaths[i])
    Dead.Archive[[i]] <- Dead.Ind
    Living.Indices <- Living.Indices[! Living.Indices %in% Dead.Ind]
    # Deaths are drawn from the binomial distribution and dead individuals are updated.
    
    Emigrants[i] <- rpois(1,Emig.Rates[i]*length(Living.Indices))
    Emig.Ind <- sample(Living.Indices,Emigrants[i])
    Emig.Archive[[i]] <- Emig.Ind
    Living.Indices <- Living.Indices[! Living.Indices %in% Emig.Ind]
    # Individuals leave the population according to the size and emigration rate. Their information is updated.
    
    Immigrants[i] <- rpois(1,Immig.Rates[i]*length(Living.Indices))
    Pop.Size[i+1] <- Pop.Size[i]+Births[i]-Deaths[i]-Emigrants[i]+Immigrants[i]
    
    if(Immigrants[i]!=0){
      New.Indices = (Highest.Used.Ind+1):(Highest.Used.Ind+Immigrants[i])
      Living.Indices <- c(Living.Indices, New.Indices)
      Immig.Archive[[i]] <- New.Indices
      New.Rows <- matrix(0,nrow = Immigrants[i],ncol = Occasions)
      Capture.Hist <- rbind(Capture.Hist,New.Rows)
      Highest.Used.Ind <- Highest.Used.Ind+Immigrants[i]
      # As with births, those joining the population are assigned an index and added.
    }
  }
  # Sample a final time.
  n.Samp <- rbinom(1,length(Living.Indices),Capture.Rates[Occasions])
  Ind.Samp <- sample(Living.Indices,n.Samp)
  Capture.Hist[Ind.Samp,Occasions] <- 1
  
  Non.Zero.ch <- Capture.Hist[rowSums(Capture.Hist) > 0,]
  Output <- as.data.frame(apply(Non.Zero.ch,1,function(row) paste(row,collapse = "")))
  colnames(Output) <- "ch"
  # Prepares the observed capture histories for use with RMark
  return(Output)
}

# Test
Population = Pop.Simulate(x0=100,Occasions=5,
                          Birth.Rates=0.1,Survival.Rates=0.9,
                          Emig.Rates = 0.05,Immig.Rates = 0.1,
                          Capture.Rates = 0.5)


#===========================#============================#===========================#============================#
# Case 1: Constant survival and capture probability (no births, emigration or immigration)

Data1 <- Pop.Simulate(x0 = 50, Occasions =  4, Survival.Rates = 0.95, Capture.Rates = 0.6)


Data1.Processed <- RMark::process.data(Data1, model = "CJS")
Data1.ddl <- RMark::make.design.data(Data1.Processed)

Case1.model <- function(){
  p.dot <- list(formula=~1)
  p.time <- list(formula=~time)
  
  phi.dot <- list(formula=~1)
  phi.time <- list(formula=~time)
  
  # Fits all possible models including capture probabilities and survival rates as both constant and time dependent.
  
  model.list <- RMark::create.model.list("CJS")
  
  Case1.results <- RMark::mark.wrapper(model.list,data = Data1.Processed,
                                       ddl = Data1.ddl,output = FALSE)
  return(Case1.Results)
}

Case1.Results <- Case1.Model()
Case1.Results
# So here the simple model is the best.
constant.p.phi <- case1.results$phi.dot.p.dot
summary(constant.p.phi)
# So phi of 0.977 and p of 0.6 is pretty close to what I simulated.

# 2nd Best model
summary(Case1.Results$Phi.time.p.dot)
# This model suggests that every unit survived in time period 2 with obscene errors.

cleanup(ask = FALSE)

#===========================#============================#===========================#============================#
# Case 2: Linearly decreasing survival probability

Data2 <- Pop.Simulate(x0 = 60, Occasions =  8, Survival.Rates = c(0.85,0.8,0.75,0.7,0.65,0.6,0.55), Capture.Rates = 0.6)

Data2.Processed <- RMark::process.data(Data2, model = "CJS")
Data2.ddl <- RMark::make.design.data(Data2.Processed)

Case2.Model <- function(){
  p.dot <- list(formula=~1)
  p.time <- list(formula=~time)
  p.Time <- list(formula=~Time)
  # Capital T Time corresponds to a linear trend in time.
  
  Phi.dot <- list(formula=~1)
  Phi.time <- list(formula=~time)
  Phi.Time <- list(formula=~Time)
  
  model.list <- RMark::create.model.list("CJS")
  
  Case2.results <- RMark::mark.wrapper(model.list,data = Data2.Processed,
                                       ddl = Data2.ddl,output = FALSE)
  return(Case2.results)
}

Case2.results <- Case2.model()
Case2.results

# Here the two best models are Phi(1)p(Time) and Phi(Time)p(1) with very close AICc. Could be a detection problem

summary(Case2.Results$phi.dot.p.Time)
# Estimate of 0.87 for survival and 0.6,0.4,0.2 for capture.

summary(Case2.Results$phi.Time.p.dot)
# This is the simulated model, Phi estimates given as 0.8,0.72 and 0.0.46 with capture 0.53.

cleanup(ask = FALSE)

#===========================#============================#===========================#============================#
# Case3: Flaws of CJS - survival probability is only 'apparent' and can be confounded with both permanent emigration and mark loss.

Flaw.Test <- Pop.Simulate(x0=100,Occasions = 5, Survival.Rates = 0.8,Emig.Rates = 0.2, Capture.Rates = 0.7)

Flaw.Processed <- process.data(Flaw.Test, model = "CJS")
Flaw.ddl <- make.design.data(flaw.Processed)

Flaw.Model <- function(){
  p.dot <- list(formula=~1)
  p.time <- list(formula=~time)
  
  Phi.dot <- list(formula=~1)
  Phi.time <- list(formula=~time)
  
  model.list <- create.model.list("CJS")
  Flaw.Results <- mark.wrapper(model.list, data = flaw.pr, ddl = flaw.ddl, output = FALSE)
  return(Flaw.Results)
}
Flaw.Results <- Flaw.Model()
Flaw.Results
# So as expected a constant rate model is preferred.
summary(Flaw.Results$Phi.dot.p.dot)
# Estimate of 0.59 which is survival - emigration.

# Emphasises that this model only provides us with apparent survival.

cleanup(ask = FALSE)
#===========================#============================#===========================#============================#
# Case 4: Confounding variables Experiment

Surv.Changed <- Pop.Simulate(x0 = 100,Occasions = 2, Survival.Rates = 0.8, Capture.Rates = 0.5,0.5)
cap.changed <- Pop.Simulate(x0 = 100, occasion = 2, survival.rates = 0.5, Capture.rates = 0.5,0.8)


#confound.test <- function()
{
  p.dot <- list(formula=~1)
  p.time <- list(formula =~time)
  Phi.dot <- list(formula = ~1)
  
  surv.dot <- mark(surv.changed, model = "CJS", model.parameters = list(Phi = Phi.dot, p = p.dot))
  cap.dot <- mark(cap.changed, model = "CJS", model.parameters = list(Phi = Phi.dot, p = p.dot))
}


surv.changed2 <- pop.simulate(x0 = 100, occasions = 3, survival.rates = c(0.9,0.6), capture.rates = c(0.6,0.7,0.8))
cap.changed2 <- pop.simulate(x0 = 100, occasions = 3, survival.rates = c(0.9,0.8), capture.rates = c(0.6,0.7,0.6))

{
  p.dot <- list(formula=~1)
  p.time <- list(formula =~time)
  Phi.dot <- list(formula = ~1)
  
  surv.dot <- mark(surv.changed2, model = "CJS", model.parameters = list(Phi = Phi.dot, p = p.dot))
  surv.time <- mark(surv.changed2, model = "CJS", model.parameters = list(Phi = Phi.dot, p = p.time))
  cap.dot <- mark(cap.changed2, model = "CJS", model.parameters = list(Phi = Phi.dot, p = p.dot))
  cap.time <- mark(cap.changed2, model = "CJS", model.parameters = list(Phi = Phi.dot, p = p.time))
}
# Not so obvious here.

#===========================#============================#===========================#============================#
# Case 5: Constant birth rate, capture probability and survival rate.
# POPAN model allows for births.

Data5 <- Pop.Simulate(x0 = 70, Occasions =  4, Survival.Rates = 0.8, Capture.Rates = 0.65, Birth.Rates = 0.2)

Data5.Processed <- RMark::process.data(Data5, model = "POPAN")
Data5.ddl <- RMark::make.design.data(Data5.Processed)


Data5.Model <- function(pr,ddl) {
  Phi.time <- list(formula = ~ -1 + time, link = "logit")
  # Using -1, +1 or nothing here do the same thing
  Phi.dot <- list(formula = ~ 1, link = "logit")
  
  p.dot = list(formula = ~1, link = "logit")
  p.time <- list(formula = ~ -1 + time, link = "logit")
  
  pent.time <- list(formula = ~ -1+ time, link = "mlogit")
  pent.dot <- list(formula = ~ 1, link = "mlogit")
  
  N.dot <- list(formula = ~ 1, link = "log")
  
  # What is mlogit and why is old timer using -1 everytime
  
  # Constant survival, capture and probability of entry
  Phi.dot.p.dot.pent.dot = mark(pr,ddl, model.parameters = list(Phi = Phi.dot, p = p.dot, pent = pent.dot,N=N.dot))
  
  # Constant capture, p.entry, survival changes with time
  Phi.time.p.dot.pent.dot = mark(pr,ddl, model.parameters = list(Phi = Phi.time, p = p.dot, pent = pent.dot,N=N.dot))
  
  # Constant survival, p.entry and changing capture probability.
  Phi.dot.p.t.pent.dot = mark(pr,ddl, model.parameters = list(Phi = Phi.dot, p = p.time, pent = pent.dot,N=N.dot))
  
  # Constant survival, capture and p.entry changing with time.
  Phi.dot.p.dot.pent.t = mark(pr,ddl, model.parameters = list(Phi = Phi.dot, p = p.dot, pent = pent.time,N=N.dot))
  
  # Survival with time, capture with time, p.entry constant
  Phi.t.p.t.pent.dot = mark(pr,ddl, model.parameters = list(Phi = Phi.time, p = p.time, pent = pent.dot,N=N.dot))
  
  # Survival and p.entry with time, capture constant
  Phi.t.p.dot.pent.t = mark(pr,ddl,model.parameters = list(Phi = Phi.time,p = p.dot,pent = pent.time,N = N.dot))
  
  # Survival constant, capture and p.entry with time
  Phi.dot.p.t.pent.t = mark(pr,ddl, model.parameters = list(Phi = Phi.dot, p = p.time, pent = pent.time,N=N.dot))
  
  # All with time
  Phi.t.p.t.pent.t = mark(pr,ddl, model.parameters = list(Phi = Phi.time,p = p.time,pent = pent.time,N = N.dot))
  
  return(collect.models())
}
case3.results <- data3.model(data3.pr,data3.ddl)
case3.results

summary(case3.results$Phi.dot.p.dot.pent.dot)
# Relatively close to our simulated values.

cleanup(ask = FALSE)