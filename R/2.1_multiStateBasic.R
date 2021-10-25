# this is a basic state-space multi-state Jolly-Seber model
# used to fit the NEON small mammal data 

# this model is the time-invariant case for all parameters
# assumes that all capture occasion are equally spaced in time
    # which they absolutely are not!

# this model has a two survival probabilities for mice 
#   one for ticks absent
#   one for ticks present

# mice recovered dead do not have a tick P/A status, just dead

# mice either have ticks or not. but a lot of records have an
# "unknown" designation for tick P/A. attempting to capture this
# with the gamma parameter. this model treats the probability 
# of marking a mouse as "unknown" the same regardless of true
# tick P/A state

# True States
# 1 = alive, ticks absent
# 2 = alive, ticks present
# 3 = recently dead
# 4 = dead (absorbing state)

# Observations
# 1 = alive, ticks unknown
# 2 = alive, ticks absent
# 3 = alive, ticks present
# 4 = recovered dead
# 5 = unobserved

model.code <- " model {
  
  # priors
  phi.a ~ dunif(0, 1)    # mouse survival with ticks absent
  phi.p ~ dunif(0, 1)    # mouse survival with ticks present
  psi.ap ~ dunif(0, 1)   # transition from tick absent to tick present
  psi.pa ~ dunif(0, 1)   # transition from tick present to tick absent
  p.a ~ dunif(0, 1)      # probability of observing a mouse with ticks absent
  p.p ~ dunif(0, 1)      # probability of observing a mouse with ticks present
  p.d ~ dunif(0, 1)      # probability of observing a dead mouse
  gamma ~ dunif(0, 1)    # probability of marking an alive mouse as unknown tick status 
  
  # parameters remain constant for each trapping occasion
  for(t in 1:(n.occasions-1)){
    Sa[t] <- phi.a
    Sp[t] <- phi.p
    Fap[t] <- psi.ap
    Fpa[t] <- psi.pa
    Rd[t] <- p.d
    Ra[t] <- p.a
    Rp[t] <- p.p
    G[t] <- gamma
  }
  
  for(i in 1:n.ind){
    for(t in first[i]:(n.occasions-1)){
    
      # transition matrix - probabilities of state S(t+1) given S(t)
      ps[1,i,t,1] <- Sp[t] * (1 - Fpa[t])
      ps[1,i,t,2] <- Sp[t] * Fpa[t]
      ps[1,i,t,3] <- 1 - Sp[t]
      ps[1,i,t,4] <- 0
      ps[2,i,t,1] <- Sa[t] * Fap[t]
      ps[2,i,t,2] <- Sa[t] * (1 - Fap[t])
      ps[2,i,t,3] <- 1 - Sa[t]
      ps[2,i,t,4] <- 0
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 0
      ps[3,i,t,4] <- 1
      ps[4,i,t,1] <- 0
      ps[4,i,t,2] <- 0
      ps[4,i,t,3] <- 0
      ps[4,i,t,4] <- 1
      
      # observation matrix - probabilities of O(t) given S(t)
      po[1,i,t,1] <- Rp[t] 
      po[1,i,t,2] <- 0 
      po[1,i,t,3] <- G[t] * (1 - Rp[t])
      po[1,i,t,4] <- 0
      po[1,i,t,5] <- 1 - Rp[t]
      po[2,i,t,1] <- 0 
      po[2,i,t,2] <- Ra[t]
      po[2,i,t,3] <- G[t] * (1 - Ra[t])
      po[2,i,t,4] <- 0
      po[2,i,t,5] <- 1 - Ra[t]
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 0
      po[3,i,t,4] <- Rd[t]
      po[3,i,t,5] <- 1 - Rd[t]
      po[4,i,t,1] <- 0
      po[4,i,t,2] <- 0
      po[4,i,t,3] <- 0
      po[4,i,t,4] <- 0
      po[4,i,t,5] <- 1
      
    }
  }
  
  # likelihood
  for(i in 1:n.ind){
  
    x[i, first[i]] <- y[i,first[i]]
    # x[i, first[i]] ~ dcat(Y1)
    
    for(t in (first[i]+1):n.occasions){
    
      # state process, draw s[t] given s[t-1]
      x[i,t] ~ dcat(ps[x[i, t-1], i, t-1,])
      
      # observation process: draw o[t] given s[t]
      y[i,t] ~ dcat(po[x[i, t], i, t-1,])
    }
  }

}"

# variables to monitor
monitor <-
  c("phi.a",      # probability of mice survival with ticks absent
    "phi.p",      # probability of mice survival with ticks present
    "psi.ap",     # probability of mice transition from tick absent to tick present
    "psi.pa",     # probability of mice transition from tick present to tick absent
    "p.a",      # probability of observing a mouse with ticks absent
    "p.p",      # probability of observing a mouse with ticks present
    "p.d",      # probability of observing a dead mouse 
    "gamma"  # probability of marking alive mouse unknown tick status when ticks are absent
    )







