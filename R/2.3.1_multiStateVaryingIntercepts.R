# this is state-space multi-state Jolly-Seber model
# this model attempts to fit site affects
# used to fit the NEON small mammal data 

# this model is hierarchical
# all parameters are estimated as site-specific with an across-site mean and precision 
# i.e. varying intercepts 

# this model is the time-invariant case for all parameters
# parameters are estimated daily (i.e. daily survival rate)
# rates are then aggregated across the number of days between
# capture events

# this model has a two daily survival probabilities for mice 
# there is a site hierarchical effect for theses parameters
#   one for ticks absent + site effect
#   one for ticks present + site effect

# two daily transition parameters
#   ticks P -> A + site effect
#   ticks A -> P + site effect

# four recovery parameters - estimated by site (not hierarchical)
#   probability of observing a mouse with ticks absent
#   probability of observing a mouse with ticks present
#   probability of observing a dead mouse
#   probability of marking an alive mouse as unknown tick status

# mice recovered dead do not have a tick P/A status, just dead



# mice either have ticks or not. but a lot of records have an
# "unknown" designation for tick P/A. attempting to capture this
# with the p.u parameter. this model treats the probability 
# of marking a mouse as "unknown" the same regardless of true
# tick P/A state

# True States
# 1 = alive, ticks absent
# 2 = alive, ticks present
# 3 = recently dead
# 4 = dead (absorbing state)

# Observations
# 1 = alive, ticks present
# 2 = alive, ticks absent
# 3 = alive, ticks unknown
# 4 = recovered dead
# 5 = unobserved

model.code <- nimbleCode({
  
  # priors
  phi.a ~ dnorm(0, tau = tau.param)    # mouse survival with ticks absent
  phi.p ~ dnorm(0, tau = tau.param)    # mouse survival with ticks present
  psi.ap ~ dnorm(0, tau = tau.param)   # transition from tick absent to tick present
  psi.pa ~ dnorm(0, tau = tau.param)   # transition from tick present to tick absent
  p.a ~ dnorm(0, tau = tau.param)      # probability of observing a mouse with ticks absent
  p.p ~ dnorm(0, tau = tau.param)      # probability of observing a mouse with ticks present
  p.d ~ dnorm(0, tau = tau.param)      # probability of observing a dead mouse
  p.u ~ dnorm(0, tau = tau.param)    # probability of marking an alive mouse as unknown tick status
  
  # across-site precision priors
  tau.phi.p.site ~ dgamma(0.01, 0.01)
  tau.phi.a.site ~ dgamma(0.01, 0.01)
  tau.psi.pa.site ~ dgamma(0.01, 0.01)
  tau.psi.ap.site ~ dgamma(0.01, 0.01)
  tau.p.a.site ~ dgamma(0.01, 0.01)
  tau.p.p.site ~ dgamma(0.01, 0.01)
  tau.p.d.site ~ dgamma(0.01, 0.01)
  tau.p.u.site ~ dgamma(0.01, 0.01)
  
  for(s in 1:n.site){
    p.a.site[s] ~ dnorm(p.a, tau = tau.p.a.site)
    p.p.site[s] ~ dnorm(p.p, tau = tau.p.p.site)
    p.d.site[s] ~ dnorm(p.d, tau = tau.p.d.site)
    p.u.site[s] ~ dnorm(p.u, tau = tau.p.u.site)
     
    logit(p.a.rate[s]) <- p.a.site[s]
    logit(p.p.rate[s]) <- p.p.site[s]
    logit(p.d.rate[s]) <- p.d.site[s]
    logit(p.u.rate[s]) <- p.u.site[s]
    
    phi.p.site[s] ~ dnorm(phi.p, tau = tau.phi.p.site)
    phi.a.site[s] ~ dnorm(phi.a, tau = tau.phi.a.site)
    psi.pa.site[s] ~ dnorm(psi.pa, tau = tau.psi.pa.site)
    psi.ap.site[s] ~ dnorm(psi.ap, tau = tau.psi.ap.site)
    
    # survival and transition parameters estimated daily
    for(day in 1:n.days){
      logit(Sa.daily[s,day]) <- phi.a.site[s] 
      logit(Sp.daily[s,day]) <- phi.p.site[s]
      logit(Fap.daily[s,day]) <- psi.ap.site[s]
      logit(Fpa.daily[s,day]) <- psi.pa.site[s]
    }
    
    
    # observation matrix - probabilities of O(t) given S(t)
    po[1,1,s] <- p.p.rate[s]
    po[1,2,s] <- 0 
    po[1,3,s] <- p.u.rate[s] * (1 - p.p.rate[s])
    po[1,4,s] <- 0
    po[1,5,s] <- (1 - p.u.rate[s])*(1 - p.p.rate[s])
    po[2,1,s] <- 0 
    po[2,2,s] <- p.a.rate[s]
    po[2,3,s] <- p.u.rate[s] * (1 - p.a.rate[s])
    po[2,4,s] <- 0
    po[2,5,s] <- (1 - p.u.rate[s])*(1 - p.a.rate[s])
    po[3,1,s] <- 0
    po[3,2,s] <- 0
    po[3,3,s] <- 0
    po[3,4,s] <- 1
    po[3,5,s] <- 0
    po[4,1,s] <- 0
    po[4,2,s] <- 0
    po[4,3,s] <- 0
    po[4,4,s] <- 0
    po[4,5,s] <- 1
    
    # aggregate transition and survival between capture occasions
    for(t in 1:(n.occasions-1)){
      
      log(Sa[s, t]) <- sum(log(Sa.daily[s, capture.index[t]:capture.index[t+1]]))
      log(Sp[s, t]) <- sum(log(Sp.daily[s, capture.index[t]:capture.index[t+1]]))
      log(Fap[s, t]) <- sum(log(Fap.daily[s, capture.index[t]:capture.index[t+1]]))
      log(Fpa[s, t]) <- sum(log(Fpa.daily[s, capture.index[t]:capture.index[t+1]]))
      
      # transition matrix - probabilities of state S(t+1) given S(t)
      ps[1,1,s,t] <- Sp[s,t] * (1 - Fpa[s,t])
      ps[1,2,s,t] <- Sp[s,t] * Fpa[s,t]
      ps[1,3,s,t] <- (1 - Sp[s,t]) * p.d.rate[s]
      ps[1,4,s,t] <- (1 - Sp[s,t]) * (1 - p.d.rate[s])
      ps[2,1,s,t] <- Sa[s,t] * Fap[s,t]
      ps[2,2,s,t] <- Sa[s,t] * (1 - Fap[s,t])
      ps[2,3,s,t] <- (1 - Sa[s,t]) * p.d.rate[s]
      ps[2,4,s,t] <- (1 - Sa[s,t]) * (1 - p.d.rate[s])
      ps[3,1,s,t] <- 0
      ps[3,2,s,t] <- 0
      ps[3,3,s,t] <- 0
      ps[3,4,s,t] <- 1
      ps[4,1,s,t] <- 0
      ps[4,2,s,t] <- 0
      ps[4,3,s,t] <- 0
      ps[4,4,s,t] <- 1
    }
  }
  
  # marginalized likelihood
  for(i in 1:n.ind){
    
    y[i, (f[i]+1):n.occasions] ~ dDHMM(
      init = Y1[1:ns],
      probObs = po[1:ns, 1:no, site.ch.index[i]],
      probTrans = ps[1:ns, 1:ns, site.ch.index[i], (f[i]+1):(n.occasions - 1)],
      len = n.occasions - f[i],
      checkRowSums = 0
    )
    
  }
  
})

# variables to monitor
monitor <- c(
  # "x",
  "tau.phi.p.site",
  "tau.phi.a.site",
  "tau.psi.pa.site",
  "tau.psi.ap.site",
  "tau.p.a.site",
  "tau.p.p.site",
  "tau.p.d.site",
  "tau.p.u.site",
  "phi.p.site",
  "phi.a.site",
  "psi.pa.site",
  "psi.ap.site",
  "p.a.site",
  "p.p.site",
  "p.d.site",
  "p.u.site",
  "phi.a",
  "phi.p",
  "psi.ap",
  "psi.pa",
  "p.a",
  "p.p",
  "p.d",
  "p.u"
)



