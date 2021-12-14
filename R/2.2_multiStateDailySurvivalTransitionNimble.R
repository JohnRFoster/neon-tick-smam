# this is state-space multi-state Jolly-Seber model
# used to fit the NEON small mammal data 

# this model is the time-invariant case for all parameters
# parameters are estimated daily (i.e. daily survival rate)
# rates are then aggrigated across the number of days between
# capture events

# this model has a two survival probabilities for mice 
#   one for ticks absent
#   one for ticks present

# two transition parameters
#   ticks P -> A
#   ticks A -> P

# four recovery parameters
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
# 1 = alive, ticks unknown
# 2 = alive, ticks absent
# 3 = alive, ticks present
# 4 = recovered dead
# 5 = unobserved

model.code <- nimbleCode({
  
  # priors
  phi.a ~ dunif(0, 1)    # mouse survival with ticks absent
  phi.p ~ dunif(0, 1)    # mouse survival with ticks present
  psi.ap ~ dunif(0, 1)   # transition from tick absent to tick present
  psi.pa ~ dunif(0, 1)   # transition from tick present to tick absent
  p.a ~ dunif(0, 1)      # probability of observing a mouse with ticks absent
  p.p ~ dunif(0, 1)      # probability of observing a mouse with ticks present
  p.d ~ dunif(0, 1)      # probability of observing a dead mouse
  p.u ~ dunif(0, 1)    # probability of marking an alive mouse as unknown tick status 
  
  # survival and transition parameters estimated daily
  Sa.daily[1:n.days] <- phi.a
  Sp.daily[1:n.days] <- phi.p
  Fap.daily[1:n.days] <- psi.ap
  Fpa.daily[1:n.days] <- psi.pa
  
  # R.a[1:n.occasions] <- p.a  
  # R.p[1:n.occasions] <- p.p  
  # R.d[1:n.occasions] <- p.d
  # R.u[1:n.occasions] <- p.u  
  
  for(t in 1:(n.occasions-1)){
    
    # aggregate transition and survival between capture occasions
    log(Sa[t]) <- sum(log(Sa.daily[capture.index[t]:capture.index[t+1]]))
    log(Sp[t]) <- sum(log(Sp.daily[capture.index[t]:capture.index[t+1]]))
    log(Fap[t]) <- sum(log(Fap.daily[capture.index[t]:capture.index[t+1]]))
    log(Fpa[t]) <- sum(log(Fpa.daily[capture.index[t]:capture.index[t+1]]))

    # transition matrix - probabilities of state S(t+1) given S(t)
    ps[1,1,t] <- Sp[t] * (1 - Fpa[t])
    ps[1,2,t] <- Sp[t] * Fpa[t]
    ps[1,3,t] <- (1 - Sp[t])
    ps[1,4,t] <- 0
    ps[2,1,t] <- Sa[t] * Fap[t]
    ps[2,2,t] <- Sa[t] * (1 - Fap[t])
    ps[2,3,t] <- (1 - Sa[t])
    ps[2,4,t] <- 0
    ps[3,1,t] <- 0
    ps[3,2,t] <- 0
    ps[3,3,t] <- 0
    ps[3,4,t] <- 1
    ps[4,1,t] <- 0
    ps[4,2,t] <- 0
    ps[4,3,t] <- 0
    ps[4,4,t] <- 1
  }
  
  # observation matrix - probabilities of O(t) given S(t)
  po[1,1] <- p.p
  po[1,2] <- 0 
  po[1,3] <- p.u * (1 - p.p)
  po[1,4] <- 0
  po[1,5] <- (1 - p.u)*(1 - p.p)
  po[2,1] <- 0 
  po[2,2] <- p.a
  po[2,3] <- p.u * (1 - p.a)
  po[2,4] <- 0
  po[2,5] <- (1 - p.u)*(1 - p.a)
  po[3,1] <- 0
  po[3,2] <- 0
  po[3,3] <- 0
  po[3,4] <- p.d
  po[3,5] <- 1 - p.d
  po[4,1] <- 0
  po[4,2] <- 0
  po[4,3] <- 0
  po[4,4] <- 0
  po[4,5] <- 1

  
  # likelihood
  for(i in 1:n.ind){

        
    # x[i, first[i]] ~ dcat(Y1[1:ns])

    # for(t in (first[i]+1):n.occasions){

    y[i, (f[i]+1):n.occasions] ~ dDHMM(
      init = Y1[1:ns],
      probObs = po[1:ns, 1:no],
      probTrans = ps[1:ns, 1:ns, (f[i]+1):(n.occasions - 1)],
      len = n.occasions - f[i],
      checkRowSums = 1
    )
      
      # # state process, draw s[t] given s[t-1]
      # x[i,t] ~ dcat(ps[x[i, t-1], t-1, 1:ns])
      # 
      # # observation process: draw o[t] given s[t]
      # # y[i,t] ~ dcat(po[x[i, t], i, t-1, 1:no])
      # y[i,t] ~ dcat(po[x[i, t], 1:no])
    # }
  }
  
})

# variables to monitor
monitor <- c(
  # "x",
  "phi.a",      # probability of mice survival with ticks absent
  "phi.p",      # probability of mice survival with ticks present
  "psi.ap",     # probability of mice transition from tick absent to tick present
  "psi.pa",     # probability of mice transition from tick present to tick absent
  "p.a",      # probability of observing a mouse with ticks absent
  "p.p",      # probability of observing a mouse with ticks present
  "p.d",      # probability of observing a dead mouse
  "p.u"  # probability of marking alive mouse unknown tick status when ticks are absent)
)



