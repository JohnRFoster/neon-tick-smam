# basic tick population model - stage structured matrix model
# all demographic parameters are constant and estimated daily
# process error is estimated via MVN with precisions on the diagonal



library(nimble)
source("Functions/nimble_functions.R")

model.code <- nimbleCode({
  
  ### priors
  phi.l.mu ~ dnorm(0, tau = tau.dem) # larvae survival
  phi.n.mu ~ dnorm(0, tau = tau.dem) # nymph survival
  phi.a.mu ~ dnorm(0, tau = tau.dem)                   # adult survival
  theta.ln ~ dnorm(0, tau = tau.dem)            # larvae -> dormant nymph daily transition 
  theta.na ~ dnorm(0, tau = tau.dem)                 # larvae -> questing nymph daily transition 
  repro.mu ~ T(dnorm(30, tau = 0.01), 0, Inf)     # reproduction
  tau.l ~ dgamma(1.0E-3, 1.0E-3)
  tau.d ~ dgamma(1.0E-3, 1.0E-3)
  tau.n ~ dgamma(1.0E-3, 1.0E-3)
  tau.a ~ dgamma(1.0E-3, 1.0E-3)
  tau.obs ~ dgamma(1.0E-3, 1.0E-3)
  
  ### precision priors
  OMEGA[1,1] <- tau.l
  OMEGA[1,2] <- 0
  OMEGA[1,3] <- 0
  OMEGA[1,4] <- 0
  OMEGA[2,1] <- 0
  OMEGA[2,2] <- tau.d
  OMEGA[2,3] <- 0
  OMEGA[2,4] <- 0
  OMEGA[3,1] <- 0
  OMEGA[3,2] <- 0
  OMEGA[3,3] <- tau.n
  OMEGA[3,4] <- 0
  OMEGA[4,1] <- 0
  OMEGA[4,2] <- 0
  OMEGA[4,3] <- 0
  OMEGA[4,4] <- tau.a
  
  ### first latent process 
  for(p in 1:n.plots){
    x[1, 1, p] ~ dpois(x.init[1])
    x[2, 1, p] ~ dpois(x.init[2])
    x[3, 1, p] ~ dpois(x.init[3])
    x[4, 1, p] ~ dpois(x.init[4])
  }
  
  # convert linear to logit (daily rates)
  logit(l2n) <- theta.ln
  logit(n2a) <- theta.na
  
  ### define parameters
  for(t in 1:N){   # loop over every day in time series
    
    theta.n2a[t] <- if_else_nimble((gdd[t] <= 1000) || (gdd[t] >= 2500),n2a,0)
    lambda[t] <- if_else_nimble((gdd[t] >= 1400) && (gdd[t] <= 2500),repro.mu,0)
    l2n.quest[t] <- if_else_nimble((gdd[t] >= 400) && (gdd[t] <= 2500),1,0)
    
    A[1,1,t] <- phi.l * (1 - l2n)
    A[1,2,t] <- 0
    A[1,3,t] <- 0
    A[1,4,t] <- lambda
    A[2,1,t] <- phi.l * l2n
    A[2,2,t] <- 1 - l2n.quest[t]
    A[2,3,t] <- 0
    A[2,4,t] <- 0
    A[3,1,t] <- 0
    A[3,2,t] <- l2n.quest[t]
    A[3,3,t] <- phi.n * (1 - theta.n2a[t])
    A[3,4,t] <- 0
    A[4,1,t] <- 0 
    A[4,2,t] <- 0
    A[4,3,t] <- phi.n * theta.n2a[t]
    A[4,4,t] <- phi.a
  }
  
  ### Process Model
  for(p in 1:n.plots){
    for(t in 1:(n.occ.plot[p]-1)){
      
      P[1:ns,1:ns,seq.days[p,t,1],p] <- A[1:ns,1:ns,seq.days[p,t,1]]
      
      for(day in 2:df[t]){
        P[1:ns,1:ns,seq.days[p,t,day],p] <- P[1:ns,1:ns,seq.days[p,t,day]+1,p] %*% A[1:ns,1:ns,seq.days[p,t,day],p]
      }
      
      # expected number questing
      Ex[1:ns,t,p] <- P[1:ns,1:ns,start[t],p] %*% x[1:ns,t,p] 
      
      # process error
      px[1:ns,t,p] ~ dmnorm(mean = Ex[1:ns,t,p], prec = OMEGA[1:ns,1:ns])
      
      x[1,t+1,p] <- px[1,t,p]
      x[2,t+1,p] <- max(px[2,t,p], 0)
      x[3,t+1,p] <- px[3,t,p]
      x[4,t+1,p] <- px[4,t,p]
      
    }
    ### Data Model ###
    for(t in 1:N_est){
      
      ## fit the blended model to observed data 
      y[1,t,p] ~ T(dnorm(x[1,t,p], tau = tau.obs), 0, Inf)
      y[3,t,p] ~ T(dnorm(x[3,t,p], tau = tau.obs), 0, Inf)
      y[4,t,p] ~ T(dnorm(x[4,t,p], tau = tau.obs), 0, Inf)
      
    } # t
  }
  
  
  
  
})

monitor <- c("x",
             "phi.l.mu",
             "phi.n.mu",
             "phi.a.mu",
             "theta.ln",
             "theta.na",
             "repro.mu",
             "tau.l",
             "tau.d",
             "tau.n",
             "tau.a",
             "tau.obs")