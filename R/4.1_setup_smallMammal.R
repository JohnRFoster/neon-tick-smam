# This script is the workflow for Jolly-Seber models 
# on the NEON small mammal models


library(tidyverse)
library(rjags)

load.module("glm")

# neon small mammal data
df <- read_csv("Data/allSmallMammals.csv")

# get capture history matrix
source("Functions/capture_matrix.R")
ch <- capture_matrix(site, df)

if(!production){
  ch <- ch[1:25,1:50] # for testing
}

# the states coded in ch:
ch.alive.u <- 1      # observed mouse with unknown tick status
ch.alive.a <- 2      # observed mouse without tick attached
ch.alive.p <- 3      # observed mouse with tick attached
ch.dead.u <- 4       # observed mouse with unknown tick status
ch.dead.a <- 5       # observed mouse without tick attached
ch.dead.p <- 6       # observed mouse with tick attached
ch.unobserved <- 7   # unobserved

# move all recovered dead mice to same observation category "dead"
# need to change unobserved category too
dead <- 4
unobserved <- dead + 1
ch[ch %in% c(ch.dead.u, ch.dead.a, ch.dead.p)] <- dead
ch[ch == ch.unobserved] <- unobserved

# Compute date of first capture
get.first <- function(x) min(which(x != unobserved))
f <- apply(ch, 1, get.first)


known_state_ms <- function(ms, notseen, unknown, present, absent, dead){
  state <- ms
  state[state == notseen] <- NA
  state[state == unknown] <- NA
  state[state == absent] <- 1
  state[state == present] <- 2
  state[state == dead] <- 3
  return(state)
}

# x.init <- init_x(ch, f, ch.alive.u, unobserved)

# create data list
data <- list(
  x = known_state_ms(ch, unobserved, ch.alive.u, ch.alive.p, ch.alive.a, dead),
  y = as.matrix(ch),
  Y1 = c(0.5, 0.5, 0, 0),
  n.ind = nrow(ch),
  n.occasions = ncol(ch),
  first = f
  )

# latent state inits, still needs debugging, get jags error
x_inits <- function(ch, f, notseen, dead){
  states <- ch
  v <- which(states == notseen)
  states[-v] <- NA
  states[v] <- rbinom(length(v), 1, 0.5) + 1
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,] != notseen))
    n2 <- max(which(ch[i,] != notseen))
    n3 <- which(ch[i,] == dead)
    if(n1 < n2){
      states[i, n1:n2] <- rbinom(length(n1:n2), 1, 0.5) + 1
      states[i, n1] <- NA
      states[i, n2] <- NA
    }
    if(length(n3) == 1){
      states[i, (n3+1):nrow(ch)] <- dead + 1
    } 
    states[i, 1:f[i]] <- NA
  }
  return(states)
}

# inits for parameters
inits <- function(){list(
  # x = x_inits(ch, f, unobserved, dead),
  phi.a = runif(1, 0.9, 1),      # probability of mice survival with ticks absent
  phi.p = runif(1, 0.9, 1),      # probability of mice survival with ticks present
  psi.ap = runif(1, 0.5, 1),     # probability of mice transition from tick absent to tick present
  psi.pa = runif(1, 0.5, 1),     # probability of mice transition from tick present to tick absent
  p.a = runif(1, 0.5, 1),        # probability of observing a mouse with ticks absent
  p.p = runif(1, 0.5, 1),        # probability of observing a mouse with ticks present
  p.d = runif(1, 0, 0.1),        # probability of observing a dead mouse 
  gamma = runif(1, 0.9, 1)       # probability of marking alive mouse unknown tick status when ticks are absent
)}

j.model <- jags.model(file = textConnection(model.code),
                      data = data,
                      inits = inits,
                      n.chains = 1,
                      n.adapt = n.adapt)

for(i in 1:n.loops){
  
  jags.out <- coda.samples(j.model,
                           variable.names = monitor,
                           thin = thin,
                           n.iter = n.iter)
  
  loop.time <- Sys.time() - loop.start
  cat(n.iter, "iterations in\n")
  print(loop.time)
  
  ## split output
  out <- list(params = NULL, predict = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
  out <- list(out = out, j.model = return$j.model)
  
  iter <- paste0(site, "_", chain.num, "_", i, ".RData")
  
  save(out, 
       monitor,
       file = file.path(site.dir, iter))
  
}




