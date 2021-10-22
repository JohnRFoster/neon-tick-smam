# This script is the workflow for Jolly-Seber models 
# on the NEON small mammal models


library(tidyverse)
library(rjags)

load.module("glm")

n.adapt <- 20

# the states coded in ch:
ch.alive.u <- 1      # observed mouse with unknown tick status
ch.alive.a <- 2      # observed mouse without tick attached
ch.alive.p <- 3      # observed mouse with tick attached
ch.dead.u <- 4       # observed mouse with unknown tick status
ch.dead.a <- 5       # observed mouse without tick attached
ch.dead.p <- 6       # observed mouse with tick attached
ch.unobserved <- 7   # unobserved

# neon small mammal data
df <- read_csv("Data/allSmallMammals.csv")

# get capture history matrix
source("Functions/capture_matrix.R")
ch <- capture_matrix("HARV", df)

# move all recoverd dead mice to same observation category "dead"
# need to change unobserved category too
dead <- 4
unobserved <- dead + 1
ch[ch %in% c(ch.dead.u, ch.dead.a, ch.dead.p)] <- dead
ch[ch == ch.unobserved] <- unobserved

ch <- ch[1:25,1:50]

# Compute date of first capture
get.first <- function(x) min(which(x != unobserved))
f <- apply(ch, 1, get.first)

# get the "known" state of each mouse, used as model initial condition
# any capture event between first and last capture that are unobserved
# get inital guess of alive with equal prob. of ticks presence or absence
init_x <- function(ch, f, unknown.code, not.seen.code){
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,] != not.seen.code))
    n2 <- max(which(ch[i,] != not.seen.code))
    if(n1 < n2){ #in case only one capture
      ch[i, n1:n2] <- rbinom(length(n1:n2), 1, 0.5) + 1
    }
    alive.unknown <- which(ch[i,] == unknown.code)
    ch[i, alive.unknown] <- rbinom(length(alive.unknown), 1, 0.5) + 1
    # if(f[i] != 1){
    #   ch[i, 1:(f[i]-1)] <- NA
    # }
    ch[i, 1:(f[i])] <- NA
  }
  return(ch)
}

x.init <- init_x(ch, f, alive.u, unobserved)

# create data list
data <- list(
  y = as.matrix(ch),
  Y1 = c(0.5, 0.5, 0),
  n.ind = nrow(ch),
  n.occasions = ncol(ch),
  first = f
  )

# inits function
inits <- function(){list(
  # x = x.init,
  phi.a = runif(1, 0.9, 1),      # probability of mice survival with ticks absent
  phi.p = runif(1, 0.9, 1),      # probability of mice survival with ticks present
  psi.ap = runif(1, 0.5, 1),     # probability of mice transition from tick absent to tick present
  psi.pa = runif(1, 0.5, 1),     # probability of mice transition from tick present to tick absent
  p.a = runif(1, 0.5, 1),        # probability of observing a mouse with ticks absent
  p.p = runif(1, 0.5, 1),        # probability of observing a mouse with ticks present
  p.d = runif(1, 0, 0.1),       # probability of observing a dead mouse 
  gamma = runif(1, 0.9, 1)   # probability of marking alive mouse unknown tick status when ticks are absent
)}



# load model
source("R/2.1_multiStateBasic.R")
j.model <- jags.model(file = textConnection(model.code),
                      data = data,
                      inits = inits,
                      n.chains = 1,
                      n.adapt = n.adapt)





