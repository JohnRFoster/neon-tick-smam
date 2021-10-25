# This script is the workflow for Jolly-Seber models 
# on the NEON small mammal models


library(tidyverse)
library(rjags)

load.module("glm")

# neon small mammal data
df <- read_csv("Data/allSmallMammals.csv") %>% suppressMessages()

# get capture history matrix
source("Functions/capture_matrix.R")
ch <- capture_matrix(site, df)


# the states coded in ch:
alive.u <- ch$alive.u      # observed mouse with unknown tick status
alive.a <- ch$alive.a      # observed mouse without tick attached
alive.p <- ch$alive.p      # observed mouse with tick attached
dead.u <- ch$dead.u       # observed mouse with unknown tick status
dead.a <- ch$dead.a       # observed mouse without tick attached
dead.p <- ch$dead.p       # observed mouse with tick attached
not.seen <- ch$unobserved   # unobserved
ch <- ch$ch

if(!production){
  ch <- ch[1:25,1:50] # for testing
}

### TO DO
# 1. FIX MODEL TO MATCH CODE SWITCH
# 2. REMOVE DEATH 

# move all recovered dead mice to same observation category "dead"
# need to change unobserved category too
dead <- 4
unobserved <- dead + 1
ch[ch %in% c(dead.u, dead.a, dead.p)] <- dead
ch[ch == not.seen] <- unobserved

# Compute date of first capture
get.first <- function(x) min(which(x != unobserved))
f <- apply(ch, 1, get.first)

# are there any dead animals?
mice.found.dead <- which(apply(ch, 1, function(x) any(x == dead))) 

# supply known states
# anything coded as unobserved or unknown tick status get NA
# anything after a death gets absorbing state code
known_state_ms <- function(ms, f, notseen, unknown, dead, dead.ind){
  state <- ms
  latent.dead <- dead - 1
  state[state == notseen] <- NA
  state[state == unknown] <- NA
  state[state == dead] <- latent.dead
  if(length(dead.ind) >= 1){
    for(i in 1:length(dead.ind)){
      n3 <- which(ms[dead.ind[i],] == dead)
      state[dead.ind[i], (n3+1):ncol(ms)] <- latent.dead + 1
    }
  }
  # state[,1] <- NA
  return(state)
}

# latent state initial condition
# anything after first capture and NA in x.known gets a guess
x_inits <- function(ks, f){
  states <- ks
  v <- which(is.na(states))
  states[-v] <- NA
  states[v] <- rbinom(length(v), 1, 0.5) + 1
  for (i in 1:nrow(ks)){
    if(f[i] > 1) states[i, 1:(f[i]-1)] <- NA
  }
  return(states)
}

x.known <- known_state_ms(ch, f, unobserved, alive.u, dead, mice.found.dead)
x.inits <- x_inits(x.known, f)

# create data list
data <- list(
  x = x.known,
  y = as.matrix(ch),
  Y1 = c(0.5, 0.5, 0, 0),
  n.ind = nrow(ch),
  n.occasions = ncol(ch),
  first = f
)

# calculate some inits from ch
n.alive <- length(ch[ch %in% c(alive.a, alive.p, alive.u)]) # total captured alive
n.captured <- length(ch[ch != unobserved]) # total captured

p.a.init <- length(ch[ch == alive.a]) / n.alive # proportion captured without ticks
p.p.init <- length(ch[ch == alive.p]) / n.alive # proportion captured with ticks
p.d.init <- length(ch[ch == dead]) / n.captured # proportion captured dead
gamma.init <- length(ch[ch == alive.u]) / n.alive # proportion captured without ticks

# inits for parameters
inits <- function(){list(
  x = x_inits(x.known, f),
  phi.a = runif(1, 0.9, 1),      # probability of mice survival with ticks absent
  phi.p = runif(1, 0.9, 1),      # probability of mice survival with ticks present
  psi.ap = runif(1, 0, 1),     # probability of mice transition from tick absent to tick present
  psi.pa = runif(1, 0, 1),     # probability of mice transition from tick present to tick absent
  p.a = abs(min(jitter(p.a.init), 1)),        # probability of observing a mouse with ticks absent
  p.p = abs(min(jitter(p.p.init), 1)),        # probability of observing a mouse with ticks present
  p.d = abs(min(jitter(p.d.init), 1)),        # probability of observing a dead mouse 
  gamma = abs(min(jitter(gamma.init), 1))       # probability of marking alive mouse unknown tick status when ticks are absent
)}

message(paste0("Capture history matrix has ", nrow(ch), " rows and ", ncol(ch), " columns"))
message(paste0("Compiling JAGS with ", n.adapt, " adaptive iterations..."))

compile.start <- Sys.time()
j.model <- jags.model(file = textConnection(model.code),
                      data = data,
                      inits = inits,
                      n.chains = 1,
                      n.adapt = n.adapt)
message("JAGS model compiled. Compile time:")
print(Sys.time()-compile.start)

for(i in 1:n.loops){
  loop.start <- Sys.time()
  jags.out <- coda.samples(j.model,
                           variable.names = monitor,
                           thin = thin,
                           n.iter = n.iter)
  
  loop.time <- Sys.time() - loop.start
  message(paste0(n.iter, " iterations in"))
  print(loop.time)
  
  total.iter <- i*n.iter
  message(paste0("For a total of ", total.iter, " iterations"))
  
  total.time <- Sys.time() - start.time
  message("Total run time of")
  print(total.time)
  
  ## split output
  out <- list(params = NULL, predict = NULL)
  mfit <- as.matrix(jags.out, chains = TRUE)
  pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
  out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
  out <- list(out = out, jags.model = j.model)
  
  iter <- paste0(site, "_", chain, "_", i, ".RData")
  
  save(out, 
       monitor,
       file = file.path(site.dir, iter))
  
}




