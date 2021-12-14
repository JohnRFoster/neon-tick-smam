# This script is the workflow for Jolly-Seber models 
# on the NEON small mammal models


library(tidyverse)

# neon small mammal data
df <- read_csv("Data/allSmallMammals.csv") %>% suppressMessages()
df <- df %>%
  filter(collectYear >= 2016, # subset to training time frame 
         collectYear <= 2020)

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
capture.index <- ch$capture.index
every.day <- ch$every.day
n.days <- length(every.day)
ch <- ch$ch

if(!production){
  ind.seq <- sample(nrow(ch), ind.test, replace = FALSE) %>% sort()
  ch <- ch[ind.seq, ] # for testing
}


# move all recovered dead mice to same observation category "dead"
# need to change unobserved category too
dead <- 4
unobserved <- dead + 1
ch[ch %in% c(dead.u, dead.a, dead.p)] <- dead
ch[ch == not.seen] <- unobserved
ch[is.na(ch)] <- 5
no.captures <- which(apply(ch, 1, function(x) all(x == 5)))
if(length(no.captures) > 0){
  ch <- ch[-no.captures,]
}


# Compute date of first capture
get.first <- function(x) min(which(x != unobserved))
f <- apply(ch, 1, get.first)

if(any(f == ncol(ch))){
  drop.rows <- which(f == ncol(ch))
  ch <- ch[-drop.rows,]
  f <- apply(ch, 1, get.first)
}

if(any(f+1 == ncol(ch))){
  drop.rows <- which(f+1 == ncol(ch))
  ch <- ch[-drop.rows,]
  f <- apply(ch, 1, get.first)
}


source("Functions/known_state_ms.R")
x.known <- known_state_ms(ch, f, unobserved, alive.u, dead)

source("Functions/x_inits.R")
x.inits <- x_inits(x.known, f)

x1 <- rep(1, length(f))
for(i in 1:length(f)){
  x1[i] <- x.known[i,f[i]]  
  if(is.na(x1[i])) x1[i] <- rbinom(1, 1, 0.5) + 1
  
  # ch[i, 1:(f[i]-1)] <- NA
}

# Y1 <- c(
#   length(x1[x1 == 1]) / length(f),
#   length(x1[x1 == 2]) / length(f),
#   length(x1[x1 == 3]) / length(f),
#   0
# )

Y1 <- c(rep(1/3, 3), 0)


data <- list(
  # x = x.known,
  y = as.matrix(ch),
  Y1 = Y1
)

constants <- list(
  n.ind = nrow(ch),
  n.occasions = ncol(ch),
  n.days = n.days,
  no = 5,
  ns = 4,
  f = f,
  capture.index = capture.index
)


# mcmc.matrix <- "NIMBLE_samples.RData"
# mcmc.obj <- "NIMBLE_check.RData"
# 
# init.file <- list.files(site.dir)
# if(mcmc.matrix %in% init.file){
#   load(file.path(site.dir, mcmc.matrix))
#   init.samps <- as.matrix(samples$samples)
# } else if(init.file == mcmc.obj){
#   load(file.path(site.dir, mcmc.obj))
#   init.samps <- as.matrix(save.ls$samples)
# } else {
#   init.samps <- matrix()
# }

# init.df <- read_csv("Data/multiStateBasicParameterSamples.csv") 
# init.samps <- init.df %>% 
#   filter(Site == site,
#          model == "multiStateBasic") %>% 
#   select(-Site, -model)

max.rate <- 0.99999
  
if(FALSE){
  init.mu <- apply(init.samps, 2, mean)
  init.sd <- apply(init.samps, 2, sd)
  
  inits <- function(){list(
    # x = x_inits(x.known, f),
    phi.a = abs(min(rnorm(1, init.mu["phi.a"], init.sd['phi.a']), max.rate)),
    phi.p = abs(min(rnorm(1, init.mu["phi.p"], init.sd['phi.p']), max.rate)),
    psi.ap = abs(min(rnorm(1, init.mu["psi.ap"], init.sd['psi.ap']), max.rate)),     
    psi.pa = abs(min(rnorm(1, init.mu["psi.pa"], init.sd['psi.pa']), max.rate)),     
    p.a = abs(min(rnorm(1, init.mu["p.a"], init.sd['p.a']), max.rate)),
    p.p = abs(min(rnorm(1, init.mu["p.p"], init.sd['p.p']), max.rate)),
    p.d = abs(min(rnorm(1, init.mu["p.d"], init.sd['p.d']), max.rate)),
    p.u = abs(min(rnorm(1, init.mu["p.u"], init.sd['p.u']), max.rate))
  )}
  
} else {
  
  # calculate some inits from ch
  n.alive <- length(ch[ch %in% c(alive.a, alive.p, alive.u)]) # total captured alive
  n.captured <- length(ch[ch != unobserved]) # total captured
  
  p.a.init <- length(ch[ch == alive.a]) / n.captured # proportion captured without ticks
  p.p.init <- length(ch[ch == alive.p]) / n.captured # proportion captured with ticks
  p.d.init <- length(ch[ch == dead]) / n.captured # proportion captured dead
  p.u.init <- length(ch[ch == alive.u]) / n.alive # proportion captured without ticks
  
  
  inits <- function(){list(
    # x = x_inits(x.known, f),
    phi.a = runif(1, 0.98, 1),            
    phi.p = runif(1, 0.98, 1),            
    psi.ap = runif(1, 0.4, 0.6),       
    psi.pa = runif(1, 0.4, 0.6),       
    p.a = runif(1, 0.4, 0.6),
    p.p = runif(1, 0.4, 0.6),
    p.d = runif(1, 0.4, 0.6), 
    p.u = runif(1, 0.4, 0.6)
  )}
  
}

# inits for parameters

message(paste0("Capture history matrix has ", nrow(ch), " rows and ", ncol(ch), " columns"))
# message(paste0("Compiling JAGS with ", n.adapt, " adaptive iterations..."))

if(n.slots > 1){
  source("Functions/run_nimble_parallel2.R")
  cl <- makeCluster(n.slots)
  samples <- run_nimble_parallel(
    cl = cl,
    model = model.code,
    constants = constants,
    data = data,
    inits = inits,
    monitor = monitor,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.ens = 10000,
    thin = thin,
    check.interval = 5,
    max.iter = max.iter,
    save.states = TRUE,
    file.name = file.path(site.dir, "NIMBLE_check.RData")
  )
  stopCluster(cl)
  
  
  message("Writing output...")
  save(samples,
       file = file.path(site.dir, "NIMBLE_samples.RData"))  
  message("--- Done ---")
} else {
  source("R/2.2_multiStateDailySurvivalTransitionNimble.R")
  modelBlock <- nimbleModel(model.code,
                            constants = constants,
                            data = data,
                            inits = inits())
  modelBlock$initializeInfo()
  # 
  # model <- modelBlock$newModel()
  # 
  cModelBlock <- compileNimble(modelBlock)
  # cModel <- compileNimble(model)
  # 
  # 
  mcmcConfBlock <- configureMCMC(cModelBlock,
                                 monitors = monitor,
                                 thin = thin)
  # mcmcConf <- configureMCMC(cModel,
  #                           monitors = monitor,
  #                           thin = thin)
  # 
  # mcmcConfBlock$removeSampler(c("phi.a", "phi.p"))
  # mcmcConfBlock$addSampler(target = c("phi.a", "phi.p"), type = "RW_block")
  # 
  # mcmcConfBlock$removeSampler(c("psi.ap", "psi.pa"))
  # mcmcConfBlock$addSampler(target = c("psi.ap", "psi.pa"), type = "RW_block")
  # 
  # mcmcConfBlock$removeSampler(c("p.a", "p.p", "p.d", "p.u"))
  # mcmcConfBlock$addSampler(target = c("p.a", "p.p", "p.d", "p.u"), type = "RW_block")
  # 
  # 
  # mcmcConfBlock$printSamplers(monitor)
  # 
  # mcmcBuild <- buildMCMC(mcmcConf)
  mcmcBuildBlock <- buildMCMC(mcmcConfBlock)
  # 
  # compMCMC <- compileNimble(mcmcBuild)
  compMCMCBlock <- compileNimble(mcmcBuildBlock)
  # 
  # compMCMC$run(niter = n.iter, reset = FALSE)
  compMCMCBlock$run(niter = n.iter, reset = FALSE)
  # 
  # library(coda)
  # message("Effective Size all RW:")
  # efs <- effectiveSize(as.matrix(compMCMC$mvSamples))
  # print(efs)
  # 
  message("Effective Size with Block:")
  efs <- effectiveSize(as.matrix(compMCMCBlock$mvSamples))
  print(efs)
  # 
  # message("Acceptance Rate all RW:")
  # a.rate <- 1 - rejectionRate(as.mcmc(as.matrix(compMCMC$mvSamples)))
  # print(a.rate)
  # 
  message("Acceptance Rate with Block:")
  a.rate <- 1 - rejectionRate(as.mcmc(as.matrix(compMCMCBlock$mvSamples)))
  print(a.rate)
}





