library(tidyverse)
library(rjags)
library(nimble)
library(parallel)

n.slots <- Sys.getenv("NSLOTS") %>% as.numeric()

message(paste("N slots:", n.slots))

n.adapt <- 500
n.burnin <- 5000
n.iter <- 50000
n.ens <- 10000
n.chains <- 3

out.dir <- "Test/multiStateBasic"
if(!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

phi.a <- 0.6    # mouse survival with ticks absent
phi.p <- 0.9    # mouse survival with ticks present
psi.ap <- 0.7   # transition from tick absent to tick present
psi.pa <- 0.2   # transition from tick present to tick absent
p.a <- 0.8      # probability of observing a mouse with ticks absent
p.p <- 0.6      # probability of observing a mouse with ticks present
p.d <- 0.1      # probability of observing a dead mouse
p.u <- 0.4
n.occasions <- 30
n.states <- 4
n.obs <- 5
marked <- matrix(0, ncol = n.states, nrow = n.occasions)
marked[,1] <- rep(5, n.occasions) # Releases in study area

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time

# 1. State process matrix
totrel <- sum(marked)*(n.occasions - 1)
PSI.STATE<- array(NA, dim=c(n.states, n.states, totrel, n.occasions - 1))
for (i in 1:totrel){
  for (t in 1:(n.occasions - 1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phi.p*(1-psi.pa), phi.p*psi.pa    , (1-phi.p)*p.d, (1-phi.p)*(1-p.d),
      phi.a*psi.ap    , phi.a*(1-psi.ap), (1-phi.a)*p.d, (1-phi.a)*(1-p.d),
      0               , 0               , 0      , 1,
      0               , 0               , 0      , 1), 
      nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      p.p, 0  , p.u*(1-p.p), 0  , (1-p.u)*(1-p.p),
      0  , p.a, p.u*(1-p.a), 0  , (1-p.u)*(1-p.a),
      0  , 0  , 0          , 1  , 0,
      0  , 0  , 0          , 0  , 1),
      nrow = n.states, byrow = TRUE)
  } #t
} #i


# Execute simulation function
source("Functions/simulate_ms.R")
sim <- simulate_ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

# Compute date of first capture
get.first <- function(x) min(which(x != 5))
f <- apply(CH, 1, get.first)

ch <- CH

unobserved <- 5
alive.u <- 3
dead <- 4

source("Functions/known_state_ms.R")
x.known <- known_state_ms(ch, f, unobserved, alive.u, dead)

source("Functions/x_inits.R")
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

alive.p <- 1
alive.a <- 2
alive.u <- 3

# inits for parameters
inits <- function(){list(
  x = x_inits(x.known, f),
  phi.a = runif(1, 0, 1),    
  phi.p = runif(1, 0, 1),    
  psi.ap = runif(1, 0, 1),   
  psi.pa = runif(1, 0, 1),   
  p.a = runif(1, 0, 1),      
  p.p = runif(1, 0, 1),      
  p.d = runif(1, 0, 1),      
  p.u = runif(1, 0, 1)       
)}

dataN <- list(
  x = x.known,
  y = as.matrix(ch),
  Y1 = c(0.5, 0.5, 0, 0)
)

constantsN <- list(
  n.ind = nrow(ch),
  n.occasions = ncol(ch),
  no = n.obs,
  ns = n.states,
  first = f
)

source("R/2.1.1_multiStateBasicNimble.R")
# message("Build model")
# model <- nimbleModel(
#   model.code,
#   constants = constantsN,
#   data = dataN,
#   inits = inits()
# )
# message("Check initialization")
# model$initializeInfo()
# 
# message("Compile model")
# cModel <- compileNimble(model)
# 
# message("Configure MCMC")
# mcmcConf <- configureMCMC(cModel, monitors = monitor)
# 
# message("Build MCMC")
# mcmcBuild <- buildMCMC(mcmcConf)
# 
# message("Test run")
# mcmcBuild$run(1)

# message("Compile MCMC")
# mcmcComp <- compileNimble(mcmcBuild)


source("Functions/run_nimble_parallel2.R")
cl <- makeCluster(n.slots)
samples <- run_nimble_parallel(
  cl = cl,
  model = model.code,
  constants = constantsN,
  data = dataN,
  inits = inits,
  monitor = monitor,
  n.iter = n.iter,
  n.burnin = n.burnin,
  n.ens = n.ens,
  # max.iter = 50000,
  thin = 10,
  check.interval = 4,
  check.params.only = TRUE,
  file.name = "Test/multiStateBasic/simulatedNIMBLE.RData"
)
stopCluster(cl)


message("Writing output...")
save(samples,
     file = "Test/multiStateBasic/simulatedNIMBLE2.RData")

# 
# message("Begin sampling")
# samples <- runMCMC(mcmcComp,
#                    niter = n.iter,
#                    nchains = n.chains,
#                    samplesAsCodaMCMC = TRUE)
# 
# 
# 
# 
# reset <- FALSE
# states <- params <- list()
# for(c in 1:n.chains){
#   
#   mcmcComp$run(niter = n.iter, reset = reset)
#   samples <- as.matrix(mcmcComp$mvSamples)
#   x.cols <- grep("x[", colnames(samples), fixed = TRUE)
#   states[[c]] <- samples[,x.cols] %>% as.mcmc()
#   params[[c]] <- samples[,-x.cols] %>% as.mcmc()
# 
# }
# 
# params1 <- as.mcmc.list(params)  
# states1 <- as.mcmc.list(states)
# 
# gelman.diag(params, multivariate = FALSE)
# 
# GBR <- gelman.plot(params)
# burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>1.1,1,any)),1)+1]
# burnin
# 
# 
# 
# 
# 
# ## split output
# out <- list(params = NULL, predict = NULL)
# mfit <- as.matrix(samples, chains = TRUE)
# pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
# chain.col <- which(colnames(mfit) == "CHAIN")
# out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
# out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
# out <- list(out = out)
# 
# 
# 
# 
# 
# 
# message(paste0("Capture history matrix has ", nrow(ch), " rows and ", ncol(ch), " columns"))
# message(paste0("Compiling JAGS with ", n.adapt, " adaptive iterations..."))
# 
# source("R/2.1_multiStateBasic.R")
# compile.start <- Sys.time()
# j.model <- jags.model(file = textConnection(model.code),
#                       data = data,
#                       inits = inits,
#                       n.chains = 3,
#                       n.adapt = n.adapt)
# 
# message("Compiling finished")
# message(paste("Start mcmc sampling:", n.iter, "iterations"))
# 
# jags.out <- coda.samples(j.model,
#                          variable.names = monitor,
#                          n.iter = n.iter)
# message("Sampling complete")
# 
# ## split output
# out <- list(params = NULL, predict = NULL)
# mfit <- as.matrix(jags.out, chains = TRUE)
# pred.cols <- grep("x[", colnames(mfit), fixed = TRUE)
# chain.col <- which(colnames(mfit) == "CHAIN")
# out$predict <- ecoforecastR::mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
# out$params <- ecoforecastR::mat2mcmc.list(mfit[, -pred.cols])
# out <- list(out = out, jags.model = j.model)
# 
# message("Saving samples...")
# save(
#   list =
#     c(
#       "out",
#       "sim",
#       "phi.a",
#       "phi.p",
#       "psi.ap",
#       "psi.pa",
#       "p.a",
#       "p.p",
#       "p.d",
#       "p.u"
#     ),
#   file = file.path(out.dir, "simulatedMCMC.RData")
# )
# message("--- DONE ---")
