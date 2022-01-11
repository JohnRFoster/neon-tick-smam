library(tidyverse)
library(lubridate)

site <- "TREE"
spp <- "Ixodes scapularis"
source("Functions/tick_data.R")
tick.ls <- tick_data(site)
data <- tick.ls$data
constants <- tick.ls$constants

if(constants$n.species == 1){
  data$y <- array(tick.ls$data$y[,,,1], dim = dim(tick.ls$data$y)[-4])
  Y.init <- array(tick.ls$constants$Y.init[,,,1], dim = dim(tick.ls$data$y)[-4]) 
} else {
  Y.init <- constants$Y.init
}

# Y.init[2,,] <- Y.init[2,,] + 500
Y.init <- Y.init + 10

# get cumGDD
source("Functions/daymet_functions.R")
cumgdd <- daymet_cumGDD(site, "tick")
all.cumGDD <- tibble()
for(p in 1:constants$n.plots){
  plot.cumGDD.p <- cumgdd %>% 
    ungroup() %>% 
    filter(plotID == constants$plots[p],
           Date %in% constants$all.days) %>%
    select(-year) %>% 
    pivot_wider(names_from = "Date",
                values_from = "cumGDD")
  
  all.cumGDD <- bind_rows(all.cumGDD, plot.cumGDD.p)
}

data$cumGDD.mu <- all.cumGDD %>% 
  select(-plotID) %>% 
  as.matrix()

# transition thresholds
trans.cgdd <- read_csv("Data/cumgddThresholds.csv") 
trans.cgdd <- trans.cgdd %>% 
  filter(siteID == site,
         scientificName == spp) 
data$rho.l1 <- trans.cgdd %>% filter(lifeStage == "Larva") %>% pull(start)
data$rho.l2 <- trans.cgdd %>% filter(lifeStage == "Larva") %>% pull(end)
data$rho.n1 <- trans.cgdd %>% filter(lifeStage == "Nymph") %>% pull(start)
data$rho.n2 <- trans.cgdd %>% filter(lifeStage == "Nymph") %>% pull(end)
data$rho.a1 <- trans.cgdd %>% filter(lifeStage == "Adult") %>% pull(start)
data$rho.a2 <- trans.cgdd %>% filter(lifeStage == "Adult") %>% pull(end)

constants$tau.dem <- 0.04
constants$ns <- 4

data$x.init <- rep(10, 4)

inits <- function() {
  list(
    px = abs(jitter(Y.init[, -1, ])),
    x = abs(jitter(Y.init)),
    cumgdd = abs(jitter(data$cumGDD.mu)),
    repro.mu = jitter(50),
    phi.l.mu = jitter(3),
    phi.n.mu = jitter(4),
    phi.a.mu = jitter(5),
    theta.ln = jitter(-6),
    theta.na = jitter(-6),
    tau.obs = abs(rnorm(3, 10, 5)),
    tau.gdd = abs(rnorm(1, 5, 2)),
    sig = abs(rnorm(4, c(1e+5, 1000, 1, 20), 3))
  )
}
itest <- inits()

message("Build model")
source("R/3.1_tickStaticOneSpecies.R")
model <- nimbleModel(model.code,
                     constants = constants,
                     data = data,
                     inits = inits())

model$initializeInfo()
cModel <- compileNimble(model)
monitor.check <- monitor[monitor != "x"]
mcmcConf <- configureMCMC(cModel, monitors = monitor.check)

# mcmcConf$removeSampler(target = "phi.a.mu")
# mcmcConf$addSampler(target = "phi.a.mu", type = "slice")
# mcmcConf$removeSampler(target = "phi.l.mu")
# mcmcConf$addSampler(target = "phi.l.mu", type = "slice")
# mcmcConf$removeSampler(target = "phi.n.mu")
# mcmcConf$addSampler(target = "phi.n.mu", type = "slice")

mcmcBuild <- buildMCMC(mcmcConf)
compMCMC <- compileNimble(mcmcBuild)
compMCMC$run(niter = 2500)

library(coda)

mcmc.eff <- tibble(
  node = names(effectiveSize(as.matrix(compMCMC$mvSamples))),
  effSize_default = effectiveSize(as.matrix(compMCMC$mvSamples)) %>% round(2),
  accRate_default = 1 - rejectionRate(as.mcmc(as.matrix(compMCMC$mvSamples)) %>% round(2))
)

print(as.matrix(mcmc.eff))
print(head(as.matrix(compMCMC$mvSamples)))

