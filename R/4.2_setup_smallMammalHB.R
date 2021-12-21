


library(lubridate)
library(nimble)
library(nimbleEcology)
library(tidyverse)
library(LaplacesDemon)

df <- read_csv("Data/allSmallMammals.csv") %>% suppressMessages()
df <- df %>%
  filter(collectYear >= 2016, # subset to training time frame 
         collectYear <= 2020)


site.coord <- readr::read_csv("Data/siteLatLon.csv") %>% suppressMessages()
sites <- site.coord %>% pull(siteID)


message("Building capture history matrix...")
source("Functions/capture_matrix.R")

ch.all <- capture.days.df <- tibble()
site.row <- vector()
p.a.init <- p.p.init <- p.d.init <- p.u.init <- rep(0, length(sites))
for(i in seq_along(sites)){
  site <- sites[i]
  ch <- capture_matrix(site, df) %>% suppressMessages()
  
  alive.u <- ch$alive.u      # observed mouse with unknown tick status
  alive.a <- ch$alive.a      # observed mouse without tick attached
  alive.p <- ch$alive.p      # observed mouse with tick attached
  dead.u <- ch$dead.u       # observed mouse with unknown tick status
  dead.a <- ch$dead.a       # observed mouse without tick attached
  dead.p <- ch$dead.p       # observed mouse with tick attached
  not.seen <- ch$unobserved   # unobserved
  ch <- as_tibble(ch$ch) # the capture matrix
  
  dead <- 4
  unobserved <- dead + 1
  ch[ch == dead.u] <- dead
  ch[ch == dead.a] <- dead
  ch[ch == dead.p] <- dead
  ch[ch == not.seen] <- unobserved
  
  # calculate some inits from ch
  n.captured <- length(ch[ch != unobserved]) # total captured
  
  p.a.init[i] <- length(ch[ch == alive.a]) / n.captured # proportion captured without ticks
  p.p.init[i] <- length(ch[ch == alive.p]) / n.captured # proportion captured with ticks
  p.d.init[i] <- length(ch[ch == dead]) / n.captured # proportion captured dead
  p.u.init[i] <- length(ch[ch == alive.u]) / n.captured # proportion captured unknown ticks
  
  if(!production){
    draw <- sample.int(nrow(ch), 30, replace = FALSE)
    ch <- ch[draw,]
  }
  
  capture.days <- tibble(trap.night = colnames(ch),
                         index = 1,
                         site = site) # trapping occasions are site specific 
  
  ch.all <- bind_rows(ch.all, ch)  # bind CHs, will align by column name (date)
  ch.all[is.na(ch.all)] <- unobserved
  capture.days.df <- bind_rows(capture.days.df, capture.days)
  site.row <- c(site.row, rep(site, nrow(ch))) # vector for site indexing
  

}

# force death probability init to > 0
p.d.init <- pmax(p.d.init, 1e-5)
ch.all[is.na(ch.all)] <- 5
no.captures <- which(apply(ch.all, 1, function(x) all(x == 5)))
if(length(no.captures) > 0){
  ch.all <- ch.all[-no.captures,]
  site.row <- site.row[-no.captures]
}

# make a wide data frame so each site gets a row, and a 1 denotes an observation
# sort rows alphabetically (by site)
trap.occasions <- capture.days.df %>%
  mutate(trap.night = ymd(trap.night)) %>% 
  arrange(trap.night) %>% 
  pivot_wider(names_from = trap.night,
              values_from = index,
              values_fill = 0) %>% 
  arrange(site)

# get site index for capture history matrix sorted alphabetically
site.ch.index <- site.row %>% 
  sort() %>% 
  as.factor() %>% 
  as.numeric()

# get site index for capture occasions
# site.occ.index <- trap.occasions %>%
#   pull(site) %>% 
#   as.factor() %>% 
#   as.numeric()

n.site <- unique(site.ch.index) %>% length()

# no longer need site column
trap.occasions <- trap.occasions %>% select(-site)  
trap.days <- colnames(trap.occasions)

ps.start.index <- apply(trap.occasions, 1, function(x) min(which(x == 1)))
ps.end.index <- apply(trap.occasions, 1, function(x) max(which(x == 1)))

# need to know which day in the entire time frame (all days, not just capture days) that captures happen
# start with figuring out all days; first to last capture across all sites
start <- first(ymd(trap.days))
end <- last(ymd(trap.days))
every.day <- seq.Date(start, end, by = 1) 

# need to know the number of occasions for each site
n.occasions <- apply(trap.occasions, 1, function(x) length(which(x != 0)))

# build matrix that has the index of each trapping occasion for each site
site.occasions <- matrix(NA, n.site, max(n.occasions))
for(i in 1:nrow(trap.occasions)){
  occ <- trap.occasions %>% slice(i)
  occ.index <- which(occ != 0)
  occ.dates <- trap.days[occ.index]
  index <- match(ymd(occ.dates), every.day)
  site.occasions[i, 1:length(index)] <- index
}

capture.index <- which(every.day %in% ymd(trap.days))

ch <- ch.all
ch[is.na(ch)] <- 5
ch <- ch %>% select(all_of(trap.days))
ch <- as.matrix(ch)

# Compute date of first capture
message("Compute date of first capture")
get.first <- function(x) min(which(x != unobserved))
f <- apply(ch, 1, get.first)


if(any(f == ncol(ch))){
  drop.rows <- which(f == ncol(ch))
  ch <- ch[-drop.rows,]
  site.ch.index <- site.ch.index[-drop.rows]
  f <- apply(ch, 1, get.first)
}

if(any(f+1 == ncol(ch))){
  drop.rows <- which(f+1 == ncol(ch))
  ch <- ch[-drop.rows,]
  site.ch.index <- site.ch.index[-drop.rows]
  f <- apply(ch, 1, get.first)
}

# get the last occasion for each site
get.last <- function(x) max(which(x != unobserved))
l <- apply(ch, 1, get.last)

source("Functions/known_state_ms.R")
x.known <- known_state_ms(ch, f, unobserved, alive.u, dead)

x1 <- rep(1, length(f))
for(i in 1:length(f)){
  x1[i] <- x.known[i,f[i]]  
  if(is.na(x1[i])) x1[i] <- rbinom(1, 1, 0.5) + 1
}

Y1 <- c(
  length(x1[x1 == 1]) / length(f),
  length(x1[x1 == 2]) / length(f),
  length(x1[x1 == 3]) / length(f),
  0
)

if(sum(Y1) != 1){
  print(Y1)
  stop("Initial probabilities not equal to 1")
} 

data <- list(
  # x = x.known,
  y = as.matrix(ch),
  Y1 = Y1
)

constants <- list(
  n.ind = nrow(ch),
  n.occasions = ncol(ch),
  # n.occasions = n.occasions,
  site.start = site.occasions[,1],
  site.end = apply(site.occasions, 1, function(x) max(x, na.rm = TRUE)),
  n.days = length(every.day),
  n.site = n.site,
  site.ch.index = site.ch.index,
  no = 5,
  ns = 4,
  f = f,
  l = l,
  capture.index = capture.index,
  # capture.index = site.occasions,
  ps.start.index = ps.start.index,
  ps.end.index = ps.end.index,
  tau.param = 1 / 5^2
)

max.rate <- 0.999
min.rate <- 1e-10

load("/projectnb/dietzelab/fosterj/FinalOut/neon-tick-smam/mice/multiStateDailySurvivalTransitionHB_test/NIMBLE_check.RData")
samples <- as.matrix(save.ls$samples)
mu <- apply(samples, 2, mean)

inits <- function(){list(
  phi.a = jitter(mu["phi.a"]),            
  phi.p = jitter(mu["phi.p"]),
  psi.ap = jitter(mu["psi.ap"]),
  psi.pa = jitter(mu["psi.pa"]),
  phi.p.site = jitter(mu[grep("phi.p.site[", names(mu), fixed = TRUE)]),
  phi.a.site = jitter(mu[grep("phi.a.site[", names(mu), fixed = TRUE)]),
  psi.pa.site = jitter(mu[grep("psi.pa.site[", names(mu), fixed = TRUE)]),
  psi.ap.site = jitter(mu[grep("psi.ap.site[", names(mu), fixed = TRUE)]),
  tau.phi.p.site = jitter(mu["tau.phi.p.site"]),
  tau.phi.a.site = jitter(mu["tau.phi.a.site"]),
  tau.psi.pa.site = jitter(mu["tau.psi.pa.site"]),
  tau.psi.ap.site = jitter(mu["tau.psi.ap.site"]),
  p.a = jitter(mu[grep("p.a[", names(mu), fixed = TRUE)]),
  p.p = jitter(mu[grep("p.p[", names(mu), fixed = TRUE)]),
  p.d = jitter(mu[grep("p.d[", names(mu), fixed = TRUE)]),
  p.u = jitter(mu[grep("p.u[", names(mu), fixed = TRUE)])
)}


message(paste0("Capture history matrix has ", nrow(ch), " rows and ", ncol(ch), " columns"))
source("R/2.3_multiStateHBDailySurvivalTransition.R")

# n.slots <- 1
if(n.slots >= 2){
  
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
  
} else {
  message("Build model")
  
  modelBlock <- nimbleModel(model.code,
                            constants = constants,
                            data = data,
                            inits = inits())
  modelBlock$initializeInfo()
  cModelBlock <- compileNimble(modelBlock)
  mcmcConfBlock <- configureMCMC(cModelBlock,
                                 monitors = monitor)
  
  # mcmcConfBlock$removeSampler(c("psi.ap", "psi.pa"))
  # mcmcConfBlock$addSampler(target = c("psi.ap", "psi.pa"), type = "RW_block")
  # 
  # mcmcConfBlock$removeSampler(c("p.a", "p.p", "p.d", "p.u"))
  # mcmcConfBlock$addSampler(target = c("p.a", "p.p", "p.d", "p.u"), type = "RW_block")
  #
  # mcmcConfBlock$printSamplers(monitor)
  
  mcmcBuildBlock <- buildMCMC(mcmcConfBlock)
  
  compMCMCBlock <- compileNimble(mcmcBuildBlock)
  
  compMCMCBlock$run(niter = 100)
  
  library(coda)
  
  mcmc.eff <- tibble(
    node = names(effectiveSize(as.matrix(compMCMCBlock$mvSamples))),
    effSize_default = effectiveSize(as.matrix(compMCMCBlock$mvSamples)) %>% round(2),
    # effSize_custom = effectiveSize(as.matrix(compMCMC.block$mvSamples)) %>% round(2),
    accRate_default = 1 - rejectionRate(as.mcmc(as.matrix(compMCMCBlock$mvSamples)) %>% round(2))
    # accRate_custom = 1 - rejectionRate(as.mcmc(as.matrix(compMCMC.block$mvSamples)) %>% round(2))
  )
  print(as.matrix(mcmc.eff))
}




message("--- Done ---")




