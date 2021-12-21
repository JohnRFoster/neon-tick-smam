library(tidyverse)
library(lubridate)

site <- "TREE"

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
 
daily.ground.temp <- read_csv("Data/bioTempDaily.csv")
met.df <- daily.ground.temp %>% 
  filter(siteID == site) %>% 
  arrange(Date) %>% 
  group_by(Year) %>% 
  mutate(growingDegree = if_else(bioTempMaximum > 0, bioTempMaximum, 0)) %>% 
  mutate(cumGDD = cumsum(growingDegree))

gdd <- rep(seq(0, 3000, length.out = 365), 6)
constants$gdd <- gdd[1:constants$N]
constants$tau.dem <- 0.04
constants$ns <- 4

data$x.init <- rep(10, 4)

#### TO DO
# TRY WITH A INDEXED BY PLOT OR PLOT SPECIFIC INDICIES TO MATCH DAYS WHEN PERMUTING A MATRIX
      # WILL NEED THE A LATTER WHEN ADDING WEATHER
  
inits <- function(){list(px = 0.1 + jitter(Y.init[,-1,]),
                         x = 1 + round(jitter(Y.init)),
                            repro.mu = jitter(20),
                            phi.l.mu = jitter(0),
                            phi.n.mu = jitter(0),
                            phi.a.mu = jitter(0),
                            theta.ln = jitter(0),
                            theta.na = jitter(0),
                            tau.obs.l = jitter(2),
                            tau.obs.n = jitter(2),
                            tau.obs.a = jitter(2),
                            tau.l = jitter(2),
                            tau.d = jitter(2),
                            tau.n = jitter(2),
                            tau.a = jitter(2))}

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
mcmcBuild <- buildMCMC(mcmcConf)
compMCMC <- compileNimble(mcmcBuild)
compMCMC$run(niter = 1000)

library(coda)

mcmc.eff <- tibble(
  node = names(effectiveSize(as.matrix(compMCMC$mvSamples))),
  effSize_default = effectiveSize(as.matrix(compMCMC$mvSamples)) %>% round(2),
  accRate_default = 1 - rejectionRate(as.mcmc(as.matrix(compMCMC$mvSamples)) %>% round(2))
)
print(as.matrix(mcmc.eff))








# 
# 
# 
# 
# 
# n.days <- df$time %>% unique() %>% length()
# n.days == nrow(df)
# 
# 
# 
# 
# 
# 
# #### maybe just use plots as an effect and track a vector of nlcdclass types?
# 
# sites <- target$siteID %>% unique()
# 
# for(i in seq_along(sites)){
#   gg <- target %>% 
#     filter(siteID == sites[i]) %>% 
#     mutate(standardCount = standardCount + 1) %>% 
#     ggplot() +
#     aes(x = time, y = standardCount) +
#     geom_point(aes(color = scientificName, shape = plotID)) +
#     scale_y_log10() +
#     facet_grid(rows = vars(lifeStage),
#                cols = vars(nlcdClass)) +
#     labs(title = sites[i]) +
#     theme_bw() +
#     theme(axis.text = element_text(size = 18),
#           title = element_text(size = 20),
#           axis.title.x = element_blank(),
#           axis.title.y = element_text(size = 16))
#   print(gg)
# }
