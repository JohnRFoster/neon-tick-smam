library(tidyverse)
library(lubridate)

site <- "UKFS"

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
## TO-DO
# add missing dates into met data to get NA cumulative gdd values
# index met with tick observations

# need dates starting with Jan. 1 of the first year through last observation date
days.year <- seq.Date(ymd(paste0(min(year(constants$all.days)), "-01-01")), 
                      max(constants$all.days),
                      by = 1)

days.df <- tibble(Date = days.year)
daily.air.temp <- read_csv("Data/airTempDaily.csv")
daily.air.temp <- daily.air.temp %>% filter(siteID == site)

met.df <- left_join(days.df, daily.air.temp, by = "Date")

met.df <- met.df %>% 
  arrange(Date) %>% 
  group_by(Year) %>%
  mutate(growingDegree = if_else(tempTripleMaximum > 0, tempTripleMaximum, 0),
         cumGDD = cumsum(growingDegree)) %>% 
  filter(Date %in% constants$all.days)


cumGDD <- pull(met.df, cumGDD)
cumGDD.var <- pull(met.df, airTempMaximumVariance)
data$gdd.mu <- cumGDD
data$gdd.mu[is.na(cumGDD)] <- approx(1:length(cumGDD), cumGDD, xout = which(is.na(cumGDD)))$y

data$gdd.tau <- 1 / cumGDD.var
data$gdd.tau[is.na(cumGDD.var)] <- mean(data$gdd.tau, na.rm = TRUE)

constants$gdd.max <- max(cumGDD, na.rm = T) + 10
constants$tau.dem <- 0.04
constants$ns <- 4

data$x.init <- rep(10, 4)

inits <- function() {
  list(
    px = 0.1 + abs(jitter(Y.init[, -1, ])),
    x = 1 + abs(jitter(Y.init)),
    gdd = abs(jitter(data$gdd.mu)),
    repro.mu = jitter(20),
    phi.l.mu = jitter(3),
    phi.n.mu = jitter(4),
    phi.a.mu = jitter(5),
    theta.ln = jitter(0),
    theta.na = jitter(0),
    tau.obs = abs(rnorm(3, 10, 5)),
    sig = abs(rnorm(4, 0, 10))
  )
}

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
