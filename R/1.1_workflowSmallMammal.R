# ========================================================== #
#       Workflow script for multi-state mouse models         #
# ========================================================== #

library(tidyverse)
library(parallel)
library(nimble)
library(nimbleEcology)

# script flow
n.slots <- 4
production <- TRUE
HB <- TRUE
weather <- TRUE
run.intake <- FALSE
check.neon <- FALSE

out.dir <- "/projectnb/dietzelab/fosterj/FinalOut/neon-tick-smam/mice"
model.dir <- "multiStateHBAllSiteEffects"


# sites to run
site.coord <- readr::read_csv("Data/siteLatLon.csv") %>% suppressMessages()
sites <- site.coord %>% pull(siteID)

# done.sites <- c("LENO", "OSBS", "TALL", "TEAK", "DELA")
# sites <- sites[!(sites %in% done.sites)]
# sites <- c("HARV", "SOAP")

# jags arguments
if(production) {
  n.burnin <- 0
  thin <- 1
  n.iter <- if_else(HB, 1000, 50000)
  max.iter <- 7e6
  
} else { # testing / dev
  n.burnin <- 0
  thin <- 1
  n.iter <- 100
  max.iter <- 1000000
  ind.test <- 50
  n.occ.test <- 100
  model.dir <- paste0(model.dir, "_test")
}

# do we need to run the intake script?
if(run.intake){
  message(paste0("Running mouse intake at ", Sys.time()))
  source("R/0_intakeSmallMammal.R")
  mouse_intake(check.neon)
  message(paste0("Mouse intake complete at ", Sys.time()))
}

# run the model individually at each site
if(HB){
  # load model
  source("R/2.3_multiStateHBDailySurvivalTransition.R")
  site.dir <- file.path(out.dir, model.dir)
} else {
  # load model
  source("R/2.2_multiStateDailySurvivalTransitionNimble.R")
  
  # array.num <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number
  # site <- sites[array.num]
  site <- "HARV"
  site.dir <- file.path(out.dir, model.dir, site)
} 

# create file path for output
if(!dir.exists(site.dir)) dir.create(site.dir, recursive = TRUE)
  
start.time <- Sys.time()
message(paste("Number of slots (chains):", n.slots))
  
if(HB){
  message(paste("Fitting hierarchical mouse model:", model.dir))
  source("R/4.2_setup_smallMammalHB.R")
} else {
  message(paste("Fitting", model.dir, "at", site))
  source("R/4.1_setup_smallMammal.R")
}
  
end.time <- Sys.time()
message("Total run time:")
print(end.time-start.time)


