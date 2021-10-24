# mouse workflow 

library(tidyverse)

# load model
source("R/2.1_multiStateBasic.R")

out.dir <- "/projectnb/dietzelab/fosterj/FinalOut/neon-tick-smam/mice"
model.dir <- "multiStateBasic"

# script flow
production <- FALSE
run.intake <- FALSE
check.neon <- FALSE
run.parallel <- FALSE
array.job <- FALSE

# sites to run
site.coord <- readr::read_csv("Data/siteLatLon.csv")
sites <- site.coord %>% pull(siteID)

# jags arguments
if(production) {
  n.adapt <- 5000
  n.chains <- 1
  thin <- 1
  n.iter <- 100000
  n.loops <- 1000
} else {
  n.adapt <- 50
  n.chains <- 1
  thin <- 1
  n.iter <- 50
  n.loops <- 10
  sites <- sites[sample(length(sites), 1)]
}


# do we need to run the intake script?
if(run.intake){
  message(paste0("Running mouse intake at ", Sys.time()))
  source("R/0_intakeSmallMammal.R")
  mouse_intake(check.neon)
  message(paste0("Mouse intake complete at ", Sys.time()))
}

# run the model individually at each site
for(s in seq_along(sites)){
  site <- sites[s]
  
  site.dir <- file.path(out.dir, model.dir, site)
  
  # create file path for output
  if(!dir.exists(site.dir)) dir.create(site.dir, recursive = TRUE)
  
  start.time <- Sys.time()
  message(paste("Running mouse model at", site, start.time))
  
  source("R/4.1_setup_smallMammal.R")
  
  end.time <- Sys.time()
  message(paste("  ", site, "completed at", end.time))
  message(paste("   Total run time:", end.time-start.time))
}


