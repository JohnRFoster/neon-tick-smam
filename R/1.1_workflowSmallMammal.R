# ========================================================== #
#       Workflow script for multi-state mouse models         #
# ========================================================== #

library(tidyverse)

# load model
source("R/2.1_multiStateBasic.R")

out.dir <- "/projectnb/dietzelab/fosterj/FinalOut/neon-tick-smam/mice"
model.dir <- "multiStateBasic"

# script flow
production <- TRUE
run.intake <- FALSE
check.neon <- FALSE
run.parallel <- FALSE
array.job <- FALSE
n.chains <- 4

# sites to run
site.coord <- readr::read_csv("Data/siteLatLon.csv")
sites <- site.coord %>% pull(siteID)

# jags arguments
if(production) {
  n.adapt <- 2500
  thin <- 1
  n.iter <- 100000
  n.loops <- 1000
  
  jobs <- tibble(site = rep(sites, each = n.chains),
                 chain = rep(1:n.chains, length(sites)))
  
  # run the model individually at each site
  array.num <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number
  site <- jobs %>% slice(array.num) %>% pull(site)
  chain <- jobs %>% slice(array.num) %>% pull(chain)
  
} else { # testing / dev
  n.adapt <- 50
  thin <- 1
  n.iter <- 50
  n.loops <- 10
  site <- sites[sample(length(sites), 1)]
  chain <- 1
}


# do we need to run the intake script?
if(run.intake){
  message(paste0("Running mouse intake at ", Sys.time()))
  source("R/0_intakeSmallMammal.R")
  mouse_intake(check.neon)
  message(paste0("Mouse intake complete at ", Sys.time()))
}

# create file path for output
site.dir <- file.path(out.dir, model.dir, site)
if(!dir.exists(site.dir)) dir.create(site.dir, recursive = TRUE)
  
start.time <- Sys.time()
message(paste("Fitting mouse model at", site, start.time))
message(paste0("Chain: ", chain))
  
source("R/4.1_setup_smallMammal.R")
  
end.time <- Sys.time()
message("Total run time:")
print(end.time-start.time)


