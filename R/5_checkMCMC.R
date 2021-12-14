# ========================================================== #
#       script for checking mcmc chains on the cluster       #
# ========================================================== #

library(rjags)
library(runjags)
library(tidyverse)

override <- TRUE
max.segments <- 100

org <- "mice"
model.dir <- "multiStateBasic"

# will be able to check all sites at once
array.num <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number
# array.num <- 5

site.coord <- readr::read_csv("Data/siteLatLon.csv") %>% suppressMessages()
sites <- site.coord %>% pull(siteID)

site <- sites[array.num]
# site <- "LENO"
message("Checking MCMC...")
message(paste("Organism:", org))
message(paste("Model:", model.dir))
message(paste("Site:", site))
message(paste("Override:", override))
message(paste("Max segments to combine:", max.segments))

top.dir <- "/projectnb/dietzelab/fosterj/FinalOut/neon-tick-smam"
mcmc.dir <- file.path(top.dir, org, model.dir, site)

source("Functions/check_mcmc.R")
check_mcmc(mcmc.dir, override = override, max.segments = max.segments)

