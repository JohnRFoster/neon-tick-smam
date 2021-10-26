# ========================================================== #
#       script for checking mcmc chains on the cluster       #
# ========================================================== #

library(rjags)
library(runjags)
library(tidyverse)


# will be able to check all sites at once
array.num <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number

site.coord <- readr::read_csv("Data/siteLatLon.csv")
sites <- site.coord %>% pull(siteID)

site <- sites[array.num]

top.dir <- "/projectnb/dietzelab/fosterj/FinalOut/neon-tick-smam/mice/multiStateBasic"
mcmc.dir <- file.path(top.dir, site)

source("Functions/check_mcmc.R")
check_mcmc(mcmc.dir)

