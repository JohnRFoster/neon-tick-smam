# ========================================================== #
#       Workflow script for tick population models           #
# ========================================================== #

out.dir <- "/projectnb/dietzelab/fosterj/FinalOut/neon-tick-smam/tick"
model.dir <- "static"

library(tidyverse)
library(parallel)
library(nimble)

# script flow
n.slots <- Sys.getenv("NSLOTS")
production <- FALSE
check.model <- TRUE
HB <- FALSE
weather <- FALSE
run.intake <- FALSE
check.neon <- FALSE

aa.sites <- c("UKFS", "TALL", "OSBS", "KONZ")
ix.sites <- c("TREE", "HARV")
both.sites <- c("SERC", "SCBI", "ORNL", "LENO", "BLAN")

ix.jobs <- tibble(species = "IX",
                  site = c(ix.sites, both.sites))
aa.jobs <- tibble(species = "AA",
                  site = c(aa.sites, both.sites))

all.jobs <- bind_rows(ix.jobs, aa.jobs)

array.num <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number
if(is.na(array.num)) array.num <- 2

job <- all.jobs %>% slice(array.num)
species <- job %>% pull(species)
site <- job %>% pull(site)

if(species == "IX"){
  source("R/3.1.1_tickStaticIX.R")
} else {
  source("R/3.1.2_tickStaticAA.R")
} 

site.dir <- file.path(out.dir, model.dir, species, site)

if(production){
  n.burnin <- 0
  n.iter <- 20000
  max.iter <- 7e6
} else {
  n.burnin <- 0
  n.iter <- 5000
  max.iter <- 2e6
  site.dir <- paste0(site.dir, "_test")
}

# create file path for output
if(!dir.exists(site.dir)) dir.create(site.dir, recursive = TRUE)


message(paste("Model:", model.dir))
message(paste("Slots:", n.slots))
message(paste("Production:", production))
message(paste("HB:", HB))
message(paste("Species:", species))
message(paste("Site:", site))
message(paste("Writing output to:", site.dir))


system.time(
  source("R/4.3_setup_ticks.R")
)

message("--- DONE ---")


