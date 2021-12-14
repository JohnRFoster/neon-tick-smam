# ====================================================== #
#              NEON driver data download                 #
# ====================================================== #

library(neonstore)
library(tidyverse)
library(lubridate)

# products we want
driver.products <- c(
  "DP1.00098.001", # relative humidity 
  "DP1.00003.001", # triple aspirated air temperature
  "DP1.00005.001", # IR biological temperature (ground temp)
  "DP1.00094.001", # soil water content and water salinity
  "DP1.00041.001"  # soil temperature
)
array.num <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number
product <- driver.products[array.num]

message("Updating NEON Data")
print(product)
neon_download(product = product)
message("---- DONE ----")
