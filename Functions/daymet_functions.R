# =========================================================================== #
# functions for extracting and working with daymet
# daymet data has already been downloaded with R/0_intakeDayMet.R
# =========================================================================== #


library(tidyverse)

#' function that calculates cumulative growing degree days for each plot
#' @param site the site being modeled
#' @param org either "tick" or "smam"
daymet_cumGDD <- function(site, org){
  df.all <- read_csv("Data/daymet_maxTemperature.csv")
  df <- df.all %>% 
    filter(grepl(site, plotID),
           data == org) %>% 
    group_by(plotID, year) %>% 
    mutate(growingDegree = if_else(maxTemperature > 10, maxTemperature - 10, 0),
           cumGDD = cumsum(growingDegree)) %>% 
    select(Date, plotID, cumGDD, year)
  
  return(df)
}

