#' this function reads the small mammal csv and creates a 
#' capture history matrix for given scale
#' 0 = not captured
#' 1 = alive
#' 2 = dead
#' 
#' @param unit the unit to extract; plot ("HARV_016"), site ("HARV"), domain ("D02")
#' @param neon.smam the small mammal neon.smam frame from the csv in /neon.smam
#' @param by_bout should capture history be for each collectDate (the default) or by trapping bout?


capture_matrix <- function(unit, neon.smam, by_bout = FALSE){
  
  if(!exists("neon.smam")) {
    neon.smam <- read_csv("Data/allSmallMammals.csv")
  }
  
  if(grepl("_[[:digit:]]{3}", unit)){ # plot
    scale <- "plot"
    df <- neon.smam %>% 
      filter(plotID == unit)
  } else if(grepl("[[:alpha:]]{4}", unit)){ # site
    scale <- "site"
    df <- neon.smam %>% 
      filter(siteID == unit)
  } else if(grepl("D[[:digit:]]{2}", unit)){ # domain
    scale <- "domin"
    df <- neon.smam %>% 
      filter(domainID == unit)
  } else {
    stop("unit argument not recognized")
  }
  message("Scale = ", scale)
  
  
  # control flow for by bout
  ch <- df %>% 
    select(tagID, collectDate, animalInTrap) %>% 
    pivot_wider(names_from = collectDate,
                values_from = animalInTrap)
    
  
}