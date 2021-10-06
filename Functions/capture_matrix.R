#' this function reads the small mammal csv and creates a 
#' capture history matrix for given scale
#' 0 = not captured
#' 1 = alive
#' 2 = dead
#' 
#' @param unit the unit to extract; plot, site, domain
#' @param data the small mammal csv in /Data
#' @param by_bout should capture history be for each collectDate (the default) or by trapping bout?


capture_matrix <- function(unit, data, by_bout = FALSE){
  
  if(grepl("_[[:digit:]]{3}", unit)){ # plot
    scale <- "plot"
    df <- data %>% 
      filter(plotID == unit)
  } else if(grepl("[[:alpha:]]{4}", unit)){ # site
    scale <- "site"
    df <- data %>% 
      filter(siteID == unit)
  } else if(grepl("D[[:digit:]]{2}", unit)){ # domain
    scale <- "domin"
    df <- data %>% 
      filter(domainID == unit)
  } else {
    stop("unit argument not recognized")
  }
  
  # control flow for by bout
  
  ch <- df %>% 
    select(tagID, collectDate)
    
  
}