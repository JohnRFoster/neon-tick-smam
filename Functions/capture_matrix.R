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
  
  # capture history with plot level identifiers
  ch.by.plot <- df %>% 
    # select(tagID, collectDate, animalInTrap, plotID) %>% 
    mutate(plotNumber = as.numeric(as.factor(plotID)), # numeric plot id
           plotCaptureStatus = plotNumber + animalInTrap - 1, # if seen, plot number
           plotCaptureStatus = if_else(animalInTrap == 0, 0, plotCaptureStatus), # revert not-seen to 0
           plotCaptureStatus = if_else(animalInTrap == 2, -1, plotCaptureStatus), # dead animals to -1
           tagID = if_else(is.na(tagID) & animalInTrap == 0,
                           "noCapture", # dates without captures, need for full matrix
                           tagID)) %>% 
    arrange(collectDate) %>% 
    filter(!is.na(tagID)) # some tags are just NA - remove
  
  # need to figure out what to do about daily capture inconsistencies across plots 
  ch <- ch.by.plot %>% 
    select(tagID, collectDate, plotCaptureStatus) %>% 
    distinct() %>% 
    pivot_wider(names_from = collectDate,
                # values_fn = length,
                values_from = plotCaptureStatus,
                values_fill = 0) %>% 
    filter(tagID != "noCapture")
  
  # make sure mice that are recorded dead stay dead!
  for(i in 1:nrow(ch)){
    ind <- ch[i,]
    if(any(ind == 2)){
      tod <- which(ind == 2)
      if(sum(ind[(tod+1):length(ind)]) > 0){
        stop("Zombie mouse!")
      }
    }
  }
  
  # get animal identifiers    
  smam.ind <- pull(ch, tagID)
  ids <- df %>% 
    filter(tagID %in% smam.ind) %>% 
    select(tagID, genusName, speciesName, nlcdClass, siteID, plotID) %>% 
    distinct() 
  
  # if animals are identified differently after first capture
  # or found in more than one plot
  if(nrow(ids) != nrow(ch)){
    
  }
  
  # might want to do plots by time and individual?
    
  
  
}