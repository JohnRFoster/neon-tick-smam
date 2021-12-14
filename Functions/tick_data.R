library(tidyverse)

tick_data <- function(site){
  data <- read_csv("Data/tickTargets.csv")
  data <- data %>% 
    filter(siteID == site) %>% 
    select(-siteID, -totalCount)
  
  aa.sites <- c("UKFS", "TALL", "OSBS", "KONZ")
  ix.sites <- c("TREE", "HARV")
  both.sites <- c("SERC", "SCBI", "ORNL", "LENO", "BLAN")
  
  if(site %in% both.sites){
    data.aa <- data %>% 
      filter(scientificName %in% c("Amblyomma americanum", "AAorIX")) %>% 
      mutate(scientificName = "Amblyomma americanum")
    data.ix <- data %>% 
      filter(scientificName %in% c("Ixodes scapularis", "AAorIX")) %>% 
      mutate(scientificName = "Ixodes scapularis")
    data.spp <- bind_rows(data.aa, data.ix)
  } else if(site %in% aa.sites) {
    data.spp <- data %>% 
      filter(scientificName != "Ixodes scapularis") %>% 
      mutate(scientificName = "Amblyomma americanum")
  } else if(site %in% ix.sites){
    data.spp <- data %>% 
      filter(scientificName != "Amblyomma americanum") %>% 
      mutate(scientificName = "Ixodes scapularis")
  }
  
  df <- data.spp %>% 
    pivot_wider(names_from = lifeStage,
                values_from = standardCount,
                values_fn = {sum}) %>% 
    arrange(plotID, time)
  
  occasions.plot <- df %>% 
    group_by(plotID, scientificName) %>% 
    summarise(n.occasions = n())
  
  if(!sum(occasions.plot$n.occasions) == nrow(df)){
    stop("Duplicate days in one or more plot")
  }
  
  plots <- df$plotID %>% unique() %>% sort() # vector of plots
  n.plots <- length(plots) # number of plots
  species <- df$scientificName %>% unique() %>% sort()
  n.species <- length(species) 
  
  # make the array of tick densities by plot and species
  Y <- array(NA, dim = c(3, # life stage
                         max(occasions.plot$n.occasions), # number of drag occasions
                         n.plots, # number of plots
                         n.species)) # number of species
  
  # matrix of drag occasion dates for each plot
  plot.dates <- matrix(NA, n.plots, max(occasions.plot$n.occasions))
  diff.time <- matrix(NA, n.plots, max(occasions.plot$n.occasions)-1)
  
  # nlcd class and total number of days for each plot
  nlcd <- N <- rep(NA, n.plots)
  
  for(spp in 1:n.species){
    for(p in seq_len(n.plots)){
      tick.plot <- df %>% 
        filter(plotID == plots[p],
               scientificName == species[spp])
      n.days <- nrow(tick.plot)
      Y[1, 1:n.days, p, spp] <- pull(tick.plot, Larva)
      Y[2, 1:n.days, p, spp] <- pull(tick.plot, Nymph)
      Y[3, 1:n.days, p, spp] <- pull(tick.plot, Adult)
      
      drags <- pull(tick.plot, time)
      plot.dates[p, 1:n.days] <- as.character(drags)
      diff.time[p, 1:(n.days-1)] <- diff.Date(drags) %>% as.numeric()
      N[p] <- sum(diff.Date(drags)) %>% as.numeric()
      nlcd[p] <- tick.plot$nlcdClass[1]
      
    }  
  }

  seq.days <- array(NA, dim = c(n.plots, max(occasions.plot$n.occasions)-1, max(diff.time, na.rm = TRUE)))
  for(p in seq_len(n.plots)){
    dt.index <- c(1, cumsum(diff.time[p,]))
    dt.index <- dt.index[!is.na(dt.index)]
    n.days <- length(which(!is.na(diff.time[p,])))
    max.interval <- max(diff.time[p,], na.rm = TRUE)
    for (i in 1:n.days) {
      xx <- (dt.index[i+1]-1):dt.index[i]
      seq.days[p, i, 1:length(xx)] <- xx
    }
  }
  
  return(list(y = Y,
              plots = plots,
              n.plots = n.plots,
              species = species,
              n.species = n.species,
              plot.dates = plot.dates,
              diff.time = diff.time,
              n.occ.plot = pull(occasions.plot, n.occasions),
              N = max(N),
              seq.days = seq.days,
              occasions.plot = occasions.plot,
              nlcd = nlcd))
}
