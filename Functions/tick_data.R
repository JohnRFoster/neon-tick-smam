library(tidyverse)

tick_data <- function(site, spp.model){
  data <- read_csv("Data/tickTargets.csv")
  data <- data %>% 
    filter(siteID == site,
           time <= "2020-12-31") %>% 
    select(-siteID, -occasionID)
  
  aa.sites <- c("UKFS", "TALL", "OSBS", "KONZ")
  ix.sites <- c("TREE", "HARV")
  both.sites <- c("SERC", "SCBI", "ORNL", "LENO", "BLAN")
  
  if(site == "ORNL"){
    data <- data %>% 
      filter(plotID != "ORNL_006")
  }
  
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
  
  spp.filter <- if_else(spp.model == "IX", "Ixodes scapularis", "Amblyomma americanum")
  
  df <- data.spp %>% 
    filter(scientificName == spp.filter) %>% 
    pivot_wider(names_from = lifeStage,
                values_from = processedCount,
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
  Y <- array(NA, dim = c(4, # life stage
                         max(occasions.plot$n.occasions), # number of drag occasions
                         n.plots, # number of plots
                         n.species)) # number of species
  Y.init <- Y
  
  # matrix of drag occasion dates for each plot
  plot.dates <- area.sampled <- matrix(NA, n.plots, max(occasions.plot$n.occasions))
  diff.time <- matrix(NA, n.plots, max(occasions.plot$n.occasions)-1)
  
  # nlcd class and total number of days for each plot
  nlcd <- N <- start.date <- end.date <- rep(NA, n.plots)
  
  for(spp in 1:n.species){
    for(p in seq_len(n.plots)){
      tick.plot <- df %>% 
        filter(plotID == plots[p],
               scientificName == species[spp])
      
      n.days <- nrow(tick.plot)
      Y[1, 1:n.days, p, spp] <- pull(tick.plot, Larva)
      Y[3, 1:n.days, p, spp] <- pull(tick.plot, Nymph)
      Y[4, 1:n.days, p, spp] <- pull(tick.plot, Adult)
      Y.init[1, 1:n.days, p, spp] <- Y[1, 1:n.days, p, spp]
      Y.init[3, 1:n.days, p, spp] <- Y[3, 1:n.days, p, spp]
      Y.init[4, 1:n.days, p, spp] <- Y[4, 1:n.days, p, spp]
      for(t in 1:n.days){
        Y.init[2, t, p, spp] <- mean(c(Y[1, t, p, spp], Y[3, t, p, spp])  )
      }
      
      drags <- pull(tick.plot, time) # when drags occurred
      start.date[p] <- first(drags) %>% as.character() # first drag
      end.date[p] <- last(drags) %>% as.character() # last drag
      plot.dates[p, 1:n.days] <- as.character(drags) # matrix of all drag days for each plot
      area.sampled[p, 1:n.days] <- pull(tick.plot, totalSampledArea) # matrix sampling effort
      diff.time[p, 1:(n.days-1)] <- diff.Date(drags) %>% as.numeric() # number of days between drag events
      N[p] <- as.numeric(ymd(end.date[p]) - ymd(start.date[p]) + 1) # total number of days in time series
      nlcd[p] <- tick.plot$nlcdClass[1] # nlcd class for each plot
      
    }  
  }
  
  first.day <- min(ymd(start.date)) # first day at site level
  last.day <- max(ymd(end.date)) # last day at site level
  all.days <- seq.Date(first.day, last.day, by = 1) # sequence of all days in timeseries at site level
  # all.days <- all.days[-which(all.days == "2016-12-31")]
  plot.start <- match(ymd(start.date), all.days) # first day of each plot w.r.t. all.days
  plot.end <- match(ymd(end.date), all.days) # last day of each plot w.r.t. all.days
  
  seq.days <- array(NA, dim = c(n.plots, max(occasions.plot$n.occasions)-1, max(diff.time, na.rm = TRUE)+1))
  p.index <- matrix(NA, n.plots, max(occasions.plot$n.occasions)-1)
  for(p in seq_len(n.plots)){
    dt.index <- cumsum(c(1, diff.time[p,]))
    dt.index <- dt.index[!is.na(dt.index)] # days index
    dt.index <- dt.index + plot.start[p] - 1 # match to site level sequence
    n.days <- length(which(!is.na(diff.time[p,])))
    max.interval <- max(diff.time[p,], na.rm = TRUE)
    for (i in 1:n.days) {
      xx <- (dt.index[i+1]-1):dt.index[i]
      seq.days[p, i, 1:length(xx)] <- xx
      p.index[p,i] <- min(xx)
    }
  }
  
  data <- list(y = Y,
               area.sampled = area.sampled)
  constants <- list(plots = plots,
                    n.plots = n.plots,
                    species = species,
                    n.species = n.species,
                    plot.dates = plot.dates,
                    diff.time = diff.time,
                    all.days = all.days,
                    plot.start = plot.start,
                    plot.end = plot.end,
                    n.occ.plot = pull(occasions.plot, n.occasions),
                    N = max(N) - 1,
                    p.index = p.index,
                    seq.days = seq.days,
                    # occasions.plot = occasions.plot,
                    nlcd = nlcd,
                    Y.init = Y.init)
  
  return(list(data = data,
              constants = constants))
}
