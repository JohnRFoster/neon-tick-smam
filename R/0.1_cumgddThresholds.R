
aa.sites <- c("UKFS", "TALL", "OSBS", "KONZ")
ix.sites <- c("TREE", "HARV")
both.sites <- c("SERC", "SCBI", "ORNL", "LENO", "BLAN")
site.vec <- c(aa.sites, ix.sites, both.sites)
data.a <- read_csv("Data/tickTargets.csv")
start.end <- tibble()
for(s in seq_along(site.vec)){
  site <- site.vec[s]
  source("Functions/tick_data.R")
  tick.ls <- tick_data(site)
  data <- tick.ls$data
  constants <- tick.ls$constants
  
  data <- data.a %>% 
    filter(siteID == site,
           # time >= "2016-01-01",
           time <= "2020-12-31") %>% 
    select(-siteID, -totalCount)
  
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
  
  source("Functions/daymet_functions.R")
  cumgdd <- daymet_cumGDD(site, "tick")
  all.cumGDD <- tibble()
  for(p in 1:constants$n.plots){
    obs.dates <- data.spp %>% 
      filter(plotID == constants$plots[p]) %>% 
      pull(time) %>% 
      unique()
    
    plot.cumGDD.p <- cumgdd %>% 
      ungroup() %>% 
      filter(plotID == constants$plots[p],
             Date %in% obs.dates) 
    
    all.cumGDD <- bind_rows(all.cumGDD, plot.cumGDD.p)
  }
  
  data.spp <- data.spp %>% rename("Date" = time)
  df.gg <- left_join(data.spp, all.cumGDD, by = c("plotID", "Date"))
  
  gg <- df.gg %>% 
    mutate(standardCount = standardCount + 1) %>% 
    ggplot() +
    aes(x = cumGDD, y = standardCount, color = scientificName) +
    geom_point() +
    scale_y_log10() +
    coord_cartesian(xlim = c(0, max(df.gg$cumGDD))) +
    facet_grid(rows = vars(scientificName),
               cols = vars(lifeStage)) +
    labs(title = site) + 
    theme_bw() +
    theme(legend.position = "none")
  print(gg)
  
  start.end.s <- df.gg %>% 
    filter(standardCount > 0) %>% 
    group_by(scientificName, lifeStage) %>% 
    summarise(start = min(cumGDD),
              end = max(cumGDD)) %>% 
    mutate(siteID = site)
  
  if(!(site %in% aa.sites)){
    start.end.ix <- df.gg %>% 
      filter(standardCount > 0) %>% 
      filter(scientificName == "Ixodes scapularis" & lifeStage == "Adult") %>% 
      summarise(start = min(cumGDD[cumGDD > max(cumGDD)*0.75]),
                end = max(cumGDD[cumGDD < max(cumGDD)*0.75])) %>% 
      mutate(scientificName = "Ixodes scapularis",
             lifeStage = "Adult",
             siteID = site)
           
    start.end.s <- start.end.s %>% 
      filter(!(scientificName == "Ixodes scapularis" & lifeStage == "Adult")) %>% 
      bind_rows(start.end.ix)
  }
    
  start.end <- bind_rows(start.end, start.end.s)
}

# double check by visual and see what needs to be changed
start.end <- start.end %>% 
  mutate(end = if_else(siteID == "LENO" & scientificName == "Ixodes scapularis" & lifeStage == "Adult", 4500, end),
         start = if_else(siteID == "LENO" & scientificName == "Ixodes scapularis" & lifeStage == "Adult", 5500, start))
# start.end <- start.end %>% 
#   mutate(start = if_else(siteID == "KONZ" & scientificName == "Amblyomma americanum" & lifeStage == "Larva", 1300, start))
# start.end <- start.end %>% 
#   mutate(start = if_else(siteID == "SCBI" & scientificName == "Amblyomma americanum" & lifeStage == "Larva", 1300, start),
#          start = if_else(siteID == "SCBI" & scientificName == "Ixodes scapularis" & lifeStage == "Larva", 1300, start))

write_csv(start.end, path = "Data/cumgddThresholds.csv")





