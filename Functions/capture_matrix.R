#' this function reads the small mammal csv and creates a 
#' capture history matrix for given scale
#' 0 = not captured
#' 1 = alive
#' 2 = dead
#' 
#' @param site the unit to extract; plot ("HARV_016"), site ("HARV"), domain ("D02")
#' @param neon.smam the small mammal neon.smam frame from the csv in /neon.smam
#' @param by_bout should capture history be for each collectDate (the default) or by trapping bout?


capture_matrix <- function(site, neon.smam, by_bout = FALSE){
  
  neon.smam <- neon.smam %>% 
    filter(siteID == site)
  
  all.days <- neon.smam %>% pull(collectDate) %>% unique()
  
  df <- neon.smam %>%
    filter(genusName == "Peromyscus") 
  
  mice.days <- df %>% pull(collectDate) %>% unique()
  
  # if Peromyscus was not seen on all days, need to fill days back in
  fill.days <- tibble()
  if(length(mice.days) < length(all.days)){
    missing.days <- all.days[which(!(all.days %in% mice.days))]
    fill.days <- neon.smam %>% 
      filter(collectDate %in% missing.days) %>% 
      distinct(collectDate, .keep_all = TRUE) %>% 
      select(collectDate, tagID) %>% 
      mutate(tagID = "noCapture",
             state = 4)
  }
  
  # the states that we need are:
  # 1 = dead
  # 2 = unobserved
  # 3 = observed mouse with tick attached
  # 4 = observed mouse without tick attached
  # 5 = observed mouse with unknown tick status

  # the total number of unique tags
  total.ind <- df %>% filter(!is.na(tagID)) %>% 
    pull(tagID) %>% 
    unique() %>% 
    length()
  
  # get the state of each mouse at each trap night
  df.state <- df %>% 
    mutate(tickOn = if_else(adultTicksAttached == "Y" |                # any life stage observed
                              nymphalTicksAttached == "Y" |
                              larvalTicksAttached == "Y", 3, 0),
           tickOn = if_else(adultTicksAttached == "N" &                # all life stages not observed
                              nymphalTicksAttached == "N" &
                              larvalTicksAttached == "N", 4, tickOn),
           tickOn = if_else(adultTicksAttached == "U" &                # all life stages unknown
                              nymphalTicksAttached == "U" &
                              larvalTicksAttached == "U", 5, tickOn),
           tickOn = if_else(is.na(tickOn), 3, tickOn),                 # NAs get unknown status
           state = if_else(animalInTrap == 1 & tickOn == 1, 3, 0),     # observed animal with tick attached
           state = if_else(animalInTrap == 1 & tickOn == 0, 4, state), # observed animal without tick attached
           state = if_else(animalInTrap == 1 & tickOn == 3, 5, state), # observed animal with unknown tick status
           state = if_else(animalInTrap == 0, 2, state),               # unobserved
           state = if_else(animalInTrap == 2, 1, state),               # dead 
           tagID = if_else(is.na(tagID) & animalInTrap == 0, 
                           "noCapture", tagID)) %>%     # placeholder for trap nights without any captures           
    select(tagID, collectDate, state) %>% 
    filter(!is.na(tagID)) 
  
  # add missing days
  df.all.days <- bind_rows(df.state, fill.days) 
  
  ch <- df.all.days %>% 
    arrange(collectDate) %>% 
    # group_by(collectDate, tagID) %>%
    # distinct() %>% 
    # ungroup() %>%
    pivot_wider(names_from = collectDate,
                values_from = state,
                values_fill = 4, # unobserved state
                values_fn = {max} # conflicting states get unknown designation 
                ) %>%
    filter(tagID != "noCapture")
  
  # check dimensions are what they should be
  if(ncol(ch[,-1]) != length(all.days)) stop("Conflicting days in data vs capture matrix")
  if(nrow(ch) != total.ind) stop("Conflicting individuals in data vs capture matrix")
  
  # make sure mice that are recorded dead stay dead! 
  mice.found.dead <- which(apply(ch[,-1], 1, function(x) any(x == 1))) # mice that died
  
  # all days after found dead should have state of 4
  check.dead <- function(x){
    day.dead <- which(x == 1)
    if(all(x[(day.dead+1):length(x)] == 4)){
      return(1)
    } else {
      return(0)
    }
  }
  
  dead.test <- apply(ch[mice.found.dead,-1], 1, check.dead)
  if(all(dead.test != 1)) stop("Zombie mouse!")
  
  return(ch)
}