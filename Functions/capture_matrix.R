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
  
  neon.df <- neon.smam %>% 
    filter(siteID == site)
  
  all.days <- neon.df %>% pull(collectDate) %>% unique()
  
  df <- neon.df %>%
    filter(genusName == "Peromyscus") 
  
  mice.days <- df %>% pull(collectDate) %>% unique()
  
  # the states that we need are:
  alive.u <- 1      # observed mouse with unknown tick status
  alive.a <- 2      # observed mouse without tick attached
  alive.p <- 3      # observed mouse with tick attached
  dead.u <- 4       # observed mouse with unknown tick status
  dead.a <- 5       # observed mouse without tick attached
  dead.p <- 6       # observed mouse with tick attached
  unobserved <- 7   # unobserved
  
  
  # if Peromyscus was not seen on all days, need to fill days back in
  fill.days <- tibble()
  if(length(mice.days) < length(all.days)){
    missing.days <- all.days[which(!(all.days %in% mice.days))]
    fill.days <- neon.smam %>% 
      filter(collectDate %in% missing.days) %>% 
      distinct(collectDate, .keep_all = TRUE) %>% 
      select(collectDate, tagID) %>% 
      mutate(tagID = "noCapture",
             state = unobserved)
  }

  # the total number of unique tags
  total.ind <- df %>% filter(!is.na(tagID)) %>% 
    pull(tagID) %>% 
    unique() %>% 
    length()
  
  # get the state of each mouse at each trap night
  df.state <- df %>% 
    mutate(tickOn = if_else(adultTicksAttached == "Y" |                # any life stage observed
                              nymphalTicksAttached == "Y" |
                              larvalTicksAttached == "Y", 1, 0),
           tickOn = if_else(adultTicksAttached == "N" &                # all life stages not observed
                              nymphalTicksAttached == "N" &
                              larvalTicksAttached == "N", 2, tickOn),
           tickOn = if_else(adultTicksAttached == "U" &                # all life stages unknown
                              nymphalTicksAttached == "U" &
                              larvalTicksAttached == "U", 3, tickOn),
           tickOn = if_else(is.na(tickOn), 3, tickOn),                 # NAs get unknown status
           state = if_else(animalInTrap == 1 & tickOn == 3, alive.u, 0),     # observed animal with unknown tick status
           state = if_else(animalInTrap == 1 & tickOn == 2, alive.a, state), # observed animal without tick attached
           state = if_else(animalInTrap == 1 & tickOn == 1, alive.p, state), # observed animal with tick attached
           state = if_else(animalInTrap == 2 & tickOn == 3, dead.u, state),  # dead animal with unknown tick status
           state = if_else(animalInTrap == 2 & tickOn == 2, dead.a, state),  # dead animal without tick attached
           state = if_else(animalInTrap == 2 & tickOn == 1, dead.p, state),  # dead animal with tick attached
           state = if_else(animalInTrap == 0, unobserved, state),            # unobserved
           tagID = if_else(is.na(tagID) & animalInTrap == 0, 
                           "noCapture", tagID)) %>%     # placeholder for trap nights without any captures           
    select(tagID, collectDate, state) %>% 
    filter(!is.na(tagID)) 
  
  # add missing days
  df.all.days <- bind_rows(df.state, fill.days) %>% 
    mutate(state = if_else(state == 0, unobserved, state))
  
  ch <- df.all.days %>% 
    group_by(collectDate, tagID) %>%
    summarise(state = min(state)) %>%  # conflicting states get unknown designation
    ungroup() %>%
    arrange(collectDate) %>% 
    pivot_wider(names_from = collectDate,
                values_from = state,
                values_fill = unobserved 
                ) %>%
    filter(tagID != "noCapture") %>% 
    select(-tagID)
  
  # check dimensions are what they should be
  if(ncol(ch) != length(all.days)) stop("Conflicting days in data vs capture matrix")
  if(nrow(ch) != total.ind) stop("Conflicting individuals in data vs capture matrix")
  
  # make sure mice that are recorded dead stay dead! 
  mice.found.dead <- which(apply(ch, 1, function(x) any(x %in% c(dead.u, dead.a, dead.p)))) # mice that died
  
  # all days after found dead should be unobserved
  check.dead <- function(x){
    if(any(x %in% c(dead.u, dead.a, dead.p))){
      day.dead <- which(x %in% c(dead.u, dead.a, dead.p))
      if(all(x[(day.dead+1):length(x)] == unobserved)){
        return(1)
      } else {
        return(2)
      }
    }
  }
  
  dead.test <- apply(ch[mice.found.dead,], 1, check.dead)
  if(all(dead.test != 1)) stop("Zombie mouse!")
  
  # change the unobserved state to 0
  # ch[ch == unobserved] <- 0
  
  return(as.matrix(ch)) 
}