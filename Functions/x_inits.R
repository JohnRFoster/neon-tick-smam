#' function to supply latent state initial condition
#' anything coded as unobserved or unknown tick status get NA
#' anything after first capture and NA in x.known gets a guess
#' 
#' @param ks known states capture history matrix
#' @param f vector of time of first capture

x_inits <- function(ks, f){
  states <- ks
  v <- which(is.na(states))
  states[-v] <- NA
  states[v] <- rbinom(length(v), 1, 0.5) + 1
  for (i in 1:nrow(ks)){
    if(f[i] > 1) states[i, 1:(f[i]-1)] <- NA
  }
  return(states)
}