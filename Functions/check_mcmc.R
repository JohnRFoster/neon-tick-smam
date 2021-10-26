#' function for combining individual chains into one mcmc object
#' then check convergence diagnostics and save output
#' 
#' @param mcmc.dir file path to where mcmc segments are stored
#' @param n.mcmc the number of posterior iterations to save
#' @param min.effect.size the minimum effective sample size allowed
#' @param override override diagnostics so that output is saved - for diagnosing non convergence

library(ecoforecastR)
library(runjags)

check_mcmc <- function(mcmc.dir, n.mcmc = 10000, gbr.threshold = 1.1, 
                       min.effect.size = 3000, override = FALSE){
  
  # list mcmc output files 
  mcmc.files <- list.files(mcmc.dir)
  
  # check if already completed
  finished <- grepl("mcmc", mcmc.files)
  if(any(finished) & !override) stop("Already combined!", call. = FALSE)
  
  # remove any other files from list and figure out the number of chains
  mcmc.files <- mcmc.files[grep("[[:alpha:]]{4}_[[:digit:]]_[[:digit:]]", mcmc.files)]
  chains <- str_extract(mcmc.files, "[0-9]+") %>% 
    unique() %>% 
    as.integer()
  
  # count number of segments per chain
  split.chains <- rep(NA, length(chains))
  for(c in chains){
    pattern <- paste0("[[:alpha:]]{4}_", c)
    split.chains[c] <- length(grep(pattern, mcmc.files))
  }
  
  # check if there are enough chains to combine
  num.out <- min(split.chains)
  cat("Number of out segments:", num.out, "\n")
  if(num.out <= 1) stop("Not enough segments to combine for at least one chain!", call. = FALSE)
  
  # load the first chain
  load(file.path(mcmc.dir, mcmc.files[1]))
  mcmc.out <- as.matrix(out$out$params)
  iter.run <- nrow(mcmc.out) # number of iterations in each 'out' segment 
  total.iter <- iter.run * num.out
  
  cat("Total iterations in each chain:", total.iter, "\n")
  cat("Chains being combined:", chains, "\n")
  cat("Convergence threshold:", gbr.threshold, "\n")
  cat("Minimum effective sample size:", min.effect.size, "\n")
  
  # combine parameter chains
  params <- list()
  for(c in seq_along(chains)){
    c.load <- chains[c]
    chain.files <- grep(paste0("[[:alpha:]]{4}_", c.load), mcmc.files)
    load(file.path(mcmc.dir, mcmc.files[chain.files[1]]))

    chain <- out$out$params
    for(i in 2:num.out){
      load(file.path(mcmc.dir, mcmc.files[chain.files[i]]))
      chain <- combine.mcmc(mcmc.objects = list(chain, out$out$params))
      if(i %% 5 == 0){
        cat(i, "segments combined\n")
      }
    }
    params[[c]] <- as.mcmc(chain)
    cat("Chain", c.load, "complete\n")
  }
  
  cat("All chains combined\n")
  params <- as.mcmc.list(params)
  
  cat("Calculating PSRF\n")
  GBR.vals <- gelman.diag(params, multivariate = FALSE)
  converge <- max(GBR.vals$psrf) < gbr.threshold
  GBR.vals
  
  cat("Determining burnin\n")
  GBR <- gelman.plot(params)
  burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>gbr.threshold,1,any)),1)+1]
  if(is.na(burnin) & !override){
    stop("Model not converged!", call. = FALSE)
  } else {
    cat("Burnin after:", burnin, "iterations\n")  
  }
  
  # if saving not burned-in output, keep last 1/4 of chains
  if(override & is.na(burnin)){
    burnin <- round(total.iter * 0.75)
  }
  
  ## determine the segment that burnin occurred
  start <- ceiling(burnin / iter.run)
  params.burn <- window(params, start = (start-1)*iter.run+1)
  
  cat("Calculating effective sample size\n")
  effsize <- effectiveSize(params.burn)
  enough.samples <- min(effsize) > min.effect.size
  if(enough.samples | override){
    print(effsize)
    if(!enough.samples){
      cat("Not enough samples for at least one parameter\n")
    } else {
      cat("Minimum effective sample sized reached\n")
    }
  } else {
    print(effsize)
    cat("\n")
    stop("Not enough samples for at least one parameter", call. = FALSE)
  }
  
  ## if there is burnin and convegence and enough samples
  if((!is.na(burnin) & converge & enough.samples) | override){
    
    ## combine predict
    cat("Combining predict output\n")
    predict <- list()
    for(c in seq_along(chains)){
      c.load <- chains[c]
      chain.files <- grep(paste0("[[:alpha:]]{4}_", c.load), mcmc.files)
      load(file.path(mcmc.dir, mcmc.files[chain.files[1]]))
      
      chain <- out$out$predict
      for(i in 2:num.out){
        load(file.path(mcmc.dir, mcmc.files[chain.files[i]]))
        chain <- combine.mcmc(mcmc.objects = list(chain, out$out$predict))
        if(i %% 5 == 0){
          cat(i, "segments combined\n")
        }
      }
      predict[[c]] <- as.mcmc(chain)
      cat("Chain", c.load, "complete\n")
    }
    cat("All chains combined\n")
    
    ## remove burnin and save output
    predict.burn <- as.mcmc.list(predict)
    
    save(params.burn, predict.burn, converge, enough.samples,
         file = file.path(mcmc.dir, "mcmcOut.RData"))
    cat("Burned-in mcmc file saved\n")
    
    params.mat <- as.matrix(params.burn)
    predict.mat <- as.matrix(predict.burn)
    
    iter <- nrow(params.mat)
    thin <- seq(1, iter, by = ceiling(iter / n.mcmc))
    
    params.mat <- params.mat[thin,]
    predict.mat <- predict.mat[thin,]
    
    save(params.mat, predict.mat, converge, enough.samples,
         file = file.path(dir, "mcmcMatrix.RData"))
    cat("Burned-in thinned matrix saved\n")

  }
  
  cat("\n--- DONE ---\n")
  
}

