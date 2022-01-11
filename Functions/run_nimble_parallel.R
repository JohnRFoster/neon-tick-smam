#' function for running nimble chains in parallel
#' this function checks concergence and will keep sampling if not convereged
#' 
#' @param cl cluster built in main script
#' @param model nimble model code
#' @param constants nimble constants list
#' @param data nimble data list
#' @param inits initial values function
#' @param monitor vector of variables to monitor
#' @param file.name full file path to where you want intermediate output to be saved
#' @param use.dzip using zer-inflated Poisson?
#' @param n.ens the number of total iterations to save
#' @param check.interval before convergence, how often do you want to print/save output? 
#' @param thin thinning interval
#' @param n.iter the number of iterations to each intermediate sampling step
#' @param n.burnin the number of burnin iterations
#' @param psrf.max the maximum value allowed for convergence (R_hat threshold)
#' @param max.iter the maximum number of iterations before sampling stops
#' @param check.params.only check for convergence on parameters only (TRUE) or parameters and states (FALSE)

run_nimble_parallel <- function(cl, model, constants, data, inits, monitor, file.name = NA,
                                use.dzip = FALSE, n.ens = 5000, check.interval = 10, thin = 1,
                                n.iter = 50000, n.burnin = 5000, psrf.max = 1.1, max.iter = 3e6,
                                check.params.only = FALSE){
  library(parallel)
  library(nimble)
  library(coda)
  
  n.cores <- length(cl) # number of cores used
  
  # export everything we need to cluster nodes
  if(use.dzip){
    source("Functions/ZIP.R")
    assign("dZIP", dZIP, envir = .GlobalEnv)
    assign("rZIP", rZIP, envir = .GlobalEnv)
    clusterExport(cl, 
                  c("model", "constants", "data", "n.iter", "n.burnin",
                    "dZIP", "rZIP", "monitor", "thin", "n.cores"),
                  envir = environment())  
  } else {
    clusterExport(cl, 
                  c("model", "constants", "data", 
                    "n.iter", "n.burnin", "monitor", "thin", "n.cores"),
                  envir = environment()) 
  }
  
  # export inits to clusters
  for(j in seq_along(cl)){
    set.seed(j)
    init <- inits()
    # print(init)
    clusterExport(cl[j], "init", envir = environment())
  }
  
  message("Running mcmc")
  out <- clusterEvalQ(cl, { # sample on each cluster
    library(nimble)
    library(coda)
    model.rw <- nimbleModel(model,
                            constants = constants,
                            data = data,
                            inits = init)
    cModel.rw <- compileNimble(model.rw)
    mcmcConf <- configureMCMC(cModel.rw, 
                              monitors = monitor,
                              thin = thin)
    # mcmcConf$addSampler(target = c("tau.obs[1]", "tau.obs[2]"), type = "RW_block")
    mcmcBuild <- buildMCMC(mcmcConf)
    compMCMC <- compileNimble(mcmcBuild)
    out.1 <- runMCMC(compMCMC, niter = n.iter, nburnin = n.burnin)
    return(as.mcmc(out.1))
  })
  out.mcmc <- as.mcmc.list(out)
  
  if(check.params.only){  
    message("Checking convergence on parameters only")
    check.mcmc <- out.mcmc[,-grep("x[", colnames(out.mcmc[[1]]), fixed = TRUE), drop = TRUE]
  } else {
    message("Checking convergence on parameters and states")
    check.mcmc <- out.mcmc
  }
  
  g.diag <- gelman.diag(check.mcmc, multivariate = FALSE)$psrf
  convergence <- max(g.diag[,"Upper C.I."]) < psrf.max
  message(paste("Convergence:", convergence)) 
  if(is.na(convergence)) print(g.diag)
  
  counter <- 1
  total.iter <- counter * n.iter
  while(!convergence & total.iter < max.iter){
    counter <- counter + 1
    total.iter <- counter * n.iter
    
    message("Running mcmc")
    out2 <- clusterEvalQ(cl, {
      out1 <- runMCMC(compMCMC, niter = n.iter)
      return(as.mcmc(out1))
    })
    out.mcmc.update1 <- as.mcmc(out2)
    
    out.mcmc.bind <- list()
    for (i in seq_len(n.cores)) {
      out.mcmc.bind[[i]] <- mcmc(rbind(out.mcmc[[i]], out.mcmc.update1[[i]]))
    }
    out.mcmc <- as.mcmc.list(out.mcmc.bind)
    if(check.params.only){  
      message("Checking convergence on parameters only")
      check.mcmc <- out.mcmc[,-grep("x[", colnames(out.mcmc[[1]]), fixed = TRUE), drop = TRUE]
    } else {
      message("Checking convergence on parameters and states")
      check.mcmc <- out.mcmc
    }
    g.diag <- gelman.diag(check.mcmc, multivariate = FALSE)$psrf
    if(counter %% check.interval == 0) {
      print(g.diag)
      if(!is.na(file.name)){
        save.ls <- list(samples = as.mcmc.list(out.mcmc),
                        convergence = convergence)
        
        message("Writing output...")
        save(save.ls, file = file.name)
        message("DONE")
      }
    } 
    convergence <- max(g.diag[,"Upper C.I."]) < psrf.max
    message(paste("Convergence:", convergence))  
    message(paste("Total iterations:", total.iter))
    
  }
  
  if(convergence){ #
    GBR <- gelman.plot(check.mcmc, ask = FALSE)
    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
    message(paste("Burnin after", burnin, "iterations"))
    
    # sometimes burnin happens at the end of the chain, need to keep sampling if that happens
    if(total.iter - burnin < 1000){
      message("Additional iterations because burnin is too close to the end of the chain")
      out3 <- clusterEvalQ(cl, {
        out <- runMCMC(compMCMC, niter = n.iter/2)
        return(as.mcmc(out))
      })
      samples <- as.mcmc.list(out)
    } else{
      samples <- window(as.mcmc.list(out.mcmc), start = burnin)
    }  
    
    # return a thinned matrix (raw mcmc objects can get big)
    samples <- as.matrix(samples)
    thin.seq <- round(seq(1, nrow(samples), length.out = n.ens)) # sequence of samples to keep
    samples <- samples[thin.seq,]
    
  } else {
    samples <- as.mcmc.list(out.mcmc)
  }
  
  
  
  return(list(samples = samples,
              convergence = convergence)) 
  
}
