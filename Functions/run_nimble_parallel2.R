#' function for running nimble chains in parallel
#' this function checks convergence and will keep sampling if not converged
#' 
#' @param cl cluster built in main script
#' @param model nimble model code
#' @param constants nimble constants list
#' @param data nimble data list
#' @param inits initial values function
#' @param monitor vector of variables to monitor
#' @param file.name full file path to where you want intermediate output to be saved. If NA nothing saves.
#' @param use.dzip using zero-inflated Poisson?
#' @param n.ens the number of total iterations to save
#' @param check.interval before convergence, how often do you want to print/save output? 
#' @param thin thinning interval
#' @param n.iter the number of iterations for each intermediate sampling step
#' @param n.burnin the number of burnin iterations
#' @param psrf.max the maximum value allowed for convergence (R_hat threshold)
#' @param min.effective.size the minimum number of effective samples allowed
#' @param max.iter the maximum number of iterations before sampling stops
#' @param calculate are we going to save latent state samples? default TRUE

run_nimble_parallel <- function(cl, model, constants, data, inits, monitor, file.name = NA,
                                use.dzip = FALSE, n.ens = 5000, check.interval = 10, thin = 1,
                                n.iter = 50000, n.burnin = 5000, psrf.max = 1.1, min.effective.size = 2000,
                                max.iter = 3e6, calculate = TRUE){
  library(parallel)
  library(nimble)
  library(coda)
  source("Functions/nimble_functions.R")
  n.cores <- length(cl) # number of cores used
  
  # export everything we need to cluster nodes
  clusterExport(cl, 
                c("model", "constants", "data", "if_else_nimble", "calculate",
                  "n.iter", "n.burnin", "monitor", "thin", "n.cores"),
                envir = environment()) 
  
  if(use.dzip){
    source("Functions/ZIP.R")
    assign("dZIP", dZIP, envir = .GlobalEnv)
    assign("rZIP", rZIP, envir = .GlobalEnv)
    clusterExport(cl, 
                  c("dZIP", "rZIP"),
                  envir = environment())  
  }
  
  # export inits to clusters
  for(j in seq_along(cl)){
    set.seed(j)
    init <- inits()
    # print(init)
    clusterExport(cl[j], "init", envir = environment())
  }
  
  message("Running mcmc...")
  start.time <- Sys.time()
  out <- clusterEvalQ(cl, { # sample on each cluster
    library(nimble)
    library(nimbleEcology)
    library(coda)
    model.rw <- nimbleModel(model,
                            constants = constants,
                            data = data,
                            inits = init,
                            calculate = calculate)
    cModel.rw <- compileNimble(model.rw)
    mcmcConf <- configureMCMC(cModel.rw, 
                              monitors = monitor,
                              thin = thin,
                              useConjugacy = TRUE)
    
    # mcmcConf$removeSampler(c("phi.a", "phi.p"))
    # mcmcConf$addSampler(target = c("phi.a", "phi.p"), type = "RW_block")

    # mcmcConf$removeSampler(c("psi.ap", "psi.pa"))
    # mcmcConf$addSampler(target = c("psi.ap", "psi.pa"), type = "RW_block")

    # mcmcConf$removeSampler(c("p.a", "p.p", "p.d", "p.u"))
    # mcmcConf$addSampler(target = c("p.a", "p.p", "p.d", "p.u"), type = "RW_block")
    
    mcmcBuild <- buildMCMC(mcmcConf)
    compMCMC <- compileNimble(mcmcBuild)
    compMCMC$run(niter = n.iter, nburnin = n.burnin)
    return(as.mcmc(as.matrix(compMCMC$mvSamples)))
  })
  delta.time <- Sys.time() - start.time
  message(paste(n.iter,"iterations in"))
  print(delta.time)
  out.mcmc <- as.mcmc.list(out)
  
  states.monitor <- if_else(any(grepl("x[", colnames(out.mcmc[[1]]), fixed = TRUE)),
                            TRUE, FALSE)
  if(states.monitor){
    check.mcmc <- out.mcmc[,-grep("x[", colnames(out.mcmc[[1]]), fixed = TRUE), drop = TRUE]
  } else {
    check.mcmc <- out.mcmc
  }
  
  # check for convergence
  g.diag <- gelman.diag(check.mcmc, multivariate = FALSE)$psrf
  convergence <- max(g.diag[,"Upper C.I."]) < psrf.max
  message(paste("Convergence:", convergence)) 
  if(is.na(convergence)){
    print(g.diag)
    convergence <- FALSE
  } 
  
  # effective sample size
  effect.sizes <- effectiveSize(check.mcmc)
  enough.samples <- min(effect.sizes) >= min.effective.size
  message(paste("Enough Effective Samples:", enough.samples)) 
  
  counter <- 1
  total.iter <- counter * n.iter
  resetMV <- FALSE
  while(!convergence | !enough.samples){
    counter <- counter + 1
    total.iter <- counter * n.iter
    
    message("Running mcmc...")
    start.time <- Sys.time()
    clusterExport(cl, "resetMV", envir = environment())
    out2 <- clusterEvalQ(cl, {
      compMCMC$run(n.iter, reset = FALSE, resetMV = FALSE)
      return(as.mcmc(as.matrix(compMCMC$mvSamples)))
    })
    out.mcmc <- as.mcmc.list(out2)
    
    delta.time <- Sys.time() - start.time
    message(paste(n.iter,"iterations in"))
    print(delta.time)
    
    if(states.monitor){
      check.mcmc <- out.mcmc[,-grep("x[", colnames(out.mcmc[[1]]), fixed = TRUE), drop = TRUE]
    } else {
      check.mcmc <- out.mcmc
    }
    
    g.diag <- gelman.diag(check.mcmc, multivariate = FALSE)$psrf
    convergence <- max(g.diag[,"Upper C.I."]) < psrf.max
    message(paste("Convergence:", convergence))  
    if(is.na(convergence)) convergence <- FALSE
    
    effect.sizes <- effectiveSize(check.mcmc)
    enough.samples <- min(effect.sizes) >= min.effective.size
    message(paste("Enough Effective Samples:", enough.samples))
    
    message(paste("Total iterations:", total.iter))
    if(total.iter > max.iter) break
    
    if(counter %% check.interval == 0) {
      message("Convergence check:")
      print(g.diag)
      message("Effective sample sizes:")
      print(effect.sizes)
      if(!is.na(file.name)){
        save.ls <- list(samples = as.mcmc.list(check.mcmc),
                        convergence = convergence)
        
        message("Writing output...")
        save(save.ls, file = file.name)
        message("DONE")
      }
      if(!convergence) resetMV <- TRUE
    } else {
      resetMV <- FALSE
    }
  }
  
  if(convergence & enough.samples){ #
    GBR <- gelman.plot(check.mcmc, ask = FALSE)
    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
    message(paste("Burnin after", burnin, "iterations"))
    
    # sometimes burnin happens at the end of the chain, need to keep sampling if that happens
    if(total.iter - burnin < 1000){
      message("Additional iterations because burnin is too close to the end of the chain")
      out3 <- clusterEvalQ(cl, {
        compMCMC$run(n.iter, reset = FALSE, resetMV = TRUE)
        return(as.mcmc(as.matrix(compMCMC$mvSamples)))
      })
      samples <- as.mcmc.list(out.3)
    } else{
      samples <- window(as.mcmc.list(out.mcmc), start = burnin)
    }  
    
    # return a thinned matrix (raw mcmc objects can get big)
    samples <- as.matrix(samples)
    if(n.ens > nrow(samples)){
      thin.seq <- round(seq(1, nrow(samples), length.out = n.ens)) # sequence of samples to keep
      samples <- samples[thin.seq,]  
    }
    
  } else {
    samples <- as.mcmc.list(out.mcmc)
  }
  
  return(list(samples = samples,
              convergence = convergence,
              enough.samples = enough.samples)) 
  
}
