

do_one_simulation <- function(n,
                              J,
                              distrib,
                              B_multiplier,
                              seed,
                              label,
                              n_boot,
                              # load_tinyvamp = FALSE,
                              folder_name,
                              return_dont_save = FALSE,
                              parallelize = TRUE,
                              verbose = FALSE,
                              return_variance = FALSE){

  if (verbose) {
    message(paste("n = ",n,sep = "",collapse = ""))
    message(paste("J = ",J,sep = "",collapse = ""))
    message(paste("distrib = ",distrib,sep = "",collapse = ""))
    message(paste("B_multiplier = ", B_multiplier,sep = "",collapse = ""))
    message(paste("seed = ",seed,sep = "",collapse = ""))
    message(paste("label = ",label,sep = "",collapse = ""))
  }


  W <- simulate_paper_data(n = n,
                           J = J,
                           B_multiplier = B_multiplier,
                           distrib = distrib,
                           seed = seed)

  poisson_fit <- try(fit_simulation_model(W,"Poisson"))
  reweighted_fit <- try(fit_simulation_model(W, "reweighted_Poisson",
                                             return_variance = return_variance))

  ### Do Bootstrapped LRT for both

  if (verbose) message("Bootstrapping...")
  if(is.list(poisson_fit)){
    poisson_null <- poisson_fit
    poisson_null$B[] <- 0
    poisson_null$B_fixed_indices[] <- TRUE
    poisson_lrt <-
      try(bootstrap_lrt(W = W,
                        fitted_model = poisson_fit,
                        null_param = poisson_null,
                        n_boot = n_boot,
                        parallelize = parallelize,
                        ncores = 5,
                        verbose = verbose))

    poisson_ci <- try(bootstrap_ci(W = W,
                                   fitted_model = poisson_fit,
                                   n_boot = n_boot,
                                   alpha = 0.05,
                                   parallelize = parallelize,
                                   ncores = 5,
                                   seed = seed,
                                   verbose = verbose))


  }

  if(is.list(reweighted_fit)){
    reweighted_null <- reweighted_fit
    reweighted_null$B[] <- 0
    reweighted_null$B_fixed_indices[] <- TRUE
    reweighted_lrt <-
      try(bootstrap_lrt(W,
                        fitted_model = reweighted_fit,
                        null_param = reweighted_null,
                        n_boot = n_boot,
                        parallelize = parallelize,
                        ncores = 5,
                        verbose = verbose))

    reweighted_ci <- try(bootstrap_ci(W = W,
                                      fitted_model = reweighted_fit,
                                      n_boot = n_boot,
                                      alpha = 0.05,
                                      parallelize = parallelize,
                                      ncores = 5,
                                      seed = seed,
                                      verbose = verbose))

  }

  if(!return_dont_save){
    saveRDS(list("poisson_lrt" = poisson_lrt,
                 "poisson_ci" = poisson_ci,
                 "reweighted_lrt" = reweighted_lrt,
                 "reweighted_ci" = reweighted_ci),
            paste(folder_name,"/",
                  paste(label,"n",n,"J", J, distrib,
                        "Bmult", B_multiplier,"nboot", n_boot,"sim",
                        seed,
                        sep = "_",collapse = "_"),
                  sep = "",collapse = ""))
  } else{
    return(list("poisson_lrt" = poisson_lrt,
                "poisson_ci" = poisson_ci,
                "reweighted_lrt" = reweighted_lrt,
                "reweighted_ci" = reweighted_ci))
  }
}
