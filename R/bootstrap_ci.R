#' Apply the Bayesian subsampled bootstrap to a fitted tinyvamp model
#' 
#' @import stats
#' @import parallel
#' 
#' 
#' @export
bootstrap_ci <- function(W,
                         fitted_model,
                         n_boot,
                         m = NULL,
                         alpha = 0.05,
                         parallelize = FALSE,
                         ncores = 5,
                         seed = NULL,
                         return_models = FALSE,
                         verbose = FALSE,
                         adjust = FALSE
                         
){
  n <- nrow(W)
  J <- ncol(W)
  
  if(is.null(m)){
    m <- sqrt(n)
  }
  if(is.null(seed)){
    seed <- 0
  }
  
  if(is.null(fitted_model$wts)){
    fitted_model$wts <- rep(1,n*J)
  }
  wts <- fitted_model$wts
  set.seed(seed)
  boot_seeds <- sample(1:1e8,n_boot)
  
  boot_results <- vector(n_boot,
                         mode = "list")
  
  if(!parallelize){
    for(boot_iter in 1:n_boot){
      print(boot_iter)
      set.seed(boot_seeds[boot_iter])
      boot_weights <- rgamma(n,m/n)
      boot_weights <- boot_weights/sum(boot_weights)
      boot_weights <- rep(boot_weights, each = J)
      boot_weights <- boot_weights*wts
      boot_weights <- J*boot_weights/sum(boot_weights)
      
      boot_model <- estimate_parameters(W = W,
                                        X = fitted_model$X,
                                        Z = fitted_model$Z,
                                        Z_tilde = fitted_model$Z_tilde,
                                        Z_tilde_gamma_cols =
                                          fitted_model$Z_tilde_gamma_cols,
                                        gammas =
                                          fitted_model$gammas,
                                        gammas_fixed_indices =
                                          fitted_model$gammas_fixed_indices,
                                        P = fitted_model$P,
                                        P_fixed_indices =
                                          fitted_model$P_fixed_indices,
                                        B = fitted_model$B,
                                        B_fixed_indices =
                                          fitted_model$B_fixed_indices,
                                        X_tilde = fitted_model$X_tilde,
                                        P_tilde = fitted_model$P_tilde,
                                        P_tilde_fixed_indices =
                                          fitted_model$P_tilde_fixed_indices,
                                        gamma_tilde = fitted_model$gamma_tilde,
                                        gamma_tilde_fixed_indices =
                                          fitted_model$gamma_tilde_fixed_indices,
                                        alpha_tilde =
                                          fitted_model$alpha_tilde,
                                        Z_tilde_list =
                                          fitted_model$Z_tilde_list,
                                        barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                                        barrier_scale = 10, #increments for value of barrier penalty
                                        max_barrier = 1e12, #maximum value of barrier_t
                                        initial_conv_tol = 1000,
                                        final_conv_tol = 0.1,
                                        constraint_tolerance = 1e-10,
                                        hessian_regularization = 0.01,
                                        criterion = "Poisson",
                                        profile_P = TRUE,
                                        verbose = verbose,
                                        wts = boot_weights,
                                        profiling_maxit = 25)
      
      boot_results[[boot_iter]] <- boot_model
      
    }
  }
  if(parallelize){
    boot_weights <- lapply(1:n_boot,
                           function(x){
                             # set.seed(x)
                             bwts <- rgamma(n,m/n)
                             bwts <- rep(bwts, each = J)
                             bwts <- bwts*wts
                             bwts <- J*bwts/sum(bwts)
                             return(bwts)
                           })
    
    boot_results <-
      parallel::mclapply(1:n_boot,
                         function(k)
                           do_one_boot(W = W,
                                       fitted_model =
                                         fitted_model,
                                       m = m,
                                       seed = boot_seeds[k],
                                       boot_weights = boot_weights[[k]]),
                         mc.cores = ncores,
                         mc.set.seed = TRUE)
  }
  
  boot_matrix <-
    do.call(cbind,
            lapply(1:n_boot,function(k) matrix(sqrt(m)*(
              boot_results[[k]]$varying$value - fitted_model$varying$value),
              ncol = 1)))
  
  if(adjust){
    num_nonsillyparams <- sum(fitted_model$varying$param != "gamma")
    adjust_factor <- (n*J)/(n*J - num_nonsillyparams)
    boot_matrix <- boot_matrix*sqrt(adjust_factor)
  }
  
  lower_boot_quantiles <- apply(boot_matrix,1, function(x)
    quantile(x, alpha/2))
  
  upper_boot_quantiles <- apply(boot_matrix,1,function(x)
    quantile(x,1-alpha/2))
  
  summary_df <- fitted_model$varying
  
  summary_df$lower_ci <- summary_df$value - (1/sqrt(n))*upper_boot_quantiles
  summary_df$upper_ci <- summary_df$value - (1/sqrt(n))*lower_boot_quantiles
  
  summary_df$lower_ci[summary_df$param %in% c("P","P_tilde")] <-
    pmax(summary_df$lower_ci[summary_df$param %in% c("P","P_tilde")],
         0)
  
  summary_df$upper_ci[summary_df$param %in% c("P","P_tilde")] <-
    pmin(summary_df$upper_ci[summary_df$param %in% c("P","P_tilde")],
         1)
  
  if(return_models){
    return(list("ci" = summary_df,
                "bootstrapped_models" = boot_results))
  } else{
    return(list("ci" = summary_df,
                "bootstrapped_models" = NULL))
  }
  
}
