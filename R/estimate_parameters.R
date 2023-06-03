
#' Fit tinyvamp model to HTS microbiome data
#'
#' This function fits a model to HTS microbiome data that allows for estimation of
#' detection efficiency effects as well as modeling of spurious read sources
#' (e.g., contamination).
#'
#'
#' @param W An \eqn{n \times J} matrix of numeric HTS output (e.g., read counts, coverages, etc.)
#' @param X The sample efficiency design -- an \eqn{n \times p} matrix
#' @param Z The sample-specimen design -- an \eqn{n \times K} matrix whose \eqn{ij}-th entry
#' indicates the proportional contribution of specimen \eqn{j} to sample \eqn{i}. Rows must
#' sum to 1 or be identically 0.
#' @param Z_tilde The spurious read design -- an \eqn{n \times \tilde{K}} matrix where
#' \eqn{\tilde{K}} is the number of spurious read sources modeled.
#' @param Z_tilde_gamma_cols A numeric vector containing the columns of Z_tilde which should be
#' multiplied by exp(gamma).
#' @param gammas A numeric vector of length n of starting values for read intensity parameter gamma
#' @param gammas_fixed_indices A logical vector of length n whose \eqn{i}-th entry is TRUE if the
#' \eqn{i}-th entry of gamma should be treated as fixed and known, and FALSE otherwise
#' @param P A \eqn{K \times J} numeric matrix giving initial values for the relative abundance matrix.
#' @param P_fixed_indices P_fixed_indices A \eqn{K \times J} logical matrix specifying any entries of P that are known. If known, the corresponding values from \code{P} will be treated as the fixed, known values.
#' @param B A \eqn{p \times J} numeric matrix giving initial values for the sample efficiencies.
#' @param B_fixed_indices A \eqn{p \times J} logical matrix specifying any entries of B that are known. If known, the corresponding values from \code{B} will be treated as the fixed, known values.
#' @param X_tilde A \eqn{\tilde{K} \times p} matrix giving the spurious read source efficiency design matrix
#' @param P_tilde A \eqn{\tilde{K} \times J} numeric matrix giving initial values for the spurious read source relative abundances.
#' @param P_tilde_fixed_indices A \eqn{\tilde{K} \times J} logical matrix indicating if the \eqn{(i,j)}th entry of \code{P_tilde} should be treated as fixed and known.
#' @param gamma_tilde A numeric vector of length \eqn{\tilde{K}} of starting values for spurious read intensity parameter gamma_tilde
#' @param gamma_tilde_fixed_indices A logical vector of length \eqn{\tilde{K}} whose \eqn{i}-th entry is TRUE if the
#' \eqn{i}-th entry of gamma_tilde should be treated as fixed and known, and FALSE otherwise
#' @param barrier_t Starting value of reciprocal barrier penalty coef. Defaults to 1.
#' @param barrier_scale Increments for value of barrier penalty. Defaults to 10.
#' @param max_barrier Maximum value of barrier_t. Defaults to 1e12.
#' @param constraint_tolerance The tolerance for the augmented Lagrangian algorithm. Final estimates of P are relative abundances to within \code{constraint_tolerance} of 1, i.e., abs(sum p_{kj} - 1) <  \code{constraint_tolerance}. Defaults to 1e-10.
#' @param hessian_regularization The second step of optimization involves a quadratic approximation to the likelihood, for which we use a modified Taylor series for stability. This is the constant that dampens the second term. Defaults to 0.01. 
#' @param criterion Should the algorithm return the Poisson maximum likelihood estimates or the reweighted Poisson maximum likelihood estimates? Options are "Poisson" or "reweighted_Poisson". 
#' @param profile_P Defaults to TRUE Run profiling step after barrier algorithm has run? If TRUE, this step is performed, possibly setting some estimated relative abundances in P equal to zero. If FALSE, profiling step is skipped and back-transformed log-ratio parameter estimated via barrier algorithm is returned for P.
#' @param profiling_maxit Maximum number of iterations to run profiling step in P for (default is 25).
#' @param wts Weights for reweighting the likelihood contributions. This is usually done to improve efficiency. Defaults to NULL. We compute the weights for you even if you choose \code{criterion = "reweighted_Poisson"}. 
#' @param verbose Do you want to know what I'm doing? Defaults to FALSE. 
#' @param bootstrap_failure_cutoff Defaults to NULL.
#' @param return_variance Defaults to FALSE.
#' 
#' @return A list containing TODO
#' 
#' @author David Clausen
#'
#' @import cir
#'
#' @export
estimate_parameters <- function(W,
                                X,
                                Z,
                                Z_tilde = NULL,
                                Z_tilde_gamma_cols,
                                gammas,
                                gammas_fixed_indices,
                                P,
                                P_fixed_indices,
                                B,
                                B_fixed_indices,
                                X_tilde,
                                P_tilde,
                                P_tilde_fixed_indices,
                                gamma_tilde,
                                gamma_tilde_fixed_indices,
                                alpha_tilde = NULL,
                                Z_tilde_list = NULL,
                                barrier_t = 1, 
                                barrier_scale = 10, 
                                max_barrier = 1e12, 
                                initial_conv_tol = 1000,
                                final_conv_tol = 0.1,
                                constraint_tolerance = 1e-10,
                                hessian_regularization = 0.01,
                                criterion = "Poisson",
                                profile_P = TRUE,
                                barrier_maxit = 500,
                                profiling_maxit = 25,
                                wts = NULL,
                                verbose = FALSE,
                                bootstrap_failure_cutoff = NULL,
                                return_variance = FALSE
){
  
  if(!(criterion %in% c("Poisson","reweighted_Poisson"))){
    stop("Argument `criterion` must be equal to
`Poisson` or `reweighted_Poisson`.")
  }
  
  
  n <- nrow(W)
  J <- ncol(W)
  
  ### add small amount to all estimated relative abundances to avoid zeroes
  P[!P_fixed_indices] <- P[!P_fixed_indices] + 0.1/J
  P_tilde[!P_tilde_fixed_indices] <-
    P_tilde[!P_tilde_fixed_indices] + 0.1/J
  
  min_regularization <- hessian_regularization
  
  
  if(sum(P[!P_fixed_indices] == 0)>0){
    P[!P_fixed_indices] <- P[!P_fixed_indices] + 0.01/J
  }
  
  if(sum(P_tilde[!P_tilde_fixed_indices]==0)>0){
    P_tilde[!P_tilde_fixed_indices] <-
      P_tilde[!P_tilde_fixed_indices] + 0.01/J
  }
  
  #mi initial alpha_tilde, Z_tilde, Z_tilde_list checks
  if(!is.null(alpha_tilde)){
    if(is.null(Z_tilde_list)){
      stop("If alpha_tilde is included in model, components of Z_tilde
to be multiplied by exp(alpha_tilde) must be provided in
Z_tilde_list")
    }
  }
  if(is.null(alpha_tilde)){
    if(!is.null(Z_tilde_list)){
      stop("If Z_tilde_list is provided, alpha_tilde must also be provided")
    }
  }
  if(is.null(Z_tilde)){
    if(is.null(Z_tilde_list)){
      stop("Z_tilde may only be NULL if Z_tilde_list is not NULL")
    }
  } else{
    if(!is.null(Z_tilde_list)){
      warning("Z_tilde argument is ignored when Z_tilde_list is provided.
(Z_tilde is constructed via parametrization in terms of components of
alpha_tilde and matrices in Z_tilde_list.)")
    }
  }
  
  
  if(criterion == "GMM"){
    if(verbose){message("Fitting initial Poisson estimator")}
    poisson_fit <-
      estimate_parameters(W = W,
                          X = X,
                          Z = Z,
                          Z_tilde = Z_tilde,
                          Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                          gammas = gammas,
                          gammas_fixed_indices = gammas_fixed_indices,
                          P = P,
                          P_fixed_indices = P_fixed_indices,
                          B = B,
                          B_fixed_indices = B_fixed_indices,
                          X_tilde = X_tilde,
                          P_tilde = P_tilde,
                          P_tilde_fixed_indices = P_tilde_fixed_indices,
                          gamma_tilde = gamma_tilde,
                          gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
                          alpha_tilde = alpha_tilde,
                          Z_tilde_list = Z_tilde_list,
                          barrier_t = barrier_t,
                          barrier_scale = barrier_scale,
                          max_barrier = max_barrier,
                          initial_conv_tol = initial_conv_tol,
                          final_conv_tol = final_conv_tol,
                          constraint_tolerance = constraint_tolerance,
                          hessian_regularization = hessian_regularization,
                          criterion = "Poisson",
                          profile_P = FALSE,
                          profiling_maxit = 25,
                          wts = wts)
    
    
    poisson_fit_params <- dataframes_to_parameters(fixed = poisson_fit$fixed,
                                                   varying = poisson_fit$varying)
    
    #update values of parameters
    gammas <- poisson_fit_params$gammas
    P <- poisson_fit_params$P
    B <- poisson_fit_params$B
    P_tilde <- poisson_fit_params$P_tilde
    gamma_tilde <- poisson_fit_params$gamma_tilde
    alpha_tilde <- poisson_fit_params$alpha_tilde
    
  }
  if(criterion == "reweighted_Poisson"){
    if(is.null(wts)){
      wts <- rep(1,n*J)
    }
    
    poisson_fit <-
      estimate_parameters(W = W,
                          X = X,
                          Z = Z,
                          Z_tilde = Z_tilde,
                          Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                          gammas = gammas,
                          gammas_fixed_indices = gammas_fixed_indices,
                          P = P,
                          P_fixed_indices = P_fixed_indices,
                          B = B,
                          B_fixed_indices = B_fixed_indices,
                          X_tilde = X_tilde,
                          P_tilde = P_tilde,
                          P_tilde_fixed_indices = P_tilde_fixed_indices,
                          gamma_tilde = gamma_tilde,
                          gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
                          alpha_tilde = alpha_tilde,
                          Z_tilde_list = Z_tilde_list,
                          barrier_t = barrier_t,
                          barrier_scale = barrier_scale,
                          max_barrier = max_barrier,
                          initial_conv_tol = initial_conv_tol,
                          final_conv_tol = final_conv_tol,
                          constraint_tolerance = constraint_tolerance,
                          hessian_regularization = hessian_regularization,
                          criterion = "Poisson",
                          verbose = verbose,
                          profile_P = TRUE,
                          profiling_maxit = 25,
                          wts = wts)
    
    poisson_fit_params <- dataframes_to_parameters(fixed = poisson_fit$fixed,
                                                   varying = poisson_fit$varying)
    
    #update values of parameters
    gammas <- poisson_fit_params$gammas
    P <- poisson_fit_params$P
    B <- poisson_fit_params$B
    P_tilde <- poisson_fit_params$P_tilde
    gamma_tilde <- poisson_fit_params$gamma_tilde
    alpha_tilde <- poisson_fit_params$alpha_tilde
    
    poisson_means <- meaninate(gammas = poisson_fit$gammas,
                               B = poisson_fit$B,
                               X = X,
                               Z = Z,
                               P = P,
                               X_tilde = X_tilde,
                               Z_tilde = Z_tilde,
                               Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                               P_tilde = poisson_fit$P_tilde,
                               gamma_tilde = poisson_fit$gamma_tilde,
                               alpha_tilde = poisson_fit$alpha_tilde,
                               Z_tilde_list = poisson_fit$Z_tilde_list)
    
    W_long <- lapply(1:n,function(i) as.numeric(W[i,]))
    W_long <- do.call(c,W_long)
    means_long <- lapply(1:n, function(i) as.numeric(poisson_means[i,]))
    means_long <- do.call(c,means_long)
    
    squerror_long <- (W_long - means_long)^2
    
    pre_wts <- tryCatch(cir::cirPAVA(y = squerror_long, x = means_long, wt= wts), 
                        error = function(c) { 
                          if (verbose) 
                            message("Fitted means are the same, breaking cirPAVA...\nNot to worry! Jitter and refit...\n")
                        }
    )
    
    if(is.null(pre_wts)) { # i.e., did the above error
      # add tiny amount of jitter to means_long so that all values are unique;
      # suspected bug in cir package otherwise causes errors at least
      # some of the time
      means_long <- means_long + runif(length(means_long),0,1e-8)
      pre_wts <- try(cir::cirPAVA(y = squerror_long, x = means_long,
                                  wt= wts))
      
      if(inherits(pre_wts,"try-error")){
        stop("Fatal error in cirPAVA")
      }
    }
    
    
    
    
    inv_wts <- numeric(length(W_long))
    inv_wts[order(means_long)] <- pre_wts +1
    if(return_variance){
      variance_df <- data.frame("squerror" = squerror_long,
                                "mean" = means_long,
                                "estd_var" = inv_wts)
    }
    # a <- min(inv_wts[inv_wts >0 & means_long>0])/
    #   min(means_long[inv_wts >0 & means_long>0])
    #
    inv_wts <- inv_wts/(means_long + 1)
    # inv_wts[means_long ==0] <- a
    # inv_wts <- sqrt(inv_wts)
    
    # new_wts <- 1/(W_long - means_long)^2
    
    wt_total <- sum(wts) #preserve number weights sum to
    wts <- wts/inv_wts
    # wts <- wts*new_wts
    wts <- wt_total*wts/sum(wts)
    
    criterion <- "Poisson"
    
    re_pois <- TRUE
    
    ### add small amount to all estimated relative abundances to avoid zeroes
    poisson_fit$varying$value[
      poisson_fit$varying$param %in% c("P","P_tilde")] <-
      poisson_fit$varying$value[
        poisson_fit$varying$param %in% c("P","P_tilde")] + 0.1/J
    
    poisson_fit_params <- dataframes_to_parameters(fixed = poisson_fit$fixed,
                                                   varying = poisson_fit$varying)
    
    #update values of parameters
    gammas <- poisson_fit_params$gammas
    P <- poisson_fit_params$P
    B <- poisson_fit_params$B
    P_tilde <- poisson_fit_params$P_tilde
    gamma_tilde <- poisson_fit_params$gamma_tilde
    alpha_tilde <- poisson_fit_params$alpha_tilde
    
    
  } else{
    re_pois <- FALSE
  }
  
  
  
  #determine number of barrier steps to take
  nsteps <- ceiling(log(max_barrier/barrier_t)/log(barrier_scale)) + 1
  
  #determine sequence of convergence tolerances
  tolerances <- exp(seq(log(initial_conv_tol),log(final_conv_tol), length.out = nsteps))
  ### Add relative abundance checks for P, P_tilde with fixed elements
  
  ### Add check: no RA parameters to be estimated may have initial value 0
  
  # determine whether any parameters need to be log-ratio transformed for estimation
  
  
  ### set up parameters according to
  ### those to be estimated and those fixed
  ### at predetermined values
  
  parameter_dfs <- parameters_to_dataframes(P,
                                            P_fixed_indices,
                                            P_tilde,
                                            P_tilde_fixed_indices,
                                            B,
                                            B_fixed_indices,
                                            gammas,
                                            gammas_fixed_indices,
                                            gamma_tilde,
                                            gamma_tilde_fixed_indices,
                                            alpha_tilde)
  
  
  varying_df <- parameter_dfs$varying
  fixed_df <- parameter_dfs$fixed
  varying_lr_df <- ra_to_lr(varying_df)
  
  # message("created parameter dfs")
  
  K <- max(c(varying_df$k[varying_df$param == "P"],
             fixed_df$k[fixed_df$param == "P"]))
  
  fixed_P_multipliers <- sapply(1:K, function(k)
    1 - sum(fixed_df$value[fixed_df$param == "P"&
                             fixed_df$k ==k]))
  
  K_tilde <- max(c(varying_df$k[varying_df$param == "P_tilde"],
                   fixed_df$k[fixed_df$param == "P_tilde"]))
  
  fixed_P_tilde_multipliers <- sapply(1:K_tilde, function(k)
    1 - sum(fixed_df$value[fixed_df$param == "P_tilde"&
                             fixed_df$k ==k]))
  
  # message('Allegedly, we have calculated fixed_P_tilde_multipliers. Here they are:')
  # print(fixed_P_tilde_multipliers)
  
  # message("stored K, K_tilde, calculated multipliers")
  
  
  # create matrices to track rho-P and rho_tilde-P_tilde relationships
  Ak_list <- get_Ak_list(fixed_df,
                         varying_df,
                         varying_lr_df)
  
  A_tilde_k_list <- get_A_tilde_k_list(fixed_df,
                                       varying_df,
                                       varying_lr_df)
  
  # message("created Ak_list and A_tilde_k_list")
  
  #calculate at outset of optimization  pass as argument to mean_jac_lr, etc.
  which_k_p <- sapply(1:K, function(k) ifelse(is.null(Ak_list[[k]]),
                                              NA, k))
  
  which_k_p <- which_k_p[!is.na(which_k_p)]
  
  #calculate at outset of optimization  pass as argument to mean_jac_lr, etc.
  which_k_p_tilde <- sapply(1:K_tilde,
                            function(k) ifelse(
                              is.null(A_tilde_k_list[[k]]),
                              NA,k
                            ))
  
  which_k_p_tilde <- which_k_p_tilde[!is.na(which_k_p_tilde)]
  
  # message("saved which_k_p and which_k_p_tilde")
  
  #calculate at outset of optimization
  which_B_rows <- unique(varying_df$k[varying_df$param == "B"])
  which_B_rows <- which_B_rows[order(which_B_rows)]
  
  #calculate at outset of optimization
  which_B_keep <- lapply(which_B_rows,
                         function(k) sapply(1:(J - 1),
                                            function(j)
                                              j %in% varying_lr_df$j[
                                                varying_lr_df$param == "B" &
                                                  varying_lr_df$k == k]
                         ))
  which_B_keep <- do.call(rbind,which_B_keep)
  
  # message("saved which_B_keep")
  
  which_gammas <- unique(varying_df$k[varying_df$param == "gamma"])
  
  which_gamma_tilde <- unique(varying_df$k[varying_df$param == "gamma_tilde"])
  
  which_unconstrained <- varying_lr_df$param %in% c("B","gamma","gamma_tilde",
                                                    "alpha_tilde")
  which_rho <- varying_lr_df$param %in% c("rho")
  which_rho_tilde <- varying_lr_df$param %in% c("rho_tilde")
  npar <- nrow(varying_lr_df)
  
  # iterate until some criterion is met
  if(verbose){message("evaluating criterion")}
  crit_value <- evaluate_criterion_lr(
    W = W,
    X = X,
    Z = Z,
    Z_tilde = Z_tilde,
    Z_tilde_list = Z_tilde_list,
    Z_tilde_gamma_cols = Z_tilde_gamma_cols,
    X_tilde = X_tilde,
    fixed_df = fixed_df,
    varying_df = varying_df,
    varying_lr_df = varying_lr_df,
    barrier_t = barrier_t,
    criterion = criterion,
    wts = wts,
    return_gmm_inv_weights = ifelse(criterion == "GMM",TRUE,FALSE))
  
  if(criterion == "GMM"){
    curr_gmm_inv_wts <- crit_value$inv_wts
    crit_value <- crit_value$gmm_crit
  } else{
    curr_gmm_inv_wts <- NULL
  }
  
  
  # print(crit_value)
  
  
  
  
  newton_counter <- 1
  barrier_counter <- 1
  if(criterion == "GMM"){
    barrier_t <- max_barrier
    newton_counter <- length(tolerances)
    
  }
  while((barrier_t <= max_barrier)&(barrier_counter <= barrier_maxit)){
    derivs <- list()
    derivs$grad <- Inf
    
    while(sqrt(sum(derivs$grad^2)) > tolerances[newton_counter]){
      
      derivs <-
        deriv_criterion_lr(W = W,
                           X = X,
                           Z = Z,
                           which_k_p = which_k_p,
                           which_k_p_tilde = which_k_p_tilde,
                           fixed_P_multipliers = fixed_P_multipliers,
                           fixed_P_tilde_multipliers = fixed_P_tilde_multipliers,
                           which_B_keep = which_B_keep,
                           which_B_rows = which_B_rows,
                           which_gammas = which_gammas,
                           which_gamma_tilde = which_gamma_tilde,
                           Z_tilde = Z_tilde,
                           Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                           alpha_tilde = alpha_tilde,
                           Z_tilde_list = Z_tilde_list,
                           X_tilde = X_tilde,
                           Ak_list = Ak_list,
                           A_tilde_k_list = A_tilde_k_list,
                           fixed_df = fixed_df,
                           varying_df = varying_df,
                           varying_lr_df = varying_lr_df,
                           K = K,
                           K_tilde = K_tilde,
                           barrier_t = barrier_t,
                           criterion = criterion,
                           lr_scale = TRUE,
                           include_log_penalty_derivatives = TRUE,
                           return_info = TRUE,
                           wts = wts,
                           gmm_inv_wts = curr_gmm_inv_wts)
      
      ########### inline testing (bad david)############
      # func_eval <- function(x){
      #   temp_lr <- varying_lr_df
      #   temp_lr$value <- x
      #   temp_crit <- evaluate_criterion_lr(
      #     W = W,
      #     X = X,
      #     Z = Z,
      #     Z_tilde = Z_tilde,
      #     X_tilde = X_tilde,
      #     Z_tilde_gamma_cols = Z_tilde_gamma_cols,
      #     Z_tilde_list = Z_tilde_list,
      #     fixed_df = fixed_df,
      #     varying_df = varying_df,
      #     varying_lr_df = temp_lr,
      #     barrier_t = barrier_t,
      #     criterion = criterion,
      #     wts = wts,
      #     include_log_penalty = TRUE,
      #     gmm_inv_wts = curr_gmm_inv_wts)
      # }
      #
      # numgrad <- numDeriv %in% grad(func_eval, varying_lr_df$value) ## Amy removed double colon
      # #
      # plot(asinh(numgrad), asinh(derivs$grad),pch = 4)
      # # abline(a = 0, b = 1, lty = 2,col = "grey")
      # plot(asinh(derivs$grad - numgrad))
      # # abline(v = c(1,4,9,12,21,30,34), lty = 2)
      # # varying_lr_df[c(1:4,9:12,21,30:34),]
      
      ############# end inline testing ################
      
      step_direction <- NULL
      unconstrained_regularization <-
        hessian_regularization*sqrt(sum(derivs$grad[which_unconstrained]^2))
      
      rho_regularization <-
        hessian_regularization*
        sqrt(sum(derivs$grad[which_rho|which_rho_tilde]^2))
      
      if(verbose){message("calculating step direction")}
      
      
      step_direction <- qr.solve(derivs$info  +
                                   unconstrained_regularization*
                                   diag(as.numeric(which_unconstrained)) +
                                   rho_regularization*
                                   diag(as.numeric(which_rho|which_rho_tilde)),
                                 derivs$grad,
                                 tol = 1e-20)
      
      grad <- derivs$grad
      
      stepsize <- 1
      prop_crit_value <- crit_value + 1
      c_1 <- 1e-4
      if(is.matrix(step_direction)){
        grad_amount <-
          as.numeric(matrix(-step_direction,nrow = 1)%*%derivs$grad)
      } else{
        grad_amount <-
          as.numeric(-matrix(as.numeric(derivs$grad),nrow = 1)%*%
                       step_direction)
        step_direction <- as.numeric(step_direction)
      }
      
      #######
      # crits <- numeric(0)
      # stepsizes <- numeric(0)
      
      # step_direction[28:66] <- 0 # only step in B
      # step_direction[c(1:27,58:66)] <- 0 #only step in gamma
      # step_direction[1:57] <- 0 # only step in rho
      # step_direction <- -step_direction #sign error?
      ######
      
      while((prop_crit_value > crit_value + stepsize*c_1*grad_amount)){
        prop_lr_df <- varying_lr_df
        prop_lr_df$value <- prop_lr_df$value - stepsize*step_direction
        prop_crit_value <- evaluate_criterion_lr(
          W = W,
          X = X,
          Z = Z,
          Z_tilde = Z_tilde,
          X_tilde = X_tilde,
          Z_tilde_gamma_cols = Z_tilde_gamma_cols,
          Z_tilde_list = Z_tilde_list,
          fixed_df = fixed_df,
          varying_df = varying_df,
          varying_lr_df = prop_lr_df,
          barrier_t = barrier_t,
          criterion = criterion,
          wts = wts,
          gmm_inv_wts = curr_gmm_inv_wts)
        
        if(is.nan(prop_crit_value)){
          prop_crit_value <- Inf
        }
        
        if(!is.null(bootstrap_failure_cutoff)){
          if(prop_crit_value < bootstrap_failure_cutoff){
            stop("Bootstrap failed on this iteration")
          }
        }
        
        ################################
        # crits <- c(crits, prop_crit_value)
        # stepsizes <- c(stepsizes,stepsize)
        ################################
        
        stepsize <- stepsize/2
        # print(stepsize)
        
      }
      
      if(stepsize<1e-5){
        print("small step size; increasing hessian regularization")
        hessian_regularization <- hessian_regularization*2
      } else{
        hessian_regularization <- max(hessian_regularization/2,
                                      min_regularization)
      }
      
      
      ############################
      # plot(log(stepsizes), crits)
      ############################
      
      
      varying_lr_df <- prop_lr_df
      crit_value <- prop_crit_value
      barrier_counter <- barrier_counter + 1
      
      if(verbose){
        print(crit_value)
        print(sqrt(sum(derivs$grad^2)))}
      
    }
    
    if(verbose){message(paste("Fit barrier sub-problem with t = ",
                              barrier_t,".", sep = "", collapse = ""))}
    
    barrier_t <- barrier_t*barrier_scale
    newton_counter <- newton_counter + 1
    
  }
  
  varying_df <- lr_to_ra(fixed_df,
                         varying_lr_df,
                         varying_df)
  
  temp_params <- dataframes_to_parameters(fixed_df,
                                          varying_df)
  
  temp_P_fixed_indices <- P_fixed_indices
  temp_P_fixed_indices[] <- TRUE
  temp_P_tilde_fixed_indices <- P_tilde_fixed_indices
  temp_P_tilde_fixed_indices[] <- TRUE
  temp_B_fixed_indices <- B_fixed_indices
  temp_B_fixed_indices[] <- TRUE
  temp_gammas_fixed_indices <- gammas_fixed_indices
  temp_gammas_fixed_indices[] <- TRUE
  temp_gamma_tilde_fixed_indices <- gamma_tilde_fixed_indices
  temp_gamma_tilde_fixed_indices[] <- TRUE
  
  if(!is.null(alpha_tilde)){
    Z_tilde <- construct_Z_tilde(Z_tilde_list,
                                 varying_lr_df$value[
                                   varying_lr_df$param == "alpha_tilde"])
    Z_tilde_list_archived <- Z_tilde_list
    alpha_tilde_archived <-     varying_lr_df$value[
      varying_lr_df$param == "alpha_tilde"]
    Z_tilde_list <- NULL
    alpha_tilde <- NULL
    archived_varying <- varying_df[varying_df$param == "alpha_tilde",]
    varying_df <- varying_df[varying_df$param != "alpha_tilde",]
    temp_params$alpha_tilde <- NULL
  } else{
    archived_varying <- NULL
  }
  
  if(profile_P){
    
    
    
    
    
    P_grad <- numeric(sum(varying_df$param=="P"))
    
    temp_dfs <- parameters_to_dataframes(
      P = temp_params$P,
      P_fixed_indices = temp_P_fixed_indices,
      P_tilde = temp_params$P_tilde,
      P_tilde_fixed_indices = temp_P_tilde_fixed_indices,
      B = temp_params$B,
      B_fixed_indices = temp_B_fixed_indices,
      gammas = temp_params$gammas,
      gammas_fixed_indices = temp_gammas_fixed_indices,
      gamma_tilde = temp_params$gamma_tilde,
      gamma_tilde_fixed_indices = temp_gamma_tilde_fixed_indices,
      alpha_tilde = alpha_tilde)
    crit_value <- evaluate_criterion_lr(W = W,
                                        X = X,
                                        Z = Z,
                                        Z_tilde = Z_tilde,
                                        Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                        Z_tilde_list = Z_tilde_list,
                                        X_tilde = X_tilde,
                                        fixed_df = temp_dfs$fixed,
                                        varying_df = temp_dfs$varying,
                                        varying_lr_df = NULL,
                                        barrier_t = NULL,
                                        criterion = "Poisson",
                                        lr_scale = FALSE,
                                        include_log_penalty = FALSE,
                                        wts = wts)
    
    
    old_crit_value <- crit_value + 10
    P_grad[] <- 100
    
    # if(verbose){
    #   print("About to profile")
    # }
    profile_counter <- 1
    if(length(which_k_p)>0){
      while((old_crit_value - crit_value > 1e-4)&(profile_counter <= profiling_maxit)){
        
        
        for(k in which_k_p){
          temp_P_fixed_indices[k,] <- FALSE
          temp_dfs <- parameters_to_dataframes(
            P = temp_params$P,
            P_fixed_indices = temp_P_fixed_indices,
            P_tilde = temp_params$P_tilde,
            P_tilde_fixed_indices = temp_P_tilde_fixed_indices,
            B = temp_params$B,
            B_fixed_indices = temp_B_fixed_indices,
            gammas = temp_params$gammas,
            gammas_fixed_indices = temp_gammas_fixed_indices,
            gamma_tilde = temp_params$gamma_tilde,
            gamma_tilde_fixed_indices = temp_gamma_tilde_fixed_indices)
          
          jacobian <- mean_jac(varying_df = temp_dfs$varying,
                               fixed_df = temp_dfs$fixed,
                               X = X,
                               Z = Z,
                               X_tilde = X_tilde,
                               Z_tilde = Z_tilde,
                               Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                               Z_tilde_list = Z_tilde_list)
          
          means <- meaninate(gammas = temp_params$gammas,
                             B = temp_params$B,
                             X = X,
                             Z = Z,
                             P = temp_params$P,
                             X_tilde = X_tilde,
                             Z_tilde = Z_tilde,
                             Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                             P_tilde = temp_params$P_tilde,
                             gamma_tilde = temp_params$gamma_tilde,
                             alpha_tilde = temp_params$alpha_tilde,
                             Z_tilde_list = Z_tilde_list)
          
          W_long <- lapply(1:n,function(i) as.numeric(W[i,]))
          W_long <- do.call(c,W_long)
          means_long <- lapply(1:n, function(i) as.numeric(means[i,]))
          means_long <- do.call(c,means_long)
          V <- diag(1/means_long)
          if(!is.null(wts)){
            diag(V) <- diag(V)*wts
          }
          if(!is.null(curr_gmm_inv_wts)){
            diag(V) <- diag(V)/curr_gmm_inv_wts
          }
          # V <- as(V, "sparseMatrix")
          
          p_grad <- -t(jacobian)%*%V%*%(W_long - means_long)
          p_hess <-   t(jacobian)%*%V%*%jacobian +
            diag(rep(hessian_regularization*sqrt(sum(  p_grad^2)),nrow(p_grad)))
          
          P_grad[1:J + (profile_counter - 1)*J] <- p_grad
          profile_counter <- profile_counter + 1
          temp_fn <- function(x){
            temp_dfs$varying$value <- x
            
            
            
            return(evaluate_criterion_lr(W = W,
                                         X = X,
                                         Z = Z,
                                         Z_tilde = Z_tilde,
                                         Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                         Z_tilde_list = Z_tilde_list,
                                         X_tilde = X_tilde,
                                         fixed_df = temp_dfs$fixed,
                                         varying_df = temp_dfs$varying,
                                         varying_lr_df = NULL,
                                         barrier_t = NULL,
                                         criterion = criterion,
                                         lr_scale = FALSE,
                                         include_log_penalty = FALSE,
                                         wts = wts,
                                         gmm_inv_wts = curr_gmm_inv_wts))
          }
          
          temp_dfs$varying$value <-
            simpl_opt_linesearch_fnnls(fn = temp_fn, #objective function to be minimized
                                       x = temp_dfs$varying$value, #starting values of simplex-constrained parameter
                                       xhess = p_hess, # hessian matrix of objective function at x
                                       xgrad = p_grad, # gradient of objective function at x
                                       lambda = 0, #trust penalty
                                       maxit = 100,
                                       constraint_tolerance = constraint_tolerance)
          
          
          ###
          
          temp_crit_value <- temp_fn(  temp_dfs$varying$value)
          # print(temp_crit_value)
          
          
          temp_params <- dataframes_to_parameters(fixed_df = temp_dfs$fixed,
                                                  varying_df = temp_dfs$varying)
          
          
          
          
          
          temp_P_fixed_indices[k,] <- TRUE
          
          # print(k)
          
        }
        old_crit_value <- crit_value
        crit_value <- temp_crit_value
      }
      
    }
    
    # if(verbose){
    #   print("Done profiling")
    # }
  }
  
  final_dfs <-
    parameters_to_dataframes(P = temp_params$P,
                             P_fixed_indices = P_fixed_indices,
                             P_tilde = temp_params$P_tilde,
                             P_tilde_fixed_indices = P_tilde_fixed_indices,
                             B = temp_params$B,
                             B_fixed_indices = B_fixed_indices,
                             gammas = temp_params$gammas,
                             gammas_fixed_indices = gammas_fixed_indices,
                             gamma_tilde = temp_params$gamma_tilde,
                             gamma_tilde_fixed_indices =
                               gamma_tilde_fixed_indices,
                             alpha_tilde  = temp_params$alpha_tilde)
  
  
  crit_value <- evaluate_criterion_lr(
    W = W,
    X = X,
    Z = Z,
    Z_tilde = Z_tilde,
    Z_tilde_gamma_cols = Z_tilde_gamma_cols,
    Z_tilde_list = Z_tilde_list,
    X_tilde = X_tilde,
    fixed_df = final_dfs$fixed,
    varying_df = final_dfs$varying,
    varying_lr_df = NULL,
    barrier_t = NULL,
    criterion = "Poisson",
    lr_scale = FALSE,
    include_log_penalty = FALSE,
    wts = wts)
  
  if(!is.null(archived_varying)){
    final_dfs$varying <- rbind(final_dfs$varying,
                               archived_varying)
    Z_tilde_list <- Z_tilde_list_archived
    alpha_tilde <- alpha_tilde_archived
  }
  
  if(re_pois){
    criterion <- "reweighted_Poisson"
  }
  
  if(!is.null(Z_tilde_list)){
    Z_tilde <- NULL
  }
  
  if(!return_variance){
    variance_df <- NULL
  }
  if(criterion != "reweighted_Poisson"){
    variance_df <- NULL
  }
  
  if(!profile_P){
    profile_counter <- 0
  }
  
  return(list(optimization_status = ifelse(
    (profile_counter < profiling_maxit)&
      (barrier_counter <  barrier_maxit),
    "Converged",
    ifelse((profile_counter >= profiling_maxit)&(barrier_counter>= barrier_maxit),
           "Barrier and profiling stages terminated early; maxit reached for both",
           ifelse(profile_counter >= profiling_maxit,
                  "Profiling stage terminated early; maxit reached",
                  "Barrier stage terminated early; maxit reached"))),
    objective = crit_value,
    varying = final_dfs$varying,
    fixed = final_dfs$fixed,
    W = W,
    X = X,
    Z = Z,
    Z_tilde = Z_tilde,
    Z_tilde_gamma_cols = Z_tilde_gamma_cols,
    gammas = temp_params$gammas,
    gammas_fixed_indices = gammas_fixed_indices,
    P = temp_params$P,
    P_fixed_indices = P_fixed_indices,
    B = temp_params$B,
    B_fixed_indices = B_fixed_indices,
    X_tilde = X_tilde,
    P_tilde = temp_params$P_tilde,
    P_tilde_fixed_indices = P_tilde_fixed_indices,
    gamma_tilde = temp_params$gamma_tilde,
    gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
    alpha_tilde = alpha_tilde,
    Z_tilde_list = Z_tilde_list,
    criterion = criterion,
    weights = wts,
    variance_function = variance_df))
}
