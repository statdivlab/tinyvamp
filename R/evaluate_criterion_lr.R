#' Criterion Evaluation Function 
#' 
#' @export
evaluate_criterion_lr <- function(W,
                                  X,
                                  Z,
                                  Z_tilde,
                                  Z_tilde_gamma_cols,
                                  Z_tilde_list = NULL,
                                  X_tilde,
                                  fixed_df,
                                  varying_df,
                                  varying_lr_df = NULL,
                                  barrier_t = NULL,
                                  criterion = "Poisson",
                                  lr_scale = TRUE,
                                  include_log_penalty = TRUE,
                                  wts = NULL,
                                  gmm_inv_wts = NULL,
                                  return_gmm_inv_weights = FALSE){

  if(lr_scale){
    varying_df <- lr_to_ra(fixed_df,
                           varying_lr_df,
                           varying_df)
  }

  params <- dataframes_to_parameters(fixed_df, varying_df)

  if(lr_scale){
    if(include_log_penalty){
      log_penalty <- calculate_log_penalty(varying_lr_df,
                                           fixed_df,
                                           barrier_t)
    } else{
      log_penalty <- 0
    }

  } else{
    log_penalty <- 0
  }

  means <- meaninate(gammas  = params$gammas,
                     B = params$B,
                     X = X,
                     Z = Z,
                     P = params$P,
                     X_tilde = X_tilde,
                     Z_tilde = Z_tilde,
                     Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                     Z_tilde_list = Z_tilde_list,
                     P_tilde = params$P_tilde,
                     gamma_tilde = params$gamma_tilde,
                     alpha_tilde = params$alpha_tilde,
                     return_separate = FALSE)



  if(criterion == "Poisson"){
    return(poisson_criterion(W = W,
                             means = means,
                             wts = wts) + log_penalty)

  }
  if(criterion == "GMM"){
    n <- nrow(W)
    W_long <- lapply(1:n,function(i) as.numeric(W[i,]))
    W_long <- do.call(c,W_long)
    means_long <- lapply(1:n, function(i) as.numeric(means[i,]))
    means_long <- do.call(c,means_long)

    if(is.null(gmm_inv_wts)){

    inv_wts <- get_gmm_inv_weights(W_long = W_long,
                           means_long = means_long)
    } else{
      inv_wts <- gmm_inv_wts
    }

    if(!return_gmm_inv_weights){
    return(gmm_criterion(W_long = W_long,
                             means_long = means_long,
                         inv_wts= inv_wts) + log_penalty)
    } else{
      return(list("gmm_crit" = gmm_criterion(W_long = W_long,
                           means_long = means_long,
                           inv_wts = inv_wts) + log_penalty,
             "inv_wts" = inv_wts))
    }

  }


}
