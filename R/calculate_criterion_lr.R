#' @export
#' 
calculate_criterion_lr <- function(W,
                                   X,
                                   Z,
                                   Z_tilde = NULL,
                                   Z_tilde_gamma_cols,
                                   Z_tilde_list = NULL,
                                   alpha_tilde = NULL,
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
                                   criterion = "Poisson", 
                                   lr_scale = TRUE, 
                                   wts = NULL) {
  
  stopifnot(criterion == "Poisson")
  stopifnot(lr_scale == TRUE)
  
  means <- meaninate(gammas = gammas, 
                     B = B, 
                     X = X,
                     Z = Z,
                     P = P, 
                     X_tilde = X_tilde,
                     Z_tilde = Z_tilde,
                     Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                     Z_tilde_list = Z_tilde_list,
                     P_tilde = P_tilde, 
                     gamma_tilde = gamma_tilde, 
                     alpha_tilde = alpha_tilde, 
                     return_separate = FALSE)
  
  poisson_criterion(W = W,
                    means = means,
                    wts = wts)
  
}