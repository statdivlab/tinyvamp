#################### Barrier Penalty Function #######################
#' Calculate barrier penalty to add to objective function inside
#' barrier algorithm
#'
#' @param varying_lr_df A data frame containing values of parameters that
#' are treated as unknown, with relative abundance parameters represented
#' on the log ratio scale (i.e., as phi and phi_tilde)
#' @param fixed_df A data frame containing values of parameters that are
#' treated as known
#' @param barrier_t The current value of t, the barrier penalty parameter
#'
#' @return The calculated value of the barrier penalty
calculate_log_penalty <- function(varying_lr_df,
                                  fixed_df,
                                  barrier_t){

  which_rho_k <- unique(varying_lr_df$k[varying_lr_df$param=="rho"])
  which_rho_tilde_k <- unique(varying_lr_df$k[varying_lr_df$param=="rho_tilde"])

  log_P <- lapply(which_rho_k,
                  function(k) (varying_lr_df$value[varying_lr_df$param == "rho" &
                                                     varying_lr_df$k == k]) %>%
                    (function(x) c(x,0) - logsum::sum_of_logs(c(x,0))))
  log_P_tilde <- lapply(which_rho_tilde_k,
                        function(k)
                          (varying_lr_df$value[
                            varying_lr_df$param == "rho_tilde" &
                              varying_lr_df$k == k]) %>%
                          (function(x) c(x,0) - logsum::sum_of_logs(c(x,0))))
  log_ra_sum <- do.call(sum,log_P) + do.call(sum,log_P_tilde)

  return((-1/barrier_t)*log_ra_sum)

}
