#################### Barrier Penalty Function #######################

calculate_log_penalty <- function(varying_lr_df,
                                  fixed_df,
                                  barrier_t){

  which_rho_k <- unique(varying_lr_df$k[varying_lr_df$param=="rho"])
  which_rho_tilde_k <- unique(varying_lr_df$k[varying_lr_df$param=="rho_tilde"])

  log_P <- lapply(which_rho_k,
                  function(k) (varying_lr_df$value[varying_lr_df$param == "rho" &
                                                     varying_lr_df$k == k]) %>%
                    (function(x) -c(x,0) + sum_of_logs(c(x,0))))
  log_P_tilde <- lapply(which_rho_tilde_k,
                        function(k)
                          (varying_lr_df$value[
                            varying_lr_df$param == "rho_tilde" &
                              varying_lr_df$k == k]) %>%
                          (function(x) -c(x,0) + sum_of_logs(c(x,0))))
  log_ra_sum <- do.call(sum,log_P) + do.call(sum,log_P_tilde)

  return((1/barrier_t)*log_ra_sum)

}
