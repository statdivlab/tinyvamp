
log_penalty_grad <- function(rho,barrier_t){
  # d/d_rho sum(-log(exp(c(rho, 0))/sum(exp(rho,0)))
  # = d/d_rho sum( - c(rho,0) + log(sum(exp(rho,0))))
  rho <- as.numeric(rho)
  log_bt <- log(barrier_t)

  n_rho <- length(rho)

  return(rep(-1/barrier_t,n_rho) + (n_rho +1)*exp(rho - log_bt - sum_of_logs(c(rho,0))))
}
