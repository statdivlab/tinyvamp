log_penalty_hess <- function(rho,barrier_t){
  rho <- as.numeric(rho)
  n_rho <- length(rho)
  log_bt <- log(barrier_t)

  ras <- exp(rho - sum_of_logs(c(rho,0)))
  ras_bt <- exp(rho - log_bt - sum_of_logs(c(rho,0)))

  return((n_rho + 1)*(diag(ras_bt) -
                        matrix(ras_bt,ncol = 1)%*%
                        matrix(ras,nrow = 1)))


}
