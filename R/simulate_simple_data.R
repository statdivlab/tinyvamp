

simulate_simple_data <- function(B,
                                 distrib = "Poisson",
                                 n = 5,
                                 gamma_mean = 11){
  if(nrow(B) != 1){
    stop("B must have a single row")
  }
  J <- ncol(B)
  gammas <- rnorm(n,gamma_mean)

  means <- meaninate(gammas,
            B = B,
            X = matrix(1, nrow = n,ncol = 1),
            Z = matrix(1, nrow = n, ncol = 1),
            P = matrix(1/J,nrow = 1, ncol = J),
            X_tilde = matrix(0,nrow = 1, ncol = 1),
            Z_tilde = matrix(0,nrow = n, ncol = 1),
            Z_tilde_gamma_cols = 1,
            P_tilde = matrix(1/J,nrow = 1, ncol = J),
            gamma_tilde = matrix(0, ncol = 1, nrow = 1),
                        alpha_tilde = NULL,
                        Z_tilde_list = NULL,
                        return_separate = FALSE,
                        exclude_gammas = FALSE)

  if(distrib == "Poisson"){
  W <- apply(means,c(1,2),function(x) rpois(1,x))
  }
  if(distrib == "nb10"){
    W <- apply(means,c(1,2), function(x) rnbinom(1, mu = x, size = 10))
  }

  if(distrib == "nb.5"){
    W <- apply(means,c(1,2), function(x) rnbinom(1, mu = x, size = .5))
  }

return(W)
}
