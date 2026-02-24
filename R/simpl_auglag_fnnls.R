#' @import fastnnls
#' 
#' 
simpl_auglag_fnnls <- function(x,
                               fn, #function of x to optimize
                               xhess, #hessian at x
                               xgrad, #gradient at x
                               lambda, #penalty parameters
                               nu = 1, #starting lagrangian penalty
                               mu = 1, #starting augmented lagrangian penalty
                               constraint_tolerance = 1e-10, #sum-to-one constraint tolerance
                               maxit = 100 # maximum number of iterations (outer loop)

){

  npar <- nrow(xhess)

  x0 <- matrix(x, ncol = 1)
  xgrad <- matrix(xgrad,ncol = 1)

  previous_V <- Inf
  # grad_norm <- sqrt(sum(xgrad^2))
  counter <- 0
  for(k in 1:maxit){
    counter <- counter + 1

    ATA <- 0.5*xhess + lambda*diag(rep(1,length(xgrad))) + matrix(mu,nrow = npar, ncol = npar)

    Ab <- t((-0.5)*(matrix(nu - 2*mu,nrow = 1, ncol = npar) -
                      2*lambda*t(x0) -
                      t(x0)%*%xhess + t(xgrad)))
    #
    #     eigen_ATA <- eigen(ATA)
    #
    #     sqrt_eigen_values <- sapply(eigen_ATA$values,
    #                                 function(t) ifelse(t>=0,sqrt(t),0))
    #
    #     inv_sqrt_eigen_values <- sapply(eigen_ATA$values,
    #                                     function(t) ifelse(t>0,1/sqrt(t),0))
    #
    #     A <- eigen_ATA$vectors%*%diag(sqrt_eigen_values)%*%t(eigen_ATA$vectors)
    #     A_inv <- eigen_ATA$vectors%*%diag(inv_sqrt_eigen_values)%*%t(eigen_ATA$vectors)
    #     b <- A_inv%*%Ab
    #
    #
    #     x <- nnls::nnls(A,b)$x
    
    ## BEGIN Amy July 8 2023
    # stopifnot(counter<60)
    if (counter >= 60) {
      warning(paste("In fastNNLS, the counter is", counter, "which is >=60. Potential convergence issue?"))
    }
    ## END Amy July 8 2023
    
    x <- fastnnls::fast_nnls(ZTx = Ab, ZTZ = ATA,
                             tolerance = constraint_tolerance)

    V <- abs(sum(x) - 1)

    satisfied <- V< constraint_tolerance

    if(V < constraint_tolerance){
      return(x)
    }

    if(V < 0.25*previous_V){
      nu <- nu + 2*mu*V
    } else{
      mu <- 2*mu
    }
    previous_V <- V

  }

  warning("Maximum iterations reached; returning initial value")
  return(x0)
}
