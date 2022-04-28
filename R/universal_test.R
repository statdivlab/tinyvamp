

universal_test <- function(W,
                           full_model,
                           null_model,
                           parallelize = FALSE){
  U <- crossfit_U(W = W,
                  full_model = full_model,
                  null_model = null_model,
                  parallelize = parallelize)

  return(1/U)
}
