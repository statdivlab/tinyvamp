##################### Mean Function #####################
meaninate <- function(gammas,
                      B,
                      X,
                      Z,
                      P,
                      X_tilde,
                      Z_tilde = NULL,
                      Z_tilde_gamma_cols,
                      P_tilde,
                      gamma_tilde,
                      alpha_tilde = NULL,
                      Z_tilde_list = NULL,
                      return_separate = FALSE,
                      exclude_gammas = FALSE){

  if(!is.null(alpha_tilde)){
    Z_tilde <- construct_Z_tilde(Z_tilde_list,
                                 alpha_tilde)
  }

  J <- ncol(B)
  n <- nrow(X)

  if(!exclude_gammas){
  #multiply appropriate columns of Z_tilde by exp(gamma)
  if(length(Z_tilde_gamma_cols >0)){
    for(colnum in Z_tilde_gamma_cols){
      Z_tilde[,colnum] <- exp(gammas)* Z_tilde[,colnum]
    }
  }

  sample_part <-
    (Z%*%P)*(
      exp(gammas%*%matrix(1,nrow = 1, ncol = J)
          + X%*%B))
  spurious_part <- Z_tilde%*%(P_tilde *
                                exp(gamma_tilde %*% matrix(1,
                                                           nrow = 1,
                                                           ncol = J) +
                                      X_tilde%*%B))
  } else{


    sample_part <-
      (Z%*%P)*(
        exp(X%*%B))
    spurious_part <- Z_tilde%*%(P_tilde *
                                  exp(gamma_tilde %*% matrix(1,
                                                             nrow = 1,
                                                             ncol = J) +
                                        X_tilde%*%B))

  }

  if(return_separate){
    return(list("sample" = sample_part,
                "spurious" = spurious_part))
  } else{
    return(sample_part + spurious_part)
  }


}
