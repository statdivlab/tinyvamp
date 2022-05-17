#' @inheritParams dataframes_to_parameters
#' @param X The sample efficiency design -- an \eqn{n \times p} matrix
#' @param Z The sample-specimen design -- an \eqn{n \times K} matrix whose \eqn{ij}-th entry
#' indicates the proportional contribution of specimen \eqn{j} to sample \eqn{i}. Rows must
#' sum to 1 or be identically 0.
#' @param Z_tilde The spurious read design -- an \eqn{n x \tilde{K}} matrix where
#' \eqn{\tilde{K}} is the number of spurious read sources modeled.
#' @param Z_tilde_gamma_cols A numeric vector containing the columns of Z_tilde which should be
#' multiplied by exp(gamma).
#' @param sparse Use sparsity in Jacobian to speed up computation (default is TRUE)
#' 
#' @author David Clausen
mean_jac_lr_faster <- function(fixed_df,
                               varying_lr_df,
                               varying_df,
                               which_k_p,
                               which_k_p_tilde,
                               which_B_rows,
                               which_B_keep,
                               which_gammas,
                               which_gamma_tilde,
                               params,
                               Ak_list,
                               A_tilde_k_list,
                               fixed_P_multipliers,
                               fixed_P_tilde_multipliers,
                               K,
                               K_tilde,
                               X,
                               Z,
                               X_tilde,
                               Z_tilde,
                               Z_tilde_gamma_cols,
                               Z_tilde_list = NULL,
                               sparse = TRUE,
                               proportion_scale = FALSE,
                               P_fixed_indices = NULL,
                               P_tilde_fixed_indices = NULL){

  if(proportion_scale){
    if(is.null(P_fixed_indices)|
       is.null(P_tilde_fixed_indices)){
      stop("P_fixed_indices and P_tilde_fixed_indices must be provided if
proportion_scale = TRUE")
    }
  }

  n <- nrow(X)
  J <- max(c(fixed_df$j,varying_df$j))
  if(proportion_scale){
    K_tilde_effective <-
      length(unique(varying_df$k[varying_df$param == "P_tilde"]))
  } else{
  K_tilde_effective <-
    length(unique(varying_lr_df$k[varying_df$param == "rho_tilde"]))
  }

  if(proportion_scale){
    K_max <- max(
      c(
      unique(
        varying_df$k[varying_df$param == "P"]
      ),
      unique(
        fixed_df$k[fixed_df$param == "P"]
      )
    ))
    K_tilde_max <- max(
      c(
        unique(
          varying_df$k[varying_df$param == "P_tilde"]
        ),
        unique(
          fixed_df$k[fixed_df$param == "P_tilde"]
        )
      ))

  } else{
    K_tilde_effective <-
      length(unique(varying_lr_df$k[varying_df$param == "rho"]))
  }
  if(!proportion_scale){
  K_gamma_tilde_effective <- sum(varying_lr_df$param == "gamma_tilde")
  } else{
    K_gamma_tilde_effective <- sum(varying_df$param == "gamma_tilde")
  }

  if(!proportion_scale){
  n_varying <- nrow(varying_lr_df)
  } else{
    n_varying <- nrow(varying_df)
  }

  #number gammas being estimated
  n_gammas <- sum(varying_lr_df$param == "gamma")

  jacobian <- matrix(0,
                     nrow = n*J,
                     ncol = n_varying)

  if(!proportion_scale){
  npars_rho <- sapply(1:length(Ak_list),
                      function(k) if(is.null(Ak_list[[k]])){
                        0
                      } else{
                        dim(Ak_list[[k]])[2] - 1
                      })

  npars_rho_tilde <- sapply(1:length(A_tilde_k_list),
                            function(k) if(is.null(A_tilde_k_list[[k]])){
                              0
                            } else{
                              dim(A_tilde_k_list[[k]])[2] - 1
                            })
  } else{
    npars_rho <- sapply(1:K_max,
                        function(d) sum(
                          varying_df$k[varying_df$param=="P"] ==d))

    npars_rho_tilde <-  sapply(1:K_tilde_max,
                               function(d) sum(
                                 varying_df$k[varying_df$param=="P_tilde"] ==d))

  }

  ### NB: which_B_keep reflects only rows of B that are not entirely fixed
  if(!is.null(which_B_keep)){
    npars_B <- apply(which_B_keep,1,sum)
  } else{
    npars_B <- 0
  }

  for(i in 1:n){
    # print(i)

    if(sum(npars_B) > 0){
      ### add derivatives wrt elements of B
      which_B_local <- which((X[i,] != 0 )|
                               apply(X_tilde[Z_tilde[i,] > 0,,drop = F],
                                     2,
                                     function(xtil) as.logical(sum(xtil!=0))))

      #only compute over rows of B that are not fixed
      which_B_local <- which_B_local[which_B_local %in% which_B_rows]

      #index which_B_local in terms of rows that are not fixed
      which_B_local <- sapply(which_B_local,
                              function(x) sum(which_B_rows<=x))

    if(length(which_B_local) > 0){
      for(B_row_index in which_B_local){

        jacobian[(i - 1)*J + 1:J,
                 sum(npars_B[1:nrow(params$B) < B_row_index]) +
                   1:(npars_B[B_row_index])] <-
          diag(sapply(1:J,
                      function(j)
                        ifelse(j<J,
                               ifelse(which_B_keep[B_row_index,j],
                                      mu_d_beta(i = i,
                                                j = j,
                                                q = B_row_index,
                                                gammas = params$gammas,
                                                B = params$B,
                                                X = X,
                                                Z = Z,
                                                P = params$P,
                                                X_tilde = X_tilde,
                                                Z_tilde = Z_tilde,
                                                Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                                Z_tilde_list = Z_tilde_list,
                                                alpha_tilde = params$alpha_tilde,
                                                P_tilde = params$P_tilde,
                                                gamma_tilde = params$gamma_tilde),0),
                               0)))[,which(which_B_keep[B_row_index,])]


        ############## inline testing ################
        # i <- 1
        # j <- 1
        # B_row_index <- 1
        # an_grad <- mu_d_beta(i = i,
        #           j = j,
        #           q = B_row_index,
        #           gammas = params$gammas,
        #           B = params$B,
        #           X = X,
        #           Z = Z,
        #           P = params$P,
        #           X_tilde = X_tilde,
        #           Z_tilde = Z_tilde,
        #           Z_tilde_gamma_cols = Z_tilde_gamma_cols,
        #           Z_tilde_list = Z_tilde_list,
        #           alpha_tilde = params$alpha_tilde,
        #           P_tilde = params$P_tilde,
        #           gamma_tilde = params$gamma_tilde)
        #
        # mean_func <- function(x,i,j){
        #   B_temp <- params$B
        #   B_temp[1,1] <- x
        #   means <- meaninate(gammas,
        #                                    B = B_temp,
        #                                    X,
        #                                    Z,
        #                                    P= params$P,
        #                                    X_tilde,
        #                                    Z_tilde = Z_tilde,
        #                                    Z_tilde_gamma_cols,
        #                                    P_tilde = params$P_tilde,
        #                                    gamma_tilde = params$gamma_tilde,
        #                                    alpha_tilde = params$alpha_tilde,
        #                                    Z_tilde_list = Z_tilde_list,
        #                                    return_separate = FALSE,
        #                                    exclude_gammas = FALSE)
        #   return(means[i,j])
        # }
        #
        # num_grad <- numDeriv::grad(function(x) mean_func(x,1,1),params$B[1,1])
        #
        #


        ##############                ################


      }}}


    if(n_gammas>0){
    ### add derivatives wrt gamma_i
    if(i %in% which_gammas){
    jacobian[(i - 1)*J + 1:J,
             sum(npars_B) + sum(which_gammas <= i)] <-
      mu_d_gamma_faster(i,
                        J,
                        gammas = params$gammas,
                        B = params$B,
                        X = X,
                        Z = Z,
                        P = params$P,
                        X_tilde = X_tilde,
                        Z_tilde = Z_tilde,
                        Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                        P_tilde = params$P_tilde,
                        gamma_tilde = params$gamma_tilde,
                        alpha_tilde = params$alpha_tilde,
                        Z_tilde_list = Z_tilde_list)
    }
}
    ### add derivatives wrt gamma_tilde


    if(length(which_k_p_tilde)>0){
      for(k_tilde in which_k_p_tilde){
        jacobian[(i - 1)*J + 1:J,
                 sum(npars_B) + n_gammas + k_tilde] <- sapply(1:J, function(j)
                   mu_d_gamma_tilde(i,
                                    j = j,
                                    k_tilde = k_tilde,
                                    gammas = params$gammas,
                                    B = params$B,
                                    X = X,
                                    Z = Z,
                                    P = params$P,
                                    X_tilde = X_tilde,
                                    Z_tilde = Z_tilde,
                                    Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                    Z_tilde_list = Z_tilde_list,
                                    alpha_tilde = params$alpha_tilde,
                                    P_tilde = params$P_tilde,
                                    gamma_tilde = params$gamma_tilde))
      }}

    # which rho derivs can be nonzero for sample i?
    which_rho_local <- which(Z[i,]>0)
    which_rho_local <- which_rho_local[which_rho_local %in% which_k_p]


    if(length(which_rho_local)>0){
      for(rho_index in which_rho_local){
        if(!proportion_scale){
        jacobian[(i - 1)*J + 1:J,
                 sum(npars_B) + n_gammas + K_gamma_tilde_effective +
                   sum(npars_rho[which_k_p[which_k_p < rho_index]]) + 1:(npars_rho[rho_index])] <-
          as.matrix(mu_d_rho_faster(i,
                                    J,
                                    gammas = params$gammas,
                                    B = params$B,
                                    X = X,
                                    Z = Z,
                                    k = rho_index,
                                    Ak_list = Ak_list,
                                    rho_k = varying_lr_df$value[varying_lr_df$k == rho_index &
                                                                  varying_lr_df$param == "rho"],
                                    fixed_P_multipliers = fixed_P_multipliers
          ))

        ########### testing #############
        # par(mfrow = c(4,4))
        # for(i in 9:16){
        # analgrad <- as.matrix(mu_d_rho_faster(i,
        #                                       J,
        #                                       gammas = params$gammas,
        #                                       B = params$B,
        #                                       X = X,
        #                                       Z = Z,
        #                                       k = rho_index,
        #                                       Ak_list = Ak_list,
        #                                       rho_k = varying_lr_df$value[varying_lr_df$k == rho_index &
        #                                                                     varying_lr_df$param == "rho"],
        #                                       fixed_P_multipliers = fixed_P_multipliers
        # ))
        #
        # mean_func <- function(x,index,i,j){
        #   temp_lr <- varying_lr_df
        #   temp_lr$value[index] <- x
        #   temp_params <- dataframes_to_parameters(fixed_df,
        #                                           lr_to_ra(fixed_df,
        #                                                    temp_lr,
        #                                                    varying_df))
        #   return(with(temp_params,meaninate(gammas,
        #             B,
        #             X,
        #             Z,
        #             P,
        #             X_tilde,
        #             Z_tilde = NULL,
        #             Z_tilde_gamma_cols,
        #             P_tilde,
        #             gamma_tilde,
        #             alpha_tilde = alpha_tilde,
        #             Z_tilde_list = Z_tilde_list,
        #             return_separate = FALSE,
        #             exclude_gammas = FALSE)[i,j]))
        # }
        # numgrad <- 0*analgrad
        # for(muj in 1:5){
        #   for(rhoj in 1:4){
        #     index <- rhoj + 21
        #     numgrad[muj,rhoj] <- numDeriv::grad(function(x) mean_func(x,index,i = i,muj),varying_lr_df$value[index])
        #   }
        # }
        #
        # plot(numgrad, analgrad,main = i)
        # abline( a= 0, b= 1, lty = 2)
        #
        # plot(numgrad - analgrad, main = i)
        # abline(a = 0, b=0, lty = 2)
        # abline(v = 0, lty = 2)
        # }
        #
        #
        #
        # ############## end testing #################
        } else{
          jacobian[(i - 1)*J + 1:J,
                   sum(npars_B) + n_gammas + K_gamma_tilde_effective +
                     sum(npars_rho[which_k_p[which_k_p < rho_index]]) + 1:(npars_rho[rho_index])] <-
            as.matrix(mu_d_rho_faster(i,
                                      J,
                                      gammas = params$gammas,
                                      B = params$B,
                                      X = X,
                                      Z = Z,
                                      k = rho_index,
                                      Ak_list = Ak_list,
                                      rho_k = varying_lr_df$value[varying_lr_df$k == rho_index &
                                                                    varying_lr_df$param == "rho"],
                                      fixed_P_multipliers = fixed_P_multipliers,
                                      proportion_scale = proportion_scale
            ))[,!P_fixed_indices[rho_index,]]

        }


      }
    }


    if(is.null(Z_tilde_list)){
    which_rho_tilde_local <- which(Z_tilde[i,]>0)
    } else{
      which_rho_tilde_local <-
        unique(
          do.call(c,lapply(Z_tilde_list, function(Ztil) which(Ztil[i,]>0))))
    }
    which_rho_tilde_local <-     which_rho_tilde_local[
      which_rho_tilde_local %in% which_k_p_tilde]

    if(length(which_rho_tilde_local) >0){
      for(rho_tilde_index in which_rho_tilde_local){
        if(!proportion_scale){

        jacobian[(i - 1)*J + 1:J,
                 sum(npars_B) + n_gammas + K_gamma_tilde_effective +
                   sum(npars_rho) + ifelse(rho_tilde_index ==1,0, sum(sapply(1:(rho_tilde_index-1),
                                               function(w) npars_rho_tilde[[w]]))) +
                   1:npars_rho_tilde[[rho_tilde_index]]] <-
          as.matrix(mu_d_rho_tilde_faster(i = i,
                                          J = J,
                                          gammas = params$gammas,
                                          B = params$B,
                                          rho_tilde_k =
                                            varying_lr_df$value[
                                              varying_lr_df$k == rho_tilde_index &
                                                varying_lr_df$param == "rho_tilde"],
                                          k_tilde = rho_tilde_index,
                                          A_tilde_k_list = A_tilde_k_list,
                                          fixed_P_tilde_multipliers = fixed_P_tilde_multipliers,
                                          X_tilde = X_tilde,
                                          Z_tilde = Z_tilde,
                                          Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                          gamma_tilde = params$gamma_tilde,
                                          alpha_tilde = params$alpha_tilde,
                                          Z_tilde_list = Z_tilde_list))
        } else{

          jacobian[(i - 1)*J + 1:J,
                   sum(npars_B) + n_gammas + K_gamma_tilde_effective +
                     sum(npars_rho) + ifelse(rho_tilde_index ==1,0, sum(sapply(1:(rho_tilde_index-1),
                                                                               function(w) npars_rho_tilde[[w]]))) +
                     1:npars_rho_tilde[[rho_tilde_index]]] <-
            as.matrix(
              mu_d_rho_tilde_faster(i = i,
                                    J = J,
                                    gammas = params$gammas,
                                    B = params$B,
                                    rho_tilde_k = NULL,
                                    k_tilde = rho_tilde_index,
                                    A_tilde_k_list = NULL,
                                    fixed_P_tilde_multipliers = NULL,
                                    X_tilde = X_tilde,
                                    Z_tilde = Z_tilde,
                                    Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                    gamma_tilde = params$gamma_tilde,
                                    alpha_tilde = params$alpha_tilde,
                                    Z_tilde_list = Z_tilde_list,
                                    proportion_scale = TRUE)
            )[,!P_tilde_fixed_indices[rho_tilde_index,]]

        }
      }
    }

    if(sum(varying_lr_df$param == "alpha_tilde")>0){
      for(a_tilde in 1:length(params$alpha_tilde)){
      jacobian[(i - 1)*J + 1:J,
               sum(npars_B) + n_gammas + K_gamma_tilde_effective +
                 sum(npars_rho) + sum(npars_rho_tilde) + a_tilde] <-
        mu_d_alpha_tilde(i,
                         J,
                         a_tilde = a_tilde,
                         gammas = params$gammas,
                         B = params$B,
                         X_tilde = X_tilde,
                         Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                         alpha_tilde = params$alpha_tilde,
                         Z_tilde_list = Z_tilde_list,
                         P_tilde = params$P_tilde,
                         gamma_tilde = params$gamma_tilde)

      }
    }





  }

  if(sparse){
    jacobian <- as(jacobian,"sparseMatrix")
  }


  return(jacobian)
}
