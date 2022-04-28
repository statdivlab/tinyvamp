

nullify <- function(W,
                    full_model,
                    null_model){

  J <- ncol(full_model$W)
  n <- nrow(full_model$W)

full_means <-
  meaninate(gammas = full_model$gammas,
            B = full_model$B,
            X = full_model$X,
            Z = full_model$Z,
            P = full_model$P,
            X_tilde = full_model$X_tilde,
            Z_tilde = full_model$Z_tilde,
            Z_tilde_gamma_cols = full_model$Z_tilde_gamma_cols,
            P_tilde = full_model$P_tilde,
            gamma_tilde = full_model$gamma_tilde,
            alpha_tilde = full_model$alpha_tilde,
            Z_tilde_list = full_model$Z_tilde_list,
            exclude_gammas= TRUE)

null_means <-
  meaninate(gammas = null_model$gammas,
            B = null_model$B,
            X = null_model$X,
            Z = null_model$Z,
            P = null_model$P,
            X_tilde = null_model$X_tilde,
            Z_tilde = null_model$Z_tilde,
            Z_tilde_gamma_cols = null_model$Z_tilde_gamma_cols,
            P_tilde = null_model$P_tilde,
            gamma_tilde = null_model$gamma_tilde,
            alpha_tilde = null_model$alpha_tilde,
            Z_tilde_list = null_model$Z_tilde_list,
            exclude_gammas = TRUE)

W0 <- W*(null_means/full_means)
W0[full_means == 0] <- 0

for(i in 1:nrow(W0)){
  W0[i,] <- (W0[i,]/sum(W0[i,]))*sum(W[i,])
}

return(W0)

}
