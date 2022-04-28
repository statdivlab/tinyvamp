
ra_to_lr <- function(varying_df){

  P_df <- varying_df[varying_df$param == "P",]
  P_tilde_df <- varying_df[varying_df$param == "P_tilde",]

  ### Generate rho data.frame
  P_rows <- unique(P_df$k)

  rho_df <- P_df[rep(F, nrow(P_df)),]
  for(p_row in P_rows){
    mini_P <- P_df[P_df$k == p_row,]
    ref_j <- max(mini_P$j)
    mini_rho <- mini_P[mini_P$j != ref_j,]
    mini_rho$param <- "rho"
    mini_rho$value <- log(mini_rho$value) - log(mini_P$value[mini_P$j == ref_j])
    rho_df <- rbind(rho_df, mini_rho)
  }

  ### Generate rho_tilde data.frame
  P_tilde_rows <- unique(P_tilde_df$k)

  rho_tilde_df <- P_tilde_df[rep(F, nrow(P_df)),]
  for(p_tilde_row in P_tilde_rows){
    mini_P_tilde <- P_tilde_df[P_tilde_df$k == p_tilde_row,]
    ref_j <- max(mini_P_tilde$j)
    mini_rho_tilde <- mini_P_tilde[mini_P_tilde$j != ref_j,]
    mini_rho_tilde$param <- "rho_tilde"
    mini_rho_tilde$value <- log(mini_rho_tilde$value) -
      log(mini_P_tilde$value[mini_P_tilde$j == ref_j])
    rho_tilde_df <- rbind(rho_tilde_df, mini_rho_tilde)
  }

  varying_lr_df <- rbind(varying_df[!(varying_df$param %in% c("P","P_tilde")),],
                         rho_df,
                         rho_tilde_df)

  # varying_lr_df <- rbind(varying_lr_df[varying_lr_df$param == "rho",],
  #                     varying_lr_df[varying_lr_df$param == "rho_tilde",],
  #                     varying_lr_df[varying_lr_df$param == "B",],
  #                     varying_lr_df[varying_lr_df$param == "gamma",],
  #                     varying_lr_df[varying_lr_df$param == "gamma_tilde",])

  rownames(varying_lr_df) <- 1:nrow(varying_lr_df)

  varying_lr_df <- rbind(varying_lr_df[varying_lr_df$param != "alpha_tilde",],
                         varying_lr_df[varying_lr_df$param == "alpha_tilde",])

  return(varying_lr_df)

}
