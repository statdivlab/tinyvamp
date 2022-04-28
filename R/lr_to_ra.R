
lr_to_ra <- function(fixed_df,
                     varying_lr_df,
                     varying_df){

  rho_df <- varying_lr_df[varying_lr_df$param == "rho",]
  rho_tilde_df <- varying_lr_df[varying_lr_df$param == "rho_tilde",]

  P_df <- varying_df[varying_df$param == "P",]
  P_tilde_df <- varying_df[varying_df$param == "P_tilde",]



  #update values in P_df with values from varying_lr_df
  rho_rows <- unique(rho_df$k)
  for(rho_row in rho_rows){

    ref_j <- max(varying_df$j[varying_df$param == "P" & varying_df$k == rho_row])

    rho_mini <- rho_df[rho_df$k == rho_row,]

    unscaled_ra <- exp(c(rho_mini$value, 0) - sum_of_logs(c(rho_mini$value, 0)))

    scaling <- 1 - sum(fixed_df$value[(fixed_df$k == rho_row) &(fixed_df$param == "rho")])
    scaled_ra <- unscaled_ra*scaling

    for(j_index in 1:length(c(rho_mini$j, ref_j))){
      P_df$value[
        (P_df$k == rho_row) &
          (P_df$j == c(rho_mini$j, ref_j)[j_index])
      ] <-
        scaled_ra[j_index]

    }

  }

  #update values in P_tilde_df with values from varying_lr_df
  rho_tilde_rows <- unique(rho_tilde_df$k)
  for(rho_tilde_row in rho_tilde_rows){

    ref_j <- max(varying_df$j[varying_df$param == "P_tilde" &
                                varying_df$k == rho_tilde_row])

    rho_tilde_mini <- rho_tilde_df[rho_tilde_df$k == rho_tilde_row,]

    unscaled_ra <- exp(c(rho_tilde_mini$value, 0) - sum_of_logs(c(rho_tilde_mini$value, 0)))

    scaling <- 1 - sum(fixed_df$value[(fixed_df$k == rho_tilde_row) &
                                        (fixed_df$param == "rho_tilde")])
    scaled_ra <- unscaled_ra*scaling

    for(j_index in 1:length(c(rho_tilde_mini$j, ref_j))){
      P_tilde_df$value[
        (P_tilde_df$k == rho_tilde_row) &
          (P_tilde_df$j == c(rho_tilde_mini$j, ref_j)[j_index])
      ] <-
        scaled_ra[j_index]

    }

  }

  P_df <- P_df[order(P_df$k,P_df$j),]
  P_tilde_df <- P_tilde_df[order(P_tilde_df$k,
                                 P_tilde_df$j),]

  to_return_df <-  rbind(varying_lr_df[
    !(varying_lr_df$param %in% c("rho", "rho_tilde","alpha_tilde")),
  ],
  P_df,
  P_tilde_df)

  to_return_df <- rbind(to_return_df,
                        varying_lr_df[varying_lr_df$para == "alpha_tilde",])

  to_return_df <- rbind(to_return_df[to_return_df$param == "P",],
                      to_return_df[to_return_df$param == "P_tilde",],
                      to_return_df[to_return_df$param == "B",],
                      to_return_df[to_return_df$param == "gamma",],
                      to_return_df[to_return_df$param == "gamma_tilde",],
                      to_return_df[to_return_df$param == "alpha_tilde",])

  return(to_return_df)
}
