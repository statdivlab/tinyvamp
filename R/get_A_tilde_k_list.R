get_A_tilde_k_list <- function(fixed_df,
                               varying_df,
                               varying_lr_df){

  J <- max(c(varying_df$j,fixed_df$j))
  K_tilde <- max(c(varying_df$k[varying_df$param == "P_tilde"],
                   fixed_df$k[fixed_df$param == "P_tilde"]))

  which_k_p_tilde <- unique(varying_lr_df$k[varying_lr_df$param == "rho_tilde"])
  which_k_p_tilde <- which_k_p_tilde[order(which_k_p_tilde)]

  A_tilde_k_list <- vector(mode = "list", K_tilde)
  for(k in which_k_p_tilde){
    fixed_p_tilde_k <- fixed_df[(fixed_df$param == "P_tilde")&(fixed_df$k == k),]
    varying_p_tilde_k <- varying_df[(varying_df$param == "P_tilde")&
                                      (varying_df$k == k),]
    varying_j <- varying_p_tilde_k$j
    varying_j <- varying_j[order(varying_j)]
    C_k <- J - nrow(fixed_p_tilde_k)
    A_tilde_k <- matrix(0, nrow = J,
                        ncol = C_k)
    for(jstar in 1:length(varying_j)){
      A_tilde_k[varying_j[jstar],jstar] <- 1
    }
    A_tilde_k_list[[k]] <- A_tilde_k
  }
  return(A_tilde_k_list)
}
