get_Ak_list <- function(fixed_df,
                        varying_df,
                        varying_lr_df){

  J <- max(c(varying_df$j,fixed_df$j))
  K <- max(c(varying_df$k[varying_df$param %in% c("P")],
             fixed_df$k[fixed_df$param %in% c("P")]))

  which_k_p <- unique(varying_lr_df$k[varying_lr_df$param == "rho"])
  which_k_p <- which_k_p[order(which_k_p)]

  Ak_list <- vector(mode = "list", K)
  for(k in which_k_p){
    fixed_p_k <- fixed_df[(fixed_df$param == "P")&(fixed_df$k == k),]
    varying_p_k <- varying_df[(varying_df$param == "P")&
                                (varying_df$k == k),]
    varying_j <- varying_p_k$j
    varying_j <- varying_j[order(varying_j)]
    C_k <- J - nrow(fixed_p_k)
    A_k <- matrix(0, nrow = J,
                  ncol = C_k)
    for(jstar in 1:length(varying_j)){
      A_k[varying_j[jstar],jstar] <- 1
    }
    Ak_list[[k]] <- A_k
  }
  return(Ak_list)
}
