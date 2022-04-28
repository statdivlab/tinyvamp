parameters_to_dataframes <- function(P,
                                     P_fixed_indices,
                                     P_tilde,
                                     P_tilde_fixed_indices,
                                     B,
                                     B_fixed_indices,
                                     gammas,
                                     gammas_fixed_indices,
                                     gamma_tilde,
                                     gamma_tilde_fixed_indices,
                                     alpha_tilde = NULL){

  ### extract important dimensions
  J <- ncol(P)
  K <- nrow(P)
  K_tilde <- nrow(P_tilde)
  p_B <- nrow(B)
  n <- nrow(as.matrix(gammas,ncol = 1))

  P_fixed_indices <- apply(P_fixed_indices, c(1,2), as.logical)
  P_tilde_fixed_indices <- apply(P_tilde_fixed_indices, c(1,2), as.logical)

  ### check that P_fixed_indices and P satisfy RA requirements
  for(k in 1:K){
    if(sum(!P_fixed_indices[k,]) >0){
      C_k <- sum(P_fixed_indices[k,])
      fixed_sum_k <- sum(P[k,P_fixed_indices[k,]])
      if(C_k >= J - 1){
        P_fixed_indices[k,] <- TRUE
        warning(paste("Row ", k, " of P set to fixed; ",
                      "\nprovided number of known elements via P_fixed_indices is greater ",
                      "\nthan number of taxa - 2, which implies entire row is known."))
      }
      if(fixed_sum_k >=1){
        P_fixed_indices[k,] <- TRUE
        warning(paste("Row ", k, " of P set to fixed; ",
                      "\ntotal relative abundances across taxa indicated as known",
                      "\nvia P_fixed_indices is 1, which implies entire row is known."))
      }
    }
  }

  ### check that P_tilde_fixed_indices and P_tilde satisfy RA requirements
  for(k in 1:K_tilde){
    if(sum(!P_tilde_fixed_indices[k,]) >0){
      C_k <- sum(P_tilde_fixed_indices[k,])
      fixed_sum_k <- sum(P_tilde[k,P_tilde_fixed_indices[k,]])
      if(C_k >= J - 1){
        P_tilde_fixed_indices[k,] <- TRUE
        warning(paste("Row ", k, " of P_tilde set to fixed; ",
                      "\nprovided number of known elements via P_tilde_fixed_indices is greater ",
                      "\nthan number of taxa - 2, which implies entire row is known."))
      }
      if(fixed_sum_k >=1){
        P_tilde_fixed_indices[k,] <- TRUE
        warning(paste("Row ", k, " of P_tilde set to fixed; ",
                      "\ntotal relative abundances across taxa indicated as known",
                      "\nvia P_tilde_fixed_indices is 1, which implies entire row is known."))
      }
    }
  }

  ### set up fixed parameter and varying parameter data.frames
  fixed_df <- data.frame("value" = numeric(0),
                         "param" = character(0),
                         "j" = numeric(0))

  varying_df <- data.frame("value" = numeric(0),
                           "param" = character(0),
                           "j" = numeric(0))


  ### Set up matrix to track specimen provenance
  p_k_matrix <- P
  for(k in 1:K){
    p_k_matrix[k,] <- rep(k,J)
  }

  ### Set up matrix to track taxon provenance

  p_j_matrix <- P
  for(j in 1:J){
    p_j_matrix[,j] <- rep(j,K)
  }

  known_temp_df <- data.frame("value" = P[P_fixed_indices],
                              "param" = rep("P",sum(P_fixed_indices)),
                              "k" = p_k_matrix[P_fixed_indices],
                              "j" = p_j_matrix[P_fixed_indices])

  varying_temp_df <- data.frame("value" = P[!P_fixed_indices],
                                "param" = rep("P",
                                              sum(!P_fixed_indices)),
                                "k" = p_k_matrix[!P_fixed_indices],
                                "j" = p_j_matrix[!P_fixed_indices])

  fixed_df <- rbind(fixed_df,
                    known_temp_df)

  varying_df <- rbind(varying_df,
                      varying_temp_df)

  ### Set up matrix to track spurious read specimen provenance
  p_tilde_k_matrix <- P_tilde
  for(k in 1:K_tilde){
    p_tilde_k_matrix[k,] <- rep(k,J)
  }

  ### Set up matrix to track spurious read taxon provenance

  p_tilde_j_matrix <- P_tilde
  for(j in 1:J){
    p_tilde_j_matrix[,j] <- rep(j,K_tilde)
  }

  known_temp_df <- data.frame("value" = P_tilde[P_tilde_fixed_indices],
                              "param" = rep("P_tilde",sum(P_tilde_fixed_indices)),
                              "k" = p_tilde_k_matrix[P_tilde_fixed_indices],
                              "j" = p_tilde_j_matrix[P_tilde_fixed_indices])

  varying_temp_df <- data.frame("value" = P_tilde[!P_tilde_fixed_indices],
                                "param" = rep("P_tilde",sum(!P_tilde_fixed_indices)),
                                "k" = p_tilde_k_matrix[!P_tilde_fixed_indices],
                                "j" = p_tilde_j_matrix[!P_tilde_fixed_indices])

  fixed_df <- rbind(fixed_df,
                    known_temp_df)

  varying_df <- rbind(varying_df,
                      varying_temp_df)



  ### Set up matrix to track effects
  B_k_matrix <- B
  for(k in 1:p_B){
    B_k_matrix[k,] <- rep(k,J)
  }

  ### Set up matrix to track effects by taxon

  B_j_matrix <- B
  for(j in 1:J){
    B_j_matrix[,j] <- rep(j,p_B)
  }


  known_temp_df <- data.frame("value" = B[B_fixed_indices],
                              "param" = rep("B", sum(B_fixed_indices)),
                              "k" = B_k_matrix[B_fixed_indices],
                              "j" = B_j_matrix[B_fixed_indices])

  varying_temp_df <- data.frame("value" = B[!B_fixed_indices],
                                "param" = rep("B",sum(!B_fixed_indices)),
                                "k" = B_k_matrix[!B_fixed_indices],
                                "j" = B_j_matrix[!B_fixed_indices])

  ### make sure B is ordered by k and *then* j

  varying_temp_df <- varying_temp_df[order(varying_temp_df$k,
                                           varying_temp_df$j),]

  fixed_df <- rbind(fixed_df,
                    known_temp_df)

  varying_df <- rbind(varying_df,
                      varying_temp_df)



  temp_df <- data.frame("value" = as.numeric(gammas),
                        "param" = "gamma",
                        "k" = 1:n,
                        "j" = 0)
  fixed_df <- rbind(fixed_df,temp_df[((1:n) %in%  which(gammas_fixed_indices)),])

  varying_df <- rbind(varying_df,temp_df[!((1:n) %in%  which(gammas_fixed_indices)),])

  temp_df <- data.frame("value" = as.numeric(gamma_tilde),
                        "param" = "gamma_tilde",
                        "k" = 1:nrow(gamma_tilde),
                        "j" = 0)

  K_tilde <- length(as.numeric(gamma_tilde))
  fixed_df <- rbind(fixed_df,temp_df[((1:K_tilde) %in%  which(gamma_tilde_fixed_indices)),])

  varying_df <- rbind(varying_df,temp_df[!((1:K_tilde) %in%  which(gamma_tilde_fixed_indices)),])

  if(!is.null(alpha_tilde)){
    varying_df <- rbind(varying_df,
                        data.frame("value" = alpha_tilde,
                                   "param" = "alpha_tilde",
                                   "k" = 1:length(alpha_tilde),
                                   "j" = 0))
  }

  return(list("fixed" = fixed_df,
              "varying" = varying_df))

}
