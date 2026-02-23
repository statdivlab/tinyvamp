#' Convert parameter values stored in data frame format to matrix format
#'
#' @param fixed_df A dataframe containing values of model parameters
#' treated as fixed and known (i.e. held constant at known values)
#' 
#' @param varying_df A dataframe containing current values of model parameters treated
#' as fixed and unknown (i.e., parameters to be estimated)
#' 
#' @return A list containing
#' \item{P}{Specimen relative abundance matrix (of dimension K x J)}
#' \item{P_tilde}{Spurious read source
#' relative abundance matrix (of dimension K-tilde x J)}
#' \item{B}{A matrix of detection efficiencies (of dimension p x J)}
#' \item{gammas}{An n-vector of sample-specific read intensities}
#' \item{gamma_tilde}{A-vector of spurious read source intensities
#' (of length K-tilde)}
#' 
#' @author David Clausen
#' 
dataframes_to_parameters <- function(fixed_df,
                                     varying_df){

  together_df <- rbind(fixed_df,varying_df)

  K <- max(together_df$k[together_df$param == "P"])

  K_tilde <- max(together_df$k[together_df$param == "P_tilde"])

  J <- max(together_df$j[together_df$param == "P"])

  p <- max(together_df$k[together_df$param == "B"])

  P <- matrix(ncol = J, nrow = K)

  P_df <- together_df[together_df$param == "P",]
  for(k in 1:K){
    P_row <- P_df[P_df$k ==k,]
    P[k,] <- P_row$value[order(P_row$j)]
  }

  P_tilde <- matrix(ncol = J, nrow = K_tilde)
  P_tilde_df <- together_df[together_df$param == "P_tilde",]
  for(k in 1:K_tilde){
    P_tilde_row <- P_tilde_df[P_tilde_df$k ==k,]
    P_tilde[k,] <- P_tilde_row$value[order(P_tilde_row$j)]
  }

  B_df <- together_df[together_df$param == "B",]
  B <- matrix(0,ncol = J, nrow = p)

  for(k in 1:p){
    B_row <- B_df[B_df$k ==k,]
    B[k,] <- B_row$value[order(B_row$j)]
  }

  B[,J] <- 0

  gammas_df <- together_df[together_df$param == "gamma",]

  gammas <- matrix(nrow = nrow(gammas_df), ncol = 1)

  gammas[] <- gammas_df$value[order(gammas_df$k)]

  gamma_tilde_df <- together_df[together_df$param == "gamma_tilde",]

  gamma_tilde <- matrix(nrow = nrow(gamma_tilde_df), ncol = 1)

  gamma_tilde[] <- gamma_tilde_df$value[order(gamma_tilde_df$k)]

  if(sum(varying_df$param == "alpha_tilde") > 0){
    alpha_tilde <- varying_df$value[varying_df$param == "alpha_tilde"]
  } else{
    alpha_tilde <- NULL
  }

  return(list("P" = P,
              "P_tilde" = P_tilde,
              "B" = B,
              "gammas" = gammas,
              "gamma_tilde" = gamma_tilde,
              "alpha_tilde" = alpha_tilde))
}
