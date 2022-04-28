

construct_Z_tilde <- function(Z_tilde_list,
                              alpha_tilde){

n_components <- length(Z_tilde_list)
if(length(alpha_tilde) != n_components - 1){
  stop("length(alpha_tilde) must equal length(Z_tilde_list) - 1")
}

Z_tilde <- Z_tilde_list[[1]]

for(a in 2:n_components){
  Z_tilde <- Z_tilde + Z_tilde_list[[a]]*exp(alpha_tilde[a - 1])
}
return(Z_tilde)
}
