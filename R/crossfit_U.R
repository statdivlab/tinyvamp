
crossfit_U <- function(W,
                       full_model,
                       null_model,
                       parallelize = FALSE){

n <- nrow(W)

if(!parallelize){
log_us <- numeric(n)

for(i in 1:n){
  # print(i)
  train_ind <- rep(TRUE,n)
  train_ind[i] <- FALSE
  log_us[i] <- calculate_U(W, full_model, null_model,
                           training_indicator = train_ind)
}
} else{
  log_us <- parallel::mclapply(1:n,
                     function(i){
                       train_ind <- rep(TRUE,n)
                       train_ind[i] <- FALSE
                       return(calculate_U(W,
                                          full_model,
                                          null_model,
                                          training_indicator = train_ind))

                     })
  log_us <- unlist(log_us)
}

return(exp(logsum::sum_of_logs(log_us) - log(n)))
}
