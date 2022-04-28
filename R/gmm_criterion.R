
gmm_criterion <- function(W_long,
                          means_long,
                          inv_wts){


squerror_long <- (W_long - means_long)^2

gmm_crit <- 0.5*sum(sapply(1:length(means_long),
       function(i){
         ifelse(inv_wts[i]>0,
                squerror_long[i]/inv_wts[i],
                ifelse(squerror_long[i] ==0,0,Inf))
       }
         ))

return(gmm_crit)
}

