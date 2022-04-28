get_gmm_inv_weights <- function(W_long,
                          means_long){

  squerror_long <- (W_long - means_long)^2

  pre_wts <- isoreg(means_long,squerror_long)

  means_long_ordered <- means_long[order(means_long)]
  squerror_long_ordered <- squerror_long[order(means_long)]

  # pre_wts$yf <- pre_wts$yf*(mean(1/pre_wts$yf))
  # sum(1/pre_wts$yf)

  # if any variances estimated to be zero at positive means
  # if(min(pre_wts$yf[means_long_ordered>0]) ==0){
  #   # if there are zero means, do linear interpolation
  #   if(sum(means_long_ordered == 0)<0){
  #     pre_wts$yf[means_long_ordered == 0] <- 0
  #
  #     max_index <- max(which(pre_wts$yf[means_long_ordered>0] ==0))
  #     (pre_wts$yf[means_long_ordered>0])[1:max_index] <-
  #       ((means_long_ordered[means_long_ordered>0])[1:max_index])*
  #       (pre_wts$yf[means_long_ordered>0])[max_index + 1]
  #   }
  #   #if no zero means, set min var equal to min nonzero var
  #   if(sum(means_long_ordered == 0) == 0){
  #     pre_wts$yf[pre_wts$yf ==0] <- min(pre_wts$yf[pre_wts$yf > 0])
  #   }
  #
  # }
  inv_wts <- numeric(length(W_long))
  inv_wts[order(means_long)] <- pre_wts$yf

  return(inv_wts)
}
