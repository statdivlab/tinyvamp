
poisson_criterion <- function(W,
                              means,
                              wts = NULL){
  n <- nrow(W)
  J <- ncol(W)
  wt_list <- lapply(1:n,
                    function(i) wts[1:J + (i - 1)*J])
  wts <- do.call(rbind,wt_list)

  rm(wt_list)

  if( (ncol(W) != ncol(means)) | (nrow(W) != nrow(means))){
    stop("W and means must have the same dimensions.")
  }

  if(is.null(wts)){
    return(sum(sapply(1:nrow(means), function(i) sapply(1:ncol(means),
                                                        function(j)
                                                          -logpois(W[i,j],means[i,j])))))
  }
  else{
    if((ncol(W) != ncol(wts)) | (nrow(W) != nrow(wts))){
      stop("Weight matrix wts must have same dimesions as W.")
    }
    return(sum(sapply(1:nrow(means), function(i) sapply(1:ncol(means),
                                                        function(j)
                                                          -wts[i,j]*logpois(W[i,j],means[i,j])))))

  }

}
