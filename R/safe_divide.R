##################### Division for Log Likelihood Derivatives #####################
safe_divide <- function(numer,denom, zero_divisor_penalty = 0){
  if((numer ==0 )&(denom == 0)){
    return(1)
  } else{
    if(denom == 0){
      return(zero_divisor_penalty)
    }
  }
  return(numer/denom)
}
