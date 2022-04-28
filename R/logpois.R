##################### Poisson-Type (Log) Density #####################
logpois <- function(x,intensity){

  if(!is.finite(intensity)){
    intensity <- 0
  }
  if((x==0)&(intensity==0)){
    return(0)
  } else{
    if(intensity <= 0){
      return(-Inf)
    } else{
      return(x*log(intensity) - lgamma(x + 1) - intensity)
    }
  }
}
