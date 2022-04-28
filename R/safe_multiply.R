
safe_multiply <- function(x,y){
  if(length(x)!= length(y)){
    stop("Arguments x and y must have equal length")
  }

  product <- x*y

  #define inf*0 = 1; -inf*0 = -1
  product[is.infinite(x) & y == 0] <- sign(x)[is.infinite(x) & y == 0]
  product[is.infinite(y) & x == 0] <- sign(y)[is.infinite(y) & x == 0]

  return(product)

}
