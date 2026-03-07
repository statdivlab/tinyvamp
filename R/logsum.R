#' Sum positive numbers via their logarithms
#' Given x and y, this function returns log(exp(x) + exp(y)).
#'
#' @author David Clausen
#'
#' @param x the thing to be summed.
#' @param y the thing to be added to x. If absent, just the sum over vector-valued x is taken.
#'
sum_of_logs <- function(x, y = NULL) {
  #calculates log(sum(exp(x)))

  if (is.null(y)) {
    if (length(x) == 1) {
      return(x)
    } else {
      if (sum(!is.infinite(x)) == 0) {
        if (sum(sign(x) != -1) == 0) {
          return(-Inf)
        } else {
          return(Inf)
        }
      }
      x <- x[order(x, decreasing = TRUE)]
      return(x[1] + log(1 + sum(exp(x[2:length(x)] - x[1]))))
    }
  } else {
    if (length(x) != length(y)) {
      stop("If y is provided, x and y must be the same length")
    }
    n <- length(x)
    maxes <- pmax(x, y)
    mins <- pmin(x, y)
    return(maxes + log(1 + exp(mins - maxes)))
  }
}


diff_of_logs <- function(larger, smaller) {
  return(larger + log(1 - exp(smaller - larger)))
}


signed_log_sum <- function(log_magnitudes, signs) {
  # if no nonzero summands, return zero
  if (sum(signs != 0) == 0) {
    return(c(0, 0))
  }

  # if no negative summands, return sum of positive summands
  if (sum(signs == -1) == 0) {
    return(c(sum_of_logs(log_magnitudes[signs == 1]), 1))
  }

  # if no positive summands, return sum of negative summands (with sign -1)
  if (sum(signs == 1) == 0) {
    return(c(sum_of_logs(log_magnitudes[signs == -1]), -1))
  }

  # otherwise sum positive and negative magnitudes separately
  positives <- log_magnitudes[signs == 1]
  negatives <- log_magnitudes[signs == -1]

  pos_logsum <- sum_of_logs(positives)
  neg_logsum <- sum_of_logs(negatives)

  # if positive magnitude equals negative magnitude, return zero
  if (pos_logsum == neg_logsum) {
    return(c(0, 0))
  }

  # otherwise return difference of positive and negative terms with appropriate sign
  return(c(
    diff_of_logs(max(pos_logsum, neg_logsum), min(pos_logsum, neg_logsum)),
    ifelse(pos_logsum > neg_logsum, 1, -1)
  ))
}
