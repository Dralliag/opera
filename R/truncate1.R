# truncate1 a real number The function \code{truncate1} is used to avoid
# instability in \code{R} due to very smal or very big numbers.  It truncate1s
# numbers so that they lie between \code{exp(-700)} and \code{exp(700)}.
# @param x A number to be truncate1d @return The truncated value of \code{x}
# @author Pierre Gaillard <pierre@@gaillard.me> @keywords ~kwd1 ~kwd2

truncate1 <- function(x) {
  is_sup <- max(x) > exp(700)
  is_inf <- min(x) <= exp(-700)
  
  res <- x
  
  if (is_sup) {
    res <- pmin(res, exp(700))
  }
  if (is_inf) {
    res <- pmax(res, exp(-700))
  }
  
  return(res)
} 