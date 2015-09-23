# truncate1 a real number The function \code{truncate1} is used to avoid
# instability in \code{R} due to very smal or very big numbers.  It truncate1s
# numbers so that they lie between \code{exp(-700)} and \code{exp(700)}.
# @param x A number to be truncate1d @return The truncated value of \code{x}
# @author Pierre Gaillard <pierre@@gaillard.me> @keywords ~kwd1 ~kwd2

truncate1 <- function(x) {
  pmin(pmax(x, exp(-700)), exp(700))
} 
