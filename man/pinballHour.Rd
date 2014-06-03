\name{pinballHour}
\alias{pinballHour}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Follow the penalized leader aggregation rule for pinball loss
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  The function \code{pinballHour} is a mixing aggregation rule for quantile regression.
  At each instance, it forms the mixture by performing a convex minimisation. It chooses the mixture that minimizes among the past a penalized criterion based on cumulated pinball loss.
}
\usage{
  pinballHour(y, experts, lambda, 
    href = 1, period = 1, 
    w0 = NULL, tau = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    %%     ~~Describe \code{Y} here~~
    A vector containing the observations to be predicted.
  }
  \item{experts}{
    %%     ~~Describe \code{experts} here~~
    A matrix containing the experts forecasts. Each column corresponds to the 
    predictions proposed by an expert to predict \code{Y}. It has as many columns
    as there are experts.
  }
  \item{lambda}{
    %%     ~~Describe \code{lambda} here~~
    A positive penalty coefficient.
  }
  \item{href}{
    %%     ~~Describe \code{href} here~~
    A number in \code{[1,period]} specifying the instant in the day when the aggregation rule can update its weights.
    It should lie in the interval \code{c(1,period)}.
  }
  \item{period}{
    %%     ~~Describe \code{period} here~~
    The number of instants in each day.
  }
  \item{w0}{
    %%     ~~Describe \code{w0} here~~
    A vector containing the prior weights of the experts. 
  }
  \item{tau}{
    %%     ~~Describe \code{tau} here~~
    A number in \code{[0,1]} describing the quantile to be predicted. Used only if \code{loss.type = "pinballloss"}.
  }
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{weights}{
    A matrix of dimension \code{c(T,N)}, with \code{T} the number of instances to be predicted and \code{N} the number of experts. 
    Each row contains the convex combination to form the predictions.
  }  
  \item{prediction}{
    A vector of length \code{T} that contains the quantiles predictions
    outputted by the aggregation rule.
  }
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{loss}}, \code{\link{ridge}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
