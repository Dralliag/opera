\name{ridgeCalib}
\alias{ridgeCalib}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ridge aggregation rule with automatic calibration of smoothing parameters
}
\description{
  The function \code{ridgeCalib} performs ridge aggregation rule with
  automatic calibration of the smoothing parameter \code{lambda} by performing an optimization
  on a finite grid. See \code{\link{ridge}} and \code{\link{ridgeHour}} for more details.
}
\usage{
  ridgeCalib(y, experts, 
    gridlambda = 1, 
    href = 1, period = 1, 
    w0 = NULL, 
    trace = F)
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
  \item{gridlambda}{
    A vector containing the initial possible values of the smoothing parameter \code{lambda}.
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
  \item{trace}{
%%     ~~Describe \code{trace} here~~
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
  \item{lambda}{
    A vector of length \code{T} containing the sequence of penalty coefficient chosen by the aggregation rule.
  }
  \item{grid}{
    A vector containing the final grid of potential penalty coefficient considered by the aggregation rule.
  }
  \item{loss}{
    The error suffered by the aggregation rule determined by \code{loss.type}.
    If \code{loss.type = 'squareloss'}, the \link{rmse} is computed.
  }
  \item{gridloss}{
    A vector of the same length as grid containing the errors suffered by the \code{\link{ridgeHour}} aggregation rule 
    if it had picked the fixed learning rates in \code{grid}.
  }
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ridge}},\code{\link{ridgeHour}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
