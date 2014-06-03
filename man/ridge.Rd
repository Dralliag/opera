\name{ridge}
\alias{ridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Ridge aggregation rule with fixed learning parameter
}
\description{
  The function \code{ridge} performs ridge aggregation rule by minimizing at each instance a penalized criterion.
  It forms at each instance linear combination of the experts' forecasts and can assign negative weights that not necessarily sum to one.
  It is useful if the experts are biased or correlated. It cannot be used with specialized experts.
}
\usage{
  ridge(y, experts, lambda, 
    w0 = NULL)
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
  \item{w0}{
    %%     ~~Describe \code{w0} here~~
    A vector containing the prior weights of the experts. 
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
  \code{\link{ridgeHour}},\code{\link{ridgeCalib}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
