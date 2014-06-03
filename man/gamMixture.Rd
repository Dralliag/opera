\name{gamMixture}
\alias{gamMixture}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Mixture aggregation rule using generalized additive model (beta version)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function \code{gamMixture} fits a general additive model (GAM) on the data to form weight vectors for the experts that can depend on exogeneous data.
The process is however not currently stable. We advice not using it yet. Use only at most one exogeneous variable.
}
\usage{
  gamMixture(y, experts, z, lambda, 
    nknots = 5, degree = 3, 
    loss.type = "squareloss", 
    href = 1, period = 1, 
    uniform = F, 
    knots = NULL, 
    tau = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    %%     ~~Describe \code{Y} here~~
    A vector that contains the observations to be predicted.
  }
  \item{experts}{
    %%     ~~Describe \code{experts} here~~
    A matrix containing the experts forecasts. Each column corresponds to the 
    predictions proposed by an expert to predict \code{Y}. It has as many columns
    as there are experts.
  }
  \item{z}{
    %%     ~~Describe \code{z} here~~
    A matrix containing the exogenous variables. Each column corresponds to covariable. It can currently support only one variable.
  }
  \item{lambda}{
    %%     ~~Describe \code{lambda} here~~
    A positive real giving the learning parameter. 
  }
  \item{nknots}{
    %%     ~~Describe \code{nknots} here~~
    The number of knots the splines have.
  }
  \item{degree}{
    %%     ~~Describe \code{degree} here~~
    The degree of the splines
  }
  \item{loss.type}{
    %%     ~~Describe \code{loss.type} here~~
    Loss function considered to evaluate the performance. It can be "squareloss",
    "mae", "mape", or "pinballloss". See \code{\link{loss}} for more details.
  }
  \item{href}{
    %%     ~~Describe \code{href} here~~
    A number in \code{[1,period]} specifying the instant in the daywhen the aggregation rule can update its weights.
    It should lie in the interval \code{c(1,period)}.
  }
  \item{period}{
    %%     ~~Describe \code{period} here~~
    The number of instants in each day.
  }
  \item{uniform}{
    %%     ~~Describe \code{uniform} here~~
    if we want the knots to be set uniformly (TRUE) or according to the quantiles (FALSE)
  }
  \item{knots}{
    %%     ~~Describe \code{knots} here~~
    An optional array containing the knots values.
  }
  \item{tau}{
    %%     ~~Describe \code{tau} here~~
    the quantile value we want to predict.
  }
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{weights}{
    A matrix containing weights.
  }
  \item{design}{
    The design matrix.
  }
  \item{knots}{
    An array containing the values of the knots.
  }
  \item{prediction}{
    A vector containing predictions of the mixture.
  }
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ridge}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
