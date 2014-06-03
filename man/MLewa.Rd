\name{MLewa}
\alias{MLewa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Exponentially Weighted Average aggregation rule with multiple learning rates
}
\description{
  The function \code{MLewa} performs the Exponentially Weighted Average aggregation rule with different learning rates for each expert.
  The learning rates are calibrated using theoretical values.
}
\usage{
  MLewa(y, experts, awake = NULL, 
    loss.type = "squareloss", 
    loss.gradient = TRUE, 
    period = 1, 
    href = 1)
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
  \item{awake}{
    %%     ~~Describe \code{awake} here~~
  	A matrix specifying the activation coefficients of the experts. Its entries lie in \code{[0,1]}. 
  	Needed if some experts are specialists and do not always form and suggest prediction. 
  	If the expert number \code{k} at instance \code{t} does not form any prediction of observation \code{Y_t}, 
    we can put \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in the mixture to predict \code{Y_t}.
  }
  \item{loss.type}{
    %%     ~~Describe \code{loss.type} here~~
    A string specifying the loss function considered to evaluate the performance. It can be "squareloss",
    "mae", "mape", or "pinballloss". See \code{\link{loss}} for more details.
  }
  \item{loss.gradient}{
    %%     ~~Describe \code{loss.gradient} here~~
    A boolean. If TRUE (default) the aggregation rule will not be directly applied to the loss function
    at hand but to a gradient version of it. The aggregation rule is then similar to gradient 
    descent aggregation rule.  
  }
  \item{period}{
    %%     ~~Describe \code{period} here~~
    The number of instants in each day.
  }
  \item{href}{
    %%     ~~Describe \code{href} here~~
    A number in \code{[1,period]} specifying the instant in the day when the aggregation rule can update its weights.
    It should lie in the interval \code{c(1,period)}.
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
    A vector of length \code{T} that contains the predictions
    outputted by the aggregation rule.
  }
  \item{eta}{
    A vector of length \code{T} containing the sequence of learning rates chosen by the aggregation rule.
  }
}
\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ewa}}, \code{\link{MLpol}}, \code{\link{MLprod}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
