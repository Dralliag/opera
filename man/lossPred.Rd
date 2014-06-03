\name{lossPred}
\alias{lossPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Instantaneous loss suffered by a prediction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
    The function \code{lossPred} computes the loss of a prediction \code{x} of an obstervation \code{y}
    knowing that we want to correct the error commited by the prediction \code{pred}. 
    It is used in the mixing aggregation rule (see \code{\link{ewa}} for instance) to compute the loss of the experts at each instance.
    We can choose \code{pred} to be the prediction of \code{y} outputed by the aggregation rule.

}
\usage{
  lossPred(x, y, pred = NULL, 
    loss.type = "squareloss", 
    loss.gradient = FALSE, 
    tau = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    A vector of \code{N} prediction of the observation \code{y} to be evaluated.
  }
  \item{y}{
     %%     ~~Describe \code{Y} here~~
    A number containing the observation.    
  }
  \item{pred}{
%%     ~~Describe \code{pred} here~~
    A reference prediction that the predictions in \code{x} aim to correct.
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
    A vector containing the loss suffered by the \code{N} predictions in \code{x}.
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{loss}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
