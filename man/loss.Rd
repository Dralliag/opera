\name{loss}
\alias{loss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Error suffered by a sequence of prediction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  The function \code{loss} computes  the sequence of instantaneous losses suffered by the predictions in \code{x} to predict the observation in \code{y}. 
}
\usage{
  loss(x, y, loss.type = "squareloss", 
    tau = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    A vector of length \code{T} containing the sequence of prediction to be evaluated.
  }
  \item{y}{
    %%     ~~Describe \code{Y} here~~
    A vector of length \code{T} that contains the observations to be predicted.
  }
  \item{loss.type}{
    %%     ~~Describe \code{loss.type} here~~
    A string specifying the loss function considered to evaluate the performance. It can be "squareloss",
    "mae", "mape", or "pinballloss". See \code{\link{loss}} for more details.
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
  A vector of length \code{T} containing the sequence of instantaneous losses suffered by the prediction \code{x}.
}

\author{
  %%  ~~who you are~~
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{mape}}, \code{\link{rmse}}, \code{\link{lossConv}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
