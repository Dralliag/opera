\name{mae}
\alias{mae}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Mean Absolute Error
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function \code{mae} computes the average absolute error of a sequence.
}
\usage{
  mae(x,y, awake = NULL, na.rm=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    A vector containing the sequence of predictions to be evaluated.
  }
  \item{y}{
    %%     ~~Describe \code{Y} here~~
    A vector containing the sequence of observations to be predicted.
  }
  \item{awake}{
    %%     ~~Describe \code{awake} here~~
    A vector specifying the activation coefficients of the predictions in \code{x}. Its entries lie in \code{[0,1]}. 
  }
  \item{na.rm}{
    logical. Should missing values (including NaN) be removed?
  }
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  Average absolute error of the sequence
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{rmse}}, \code{\link{mape}}, \code{\link{loss}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
