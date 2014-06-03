\name{mape}
\alias{mape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Mean Absolute Percentage Error
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  The function \code{mape} computes the mean absolute percentage error of a sequence of predictions. It is an equivalent way of running \code{mean(loss(x,y,loss.type='mape'))} (see \code{\link{loss}}).
}
\usage{
  mape(x, y, na.rm=TRUE)
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
  The mean absolute percentage error of the sequence of prediction \code{x}.
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
   \code{\link{mae}}, \code{\link{rmse}}, \code{\link{loss}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
