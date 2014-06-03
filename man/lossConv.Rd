\name{lossConv}
\alias{lossConv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Error of a convex combination
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  The function \code{lossConv} computes the loss suffered by the constant convex mixture \code{p}.
}
\usage{
  lossConv(p, y, experts, 
    awake, 
    loss.type = "squareloss")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
    %%     ~~Describe \code{p} here~~
    A non-negative vector that sum to one. It is of length \code{N}, the number of experts.
  }
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
  \item{awake}{
    %%     ~~Describe \code{awake} here~~
    A matrix specifying the activation coefficients of the experts. Its entries lie in \code{[0,1]}. 
    Needed if some experts are specialists and do not always form and suggest prediction. 
    If the expert number \code{k} at instance \code{t} does not form any prediction of observation \code{Y_t}, 
    we can put \code{awake[t,k]=0} so that the mixture does not consider expert \code{k} in	the mixture to predict \code{Y_t}.
  }
  \item{loss.type}{
    %%     ~~Describe \code{loss.type} here~~
    A string specifying the loss function considered to evaluate the performance. It can be "squareloss",
    "mae", "mape", or "pinballloss". See \code{\link{loss}} for more details.
  }
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  The average errors suffered by the mixture. Note that if \code{loss.type = "squareloss"}, the \code{rmse} is returned.
  Note also that the instance are weighted according to the number of activated experts. 
}

\note{
%%  ~~further notes~~
The function \code{lossConv} is for instance used to compute the best convex combination in hindsight.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as 
  \code{\link{loss}}, \code{\link{bestConvex}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
