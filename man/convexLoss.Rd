\name{convexLoss}
\alias{convexLoss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  agerage error suffered by convex combination
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function \code{convexLoss} computes the loss incured by a fixed convex combination of experts.
}
\usage{
  convexLoss(weights, y, experts, awake)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weights}{
%%     ~~Describe \code{weights} here~~
   vector of weights assigned to each experts. It must be the same length as the
   number of columns of \code{experts}. It entries must be non-negative and sum to one.
}
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
    we can put \code{awake[t,k]=0} so that the mixture does not    consider expert \code{k} in	the mixture to predict \code{Y_t}.
  }
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A dataframe containing several measures of losses suffered by the mixture
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{loss}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
