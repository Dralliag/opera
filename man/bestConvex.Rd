\name{bestConvex}
\alias{bestConvex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  best convex oracle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
	The function \code{bestConvex} computes the best convex combination oracle, or in other words the fixed convex combination of experts in \code{experts} that would have perfomed the best to predict the observations in \code{y}.
   
}
\usage{
  bestConvex(y, experts, awake = NULL, 
    loss.type = "squareloss", 
    niter = 1, method = "L-BFGS-B", 
    control = list(maxit = 100, trace = T))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    %%     ~~Describe \code{y} here~~
    vector that contains the observations to be predicted.
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
  \item{loss.type}{
    %%     ~~Describe \code{loss.type} here~~
    A string specifying the loss function considered to evaluate the performance. It can be "squareloss",
    "mae", "mape", or "pinballloss". See \code{\link{loss}} for more details.
  }
  \item{niter}{
%%     ~~Describe \code{niter} here~~
	  Number of times the convex optimisation process is repeated with different 
	  initial values. Note that if the experts are always activated, it does not need
	  to be greater than 1.
  }
  \item{method}{
    The optimization method to be used. See \code{\link{optim}}.
  }
  \item{control}{
%%     ~~Describe \code{control} here~~
	  Additional parameters that are passed to \code{optim} function.
}
}
\value{
  \item{loss}{ 
    The loss suffered by the best fixed convex combination.
  }
  \item{weights}{
    A vector containing the best fixed convex combination chosen in hindsight.
  }
  \item{prediction}{
    A vecor containing the predictions of the best fixed convex combination.
  }
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{bestLinear}}, \code{\link{bestShifts}}, \code{\link{loss}}, \code{\link{lossConv}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
