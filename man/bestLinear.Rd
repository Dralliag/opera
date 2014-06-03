\name{bestLinear}
\alias{bestLinear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
best linear oracle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function \code{bestLinear} computes the best linear combination oracle or in other words the fixed linear combination of experts in \code{experts} that would have 
perfomed the best to predict the observations in \code{y}.
}
\usage{
  bestLinear(y, experts, lambda = 1, 
    awake = NULL, 
    loss.type = "squareloss")
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
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
	Smoothing parameter. No smoothing corresponds to \code{lambda = 0}, however
	the linearn system might be singular.  
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
}
\value{
  \item{loss}{ 
    The loss suffered by the best fixed linear combination.
  }
  \item{u}{
    A vector containing the best fixed linear combination chosen in hindsight.
  }
  \item{prediction}{
    A vecor containing the predictions of the best fixed linear combination.
  }
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{bestConvex}}, \code{\link{bestShifts}}, \code{\link{loss}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
