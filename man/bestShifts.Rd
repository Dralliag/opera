\name{bestShifts}
\alias{bestShifts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
best sequence of experts oracle
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function \code{bestShifts} computes for all number m of stwitches the sequence of experts
with at most $m$ shifts that would have performed the best to predict the sequence of
observations in \code{y}.
}
\usage{
bestShifts(y, experts, awake, 
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
 %%  Returns a matrix of dimension \code{c(T,N,3)} where \code{T} is the number of instance to be predicted (i.e., the length of the sequence \code{y}) and \code{N} is the number of experts.
  A matrix \code{L} of dimension \code{c(T,N,3)} where the third dimension is the type of loss (1:squareloss, 2:mae, 3:mape), and the value of $L(m,k,l)$ is the loss (determined by \code{l}) suffered by the best sequence of expert with at most $m-1$ shifts and finishing with expert number $k$.     
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{bestConvex}}, \code{\link{bestShiftsDay}}, \code{\link{loss}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
