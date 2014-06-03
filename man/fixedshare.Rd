\name{fixedshare}
\alias{fixedshare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fixed-share aggregation rule
}
\description{
  The function \code{fixedshare} performs the Fixed-share aggregation rule. 
  It considers a sequence \code{y} of observations to be predicted sequentially with the help of experts advices \code{x}. 
  The forms at each instance \code{t} a prediction by assigning weight to the experts advices and combining them.
  The main difference with \code{\link{ewa}} aggregation rule rely in the mixing rate \code{alpha}\eqn{\in [0,1]} wich considers at each instance a small probability \code{alpha} to have a rupture in the sequence and that the best expert may change. Fixed-share aggregation rule can thus compete with the best sequence of experts that can change a few times (see \code{\link{bestShifts}}), while \code{\link{ewa}} can only compete with the best fixed expert.
}
\usage{
  fixedshare(y, experts, eta, alpha, 
    awake = NULL, 
    loss.type = "squareloss", 
    loss.gradient = TRUE, 
    w0 = NULL)
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
    we can put \code{awake[t,k]=0} so that the mixture does not    consider expert \code{k} in	the mixture to predict \code{Y_t}.
  }
  \item{eta}{
    %%     ~~Describe \code{eta} here~~
    A positive learning rate. The bigger it is the faster the aggregation rule will learn from
   observations and experts performances. However too hight values lead to unstable
   weight vectors and thus unstable predictions. See \code{\link{ewaCalib}} to 
   automatically and sequentially tun this parameter.
}
  \item{alpha}{
%%     ~~Describe \code{eta} here~~
  A mixing rate in \code{[0,1]}. A parameter that evaluates the probability of rupture
   at each instance. Usually it should be close to 0 to perform well.
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
  \item{w0}{
   %%     ~~Describe \code{w0} here~~
   prior weights vector of the experts. 
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  \item{weights }{a matrix of dimension \code{c(T,N)}, with \code{T} the number
   of instances to be predicted and \code{N} the number of experts. Each row 
   contains the convex combination to form the predictions}  \item{prediction }{
    A vector of length \code{T} that contains the predictions
    outputted by the aggregation rule.
  }
   \item{cumulatedloss }{the cumulated loss suffered by the aggregation rules}
   \item{lastweights }{the last weights returned by the aggregation rule if all experts were
   active. It is usefull if we want to relaunch the aggregation rule with this weight vector
   as prior.}
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{fixedshareHour}}, \code{\link{fixedshareCalib}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
