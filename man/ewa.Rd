\name{ewa}
\alias{ewa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Exponentially Weighted Average aggregation rule
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  The function \code{ewa} performs the exponentially weighted average aggregation rule. 
  It considers a sequence \code{y} of observations to be predicted sequentially with the help of experts advices \code{x}.   
  The forms at each instance \code{t} a prediction by assigning weight to the experts advices and combining them.
}
\usage{
  ewa(y, experts, eta, awake = NULL, 
    loss.type = "squareloss", 
    loss.gradient = TRUE, 
    w0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    %%     ~~Describe \code{Y} here~~
    A vector that contains the observations to be predicted.
  }
  \item{experts}{ 
    %%     ~~Describe \code{X} here~~
    A matrix containing the experts forecasts. Each column corresponds to the 
    predictions proposed by an expert to predict \code{Y}. It has as many columns
    as there are experts.
  }
  \item{awake}{
    %%     ~~Describe \code{awake} here~~
    A matrix specifying the activation coefficients of the experts. Its entries lie in \code{[0,1]}. 
    Needed if some experts are specialists and do not always form and suggest prediction. 
    If the expert number \code{k} at instance \code{t} does not form any prediction of observation \code{Y_t}, 
    we can put \code{awake[t,k]=0} so that the mixture does not    consider expert \code{k} in  the mixture to predict \code{Y_t}.
  }
  \item{eta}{
    %%     ~~Describe \code{eta} here~~
    A positive learning rate. The bigger it is the faster the aggregation rule will learn from observations and experts performances. However too hight values lead to unstable
    weight vectors and thus unstable predictions. See \code{\link{ewaCalib}} to 
    automatically and sequentially tun this parameter.
}
  \item{loss.type}{
    %%     ~~Describe \code{loss.type} here~~
    A string specifying the loss function considered to evaluate the performance. It can be "squareloss", "mae", "mape", or "pinballloss". See \code{\link{loss}} for more details.
  }
  \item{loss.gradient}{
    %%     ~~Describe \code{loss.gradient} here~~
    A boolean. If TRUE (default) the aggregation rule will not be directly applied to the loss function at hand but to a gradient version of it. The aggregation rule is then similar to gradient descent aggregation rule.  
  }
  \item{w0}{
    %%     ~~Describe \code{w0} here~~
    A vector containing prior weights vector of the experts. 
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{weights }{
    A matrix of dimension \code{c(T,N)}, with \code{T} the number of instances to be predicted and \code{N} the number of experts. 
    Each row contains the convex combination to form the predictions.
  }  
  \item{prediction }{
    A vector of length \code{T} that contains the predictions
    outputted by the aggregation rule.
  }
  \item{cumulatedloss }{
    The cumulated loss suffered by the aggregation rule.
  }
  \item{regret }{
    An array that contains the cumulated regret suffered by the aggregation rule against each expert.
  }
%% ...
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ewaHour}}, \code{\link{ewaCalib}}
}
\examples{
library('opera')              # load the package
set.seed(1)                   

T = 100                       # number of instances
t = 1:T                       # instances
Y = cos(5*2*pi*t / T)         # sequence to be predicted

X1 = Y + 0.1*rnorm(T)         # first expert (with small average error)
X2 = Y + 0.3*rnorm(T)         # second expert
awake1 = rep(c(rep(1,9),0),T/10) # the first expert is not always available
awake2 = rep(1,T)             # the second expert is always available

X = cbind(X1,X2)              # matrix of experts
awake = cbind(awake1,awake2)  # activation matrix

matplot(X, type='l', col=2:3) # plot experts' predictions
lines(Y)                      # plot observations

# Performance of the experts
cat('Expert 1, rmse :', rmse(X1,Y,awake=awake1), '\n')
cat('Expert 2, rmse :', rmse(X2,Y,awake=awake2), '\n')

# Performance of taking expert 1 if available, expert 2 otherwise
X3 = X1 * awake[,1] + X2 * (1-awake[,1])
cat("Best sequence of experts in hindsight, rmse :", rmse(X3,Y), '\n\n')


# EWA with fixed learning rate
mixture = ewa(y=Y, experts=X, eta=1, awake=awake, loss.type='squareloss', loss.gradient=FALSE) 
# plot weights assigned to both experts (when an expert is not available its weight is 0)
matplot(mixture$weights, type='l', main='EWA with fixed learning rate', col=2:3) 
cat('EWA mixture, rmse :', rmse(mixture$prediction,Y), '\n')

# ewa algorithm with gradient loss function
mixture = ewa(y=Y, experts=X, eta=1, awake=awake, loss.type='squareloss', loss.gradient=TRUE) 
matplot(mixture$weights, type='l', main='EWA with gradient losses', col=2:3) 
cat('EWA mixture with gradient losses, rmse :', rmse(mixture$prediction,Y), '\n')

# ewa algorithm with automatic calibration of the learning parameter
mixture = ewaCalib(y=Y, experts=X, awake = awake, loss.type='squareloss', loss.gradient=TRUE) 
matplot(mixture$weights, type='l', main = 'Automatic EWA', col=2:3) 
cat('EWA mixture with automatic tuning, rmse :', rmse(mixture$prediction,Y), '\n')

# MLpol aggregation rule
mixture = MLpol(y=Y, experts=X, awake= NULL, loss.gradient=TRUE)
mixture$prediction = apply(mixture$weights*X, 1, sum)
matplot(mixture$weights, type='l', main = 'MLpol mixture', col=2:3, ylim = c(0,1))
cat('MLpol mixture, rmse :', rmse(mixture$prediction,Y), '\n')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
