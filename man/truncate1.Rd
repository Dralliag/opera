\name{truncate1}
\alias{truncate1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  truncate1 a real number 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  The function \code{truncate1} is used to avoid instability in \code{R} due to very smal or very big numbers.
  It truncate1s numbers so that they lie between \code{exp(-700)} and \code{exp(700)}.
}
\usage{
  truncate1(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    %%     ~~Describe \code{x} here~~
    A number to be truncate1d
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  The truncate1d value of \code{x}
}

\author{
  Pierre Gaillard <pierre-p.gaillard@edf.fr>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
