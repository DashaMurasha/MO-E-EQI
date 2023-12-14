\name{q_mean}
\alias{q_mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The quantile calculation
}
\description{
Calculates current quantile at level beta for a given emulator
}
\usage{
q_mean(x, model, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A point, a vector or a matrix of points in the parameters' space to evaluate quantiles at
}
  \item{model}{
Current emulator model of \code{\link[DiceKriging]{km}} type
}
%  \item{beta}{
%Quantile level between [0.5,1)
%}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Daria Semochkina
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
