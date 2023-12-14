\name{pred_Q}
\alias{pred_Q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicts quantile values at the next time step
}
\description{
Calculates the predicted quantile values at new potential design points x and new noise level tau_new
}
\usage{
pred_Q(newdata, model, beta, tau_new)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{newdata}{
A point, a vector or a matrix of new points in the parameters' space to evaluate quantiles at.
}
  \item{model}{
Current emulator model
}
  \item{beta}{
Quantile level between [0.5,1)
}
  \item{tau_new}{
The new level of noise at location(s) x at the next time step
}
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
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (x)
%{
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
