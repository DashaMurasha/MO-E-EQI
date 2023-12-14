\name{mult_EQI}
\alias{mult_EQI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multi-Objective EQI
%%  ~~function to do ... ~~
}
\description{
Calculates expected quantile improvement for two objectives
}
\usage{
Multi-Objective EQI(newdata,design_X, model_f1, model_f2, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{newdata}{
n_new x k data frame of new design points to test for EQI, where n_new is the number of new sample points
}
  \item{design_X}{
n x k data frame of sample locations for each objective (same locations for both objectives at the moment)
}
  \item{model_f1}{
Emulator for the first objective
}
  \item{model_f2}{
Emulator for the second objective
}
  \item{beta}{
Quantile level
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Outputs:
  \item{metric}{ - -EQI or -log10(EQI), depending on the Option (So far only implemented Option=='NegLogExpImp')
  }
  \item{Pq1,Pq2}{ - non-dominated objective function values on the Pareto front}
  \item{PX}{ - locations of non-dominated solutions}
  \item{s_Q1, s_Q2}{ - standard deviations at the proposed (newdata) points}
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
