\name{tau_new_func}
\alias{tau_new_func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Future noise estimate
}
\description{
This is a function of a computational budget, evaluating how noise reduces with increasing computational budget.
}
\usage{
tau_new_func(n_rep, sd_noise)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_rep}{
Numer of repetitions at each new location.
}
  \item{sd_noise}{
A vector of length 2 providing the default noise in both objectives of the process if only 1 run of the model is performed.
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
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line