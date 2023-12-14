\name{tau_eq_sqrd}
\alias{tau_eq_sqrd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function calculates the joint variance of a repeated observation.
}
\description{
This function calculates the variance of a repeated observation, where the first observation's variance is the first input and the variance of the added observation is the second input. Both observations are at the same input.
}
\usage{
tau_eq_sqrd(tau1sqrd, tau2sqrd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau1sqrd}{
The variance of the first observation.
}
  \item{tau2sqrd}{
The variance of the second observation.
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
mean_obs
}
\examples{

tau1sqrd <- 1e2
tau2sqrd <- 1e1

tau_eq_sqrd(tau1sqrd, tau2sqrd) 

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
