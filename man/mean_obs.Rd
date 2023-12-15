\name{mean_obs}
\alias{mean_obs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function calculates the harmonic mean of two observations.
}
\description{
This function calculates the harmonic mean of two observations with corresponding variances.}
\usage{
mean_obs(y1, y2, tau1sqrd, tau2sqrd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1}{
The value of the first observation.
}
  \item{y2}{
The value of the second observation.
}
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
tau_eq_sqrd
}

\examples{
y1 <- 3
y2 <- 3.1

tau1sqrd <- 1e2
tau2sqrd <- 1e1

mean_obs(y1, y2, tau1sqrd, tau2sqrd) 

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.