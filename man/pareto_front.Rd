\name{pareto_front}
\alias{pareto_front}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function finds the Pareto front for a 2-D objective space.
}
\description{
This function identifies the non-dominated solutions for the given set of solutions for two objectives.
}
\usage{
pareto_front(y1_obs, y2_obs, design)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1_obs}{
The observations for the first objective.
}
  \item{y2_obs}{
The observations for the second objective.
}
  \item{design}{
The inputs combinations corresponding to each of the observations.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Noise distribution parameter
a=0.5

# Test functions
f1 <- function(x){
  x1 <- unlist(x[,1])
  x2 <- unlist(x[,2])
  1-sin(x1)+x2/10+a*cos(0.3)
}
f2 <- function(x){
  x1 <- unlist(x[,1])
  x2 <- unlist(x[,2])
  1-cos(x1)+x2/3+a*sin(0.3)

}


design <- matrix(runif(100,0,1), 50, 2)
                      
y1_obs <- f1(design)
y2_obs <- f2(design)

pareto_front(y1_obs,y2_obs,design)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.