\name{design_repetitions}
\alias{design_repetitions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matches rows in matrices or data.frames
}
\description{
This function matches the rows of a matrix or data.frame to rows of another matrix or data.frame. Thus the function returns a matrix with the row numbers of (first) matches of all its first argument in its second.
}
\usage{
design_repetitions(x, design)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix or data.frame
}
  \item{design}{
Matrix or data.frame
}
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

x <- expand.grid(x.1 = sample(1:2, 100,replace = T),
                       x.2 = sample(3:4, 100,replace = T))
                       
design <- matrix(c(1,3,
                   2,3,
                   1,4,
                   2,4),4,2, byrow=T)

design_repetitions(x, design)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
