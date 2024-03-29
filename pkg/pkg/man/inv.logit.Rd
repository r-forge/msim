\name{inv.logit}
\alias{inv.logit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inverse Logit Function
}
\description{
Calculate the inverse logit function: f(x)=exp(x)/(1+exp(x))
}
\usage{
inv.logit(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a real valued scalar or vector
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
Lindsey Dietz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x<-c(1,2,3,4)
inv.logit(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
