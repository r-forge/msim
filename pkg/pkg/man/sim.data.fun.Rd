\name{sim.data.fun}
\alias{sim.data.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Data from logistic mixed model or poisson-normal model
}
\description{
This function will simulate the data for the logistic mixed model: logit(y_ij)=mu+sigma* theta_i or the poisson-normal model: log(y_ij)=mu+sigma* theta_i where theta_i is N(0,1) in either case.
}
\usage{
sim.data.fun(num_subs = NULL, obs_per_sub = NULL, true.mu = NULL, true.sigma = NULL,family="binomial", set.seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{num_subs}{
Index of i (number of subjects); default is NULL}
  \item{obs_per_sub}{
Vector of length num_subs; Index of j (number of observations per subject); default is NULL
}
  \item{true.mu}{
True value of mu; default is NULL
}
  \item{true.sigma}{
True value of sigma; default is NULL
}
  \item{family}{
Exponential family from which to generate top level of hierarchy; default is "binomial"
}
  \item{set.seed}{
Random seed start value for reproducibility; default is NULL
}
}
\details{
The user must enter values for num_subs, obs_per_sub, true.mu and true.sigma or the function will produce an error. The set.seed argument is optional.
}
\value{
\item{y}{a matrix of generated y_ij of dimension num_subs* max(obs_per_sub)}
\item{y.i}{a vector of the sums over j of y_ij; this will have length num_subs}
\item{n.i}{a vector of the number of observations for each subject; this should exactly match obs_per_sub}

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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
