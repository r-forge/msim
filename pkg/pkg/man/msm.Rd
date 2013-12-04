\name{msm}
\alias{msm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Method of Simulated Moments
}
\description{
Function to produce averaged estimates of multiple runs of method of simulated moments.
}
\usage{
msm(family = "binomial", nsim = 1, num_MC_sims = 10000, num_subs = NULL, obs_per_sub = NULL, y.i = NULL, 
start = c(0, 1), message = TRUE, true.mu = NULL, true.sigma = NULL, method = "nleqslv")
}

\arguments{
  \item{family}{
Exponential family to draw from; currently only accepts "binomial" or "poisson"
}
  \item{nsim}{
Number of simulations of MSM estimates; default is 1
}
  \item{num_MC_sims}{
Number of values used to produce one Monte Carlo estimate for MSM; default is 10000
}
  \item{num_subs}{
Index of i (number of subjects); default is NULL
}
  \item{obs_per_sub}{
Vector of length num_subs; Index of j (number of observations per subject); default is NULL}
  \item{y.i}{
Sums over j of the y_ij, produced by simulate.fun or provided by user; default is NULL
}
  \item{start}{
Vector of starting values for (mu,sigma); default values are (0,1)
}
  \item{message}{
Logical; allows user to display message regarding standard error calculations}
  \item{true.mu}{
True value of mu; default is NULL
}
  \item{true.sigma}{
True value of sigma; default is NULL
}
  \item{method}{
One of ("multiroot","optim","nleqslv"); default is multiroot. This determines the solver utilized within the MSM. If multiroot is selected, the function will use the multiroot{rootSolve} function. If optim is selected, the function will use the optim{base} function to minimize the Euclidean norm of the system. If nleqslv is chosen, nleqslv{nleqslv} will solve the system of equations using the Newton method.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{mu}{Averages nsim estimates of mu}
\item{mu.se}{Averages nsim estimates of root mean squared error mu}
\item{sigma}{Averages nsim estimates of sigma}
\item{sigma.se}{Averages nsim estimates of root mean squared error sigma}
\item{sigma2}{Averages nsim estimates of sigma2 based on squaring the estimate of sigma}	
\item{sigma2.se}{Averages nsim estimates of root mean squared error sigma2}	
}
\references{
Jiang, J. (1998). Consistent Estimators in Generalized Linear Mixed Models. Journal of the American Statistical Association, 93, 720–729.

Jiang, J. and Zhang, W. (2001). Robust estimation in generalized linear mixed models. Biometrika, 88, 753–765.}
\author{
Lindsey Dietz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{optim}}
\code{\link{multiroot}}
\code{\link{nleqslv}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
