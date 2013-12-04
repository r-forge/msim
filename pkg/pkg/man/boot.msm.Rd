\name{boot.msm}
\alias{boot.msm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap Correction for Method of Simulated Moments
}
\description{
boot.msm utilizes the estimates produced from msm in order to hopefully eliminate some of the bias in the MSM estimates. Currently, the running time for the function is fairly long, and the default is 10 which seems small, but seems to produce better estimates in simulations.}
\usage{
boot.msm(msm.est = NULL, boot.message = TRUE, boot.sim = 10, family = "binomial", 
nsim = 1, num_MC_sims = 10000, num_subs = NULL, obs_per_sub = NULL, start = c(0, 1), method = "nleqslv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{msm.est}{
Estimates produced by the msim \{msm\} function; default is NULL}
  \item{boot.message}{
Logical; allows user to display message regarding standard error calculations
}
  \item{boot.sim}{
Number of bootstrap samples required; default is 10
}
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
Vector of length num_subs; Index of j (number of observations per subject); default is NULL
}
\item{start}{
Vector of starting values for (mu,sigma); default values are (0,1)
}
  \item{method}{
One of ("multiroot","optim","nleqslv"); default is multiroot. This determines the solver utilized within the MSM. If multiroot is selected, the function will use the multiroot{rootSolve} function. If optim is selected, the function will use the optim{base} function to minimize the Euclidean norm of the system. If nleqslv is chosen, nleqslv{nleqslv} will solve the system of equations using the Newton method.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{boot.mu}{Averages boot.sim estimates of mu}
 \item{boot.mu.se}{Averages boot.sim estimates of root mean squared error mu}
 \item{boot.sigma}{Averages boot.sim estimates of sigma}
 \item{boot.sigma.se}{Averages boot.sim estimates of root mean squared error sigma}
 \item{boot.sigma2}{Averages boot.sim estimates of sigma2 based on squaring estimates of sigma}
 \item{boot.sigma2.se}{Averages boot.sim estimates of root mean squared error sigma2}
}
\references{
Jiang, J. (1998). Consistent Estimators in Generalized Linear Mixed Models. Journal of the American Statistical Association, 93, 720–729.

Jiang, J. and Zhang, W. (2001). Robust estimation in generalized linear mixed models. Biometrika, 88, 753–765.
}
\author{
Lindsey Dietz}
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
