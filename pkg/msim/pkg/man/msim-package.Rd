\name{msim-package}
\alias{msim-package}
\alias{msim}
\docType{package}
\title{
Computes the Method of Simulated Moments

}
\description{
MSM is methodology which is computationally feasible and consistent. In applying the usual method of moments, one first identifies a set of sufficient statistics. A set of estimating equations is obtained by equating sample moments of the sufficient statistics to their expectations. Such expectations typically involve integrals, the highest dimension of which equals the number of sources of random effects. Expectations are then simulated. Finally, parameters are estimated by utilizing a Newton Raphson algorithm to solve the nonlinear system of equations.

In practice, these estimates can be rather slow to converge to the true parameter values and im- plementation indicates a large bias in practice. There also could be selection bias by the authors in the examples presented in their papers. In order to alleviate some of this bias, we have implemented a parametric bootstrap bias-correction method which in practice has usually produced more reasonable results than the original MSM estimates.
}
\details{
\tabular{ll}{
Package: \tab msim\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-13\cr
License: \tab GPL (>=2)\cr
}
The most important functions in the package are msm and boot.msm which will actually compute MSM estimates. The package is accompanied by a vignette (msimVignette) as well as a design document (msim_design_v3) which explains much of the theory and algorithms in the package.}
\author{
Lindsey Dietz diet0146@umn.edu
}
\references{
Jiang, J. (1998). Consistent Estimators in Generalized Linear Mixed Models. Journal of the American Statistical Association, 93, 720–729.

Jiang, J. and Zhang, W. (2001). Robust estimation in generalized linear mixed models. Biometrika, 88, 753–765.
}
\keyword{ package }
\seealso{
\code{\link{optim}}
\code{\link{multiroot}}
\code{\link{nleqslv}}
}
\examples{

}
