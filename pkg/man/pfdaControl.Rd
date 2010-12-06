\name{pfdaControl}
\alias{pfdaControl}

\title{Control Setting for pfda}
\description{
This gives the settings and defaults for fitting the EM algorithms in the \code{\link{pfda}} function
}
\usage{
pfdaControl(..., penalty.method = c("AIC", "CV"), minimum.variance = 1e-04, convergence.tolerance = 0.01, max.iterations = 10000, nfolds = 10,  binary.k0 = 100, binary.kr = 10, binary.burnin = 100, nknots=11)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{any extra arguments to be passed in.}
  \item{penalty.method}{the method for optimizing penalties. Either CV for cross validation or AIC for Akaike's an information criteria.}
  \item{minimum.variance}{The minimum allowable as a variance.}
  \item{convergence.tolerance}{The tolerance for determining convergence of the EM algorithm.}
  \item{max.iterations}{The maximum number of Iterations before determining the algorithm failed.}
  \item{nfolds}{Number of folds in cross validation.}
  \item{binary.k0}{The number of simulated runs for the initial steps in the stocastic approximation involved with the binary drivers.}
  \item{binary.kr}{the ongoing number of simulation runs for the stocastic approximation for binary drivers.}
  \item{binary.burnin}{The length of the burn in period for stochasitc approximation.}
  \item{nknots}{the number of internal knots for the splines. Ignored if knots are specified.}
}
\details{
	This is a convenience function for specifying a list of control parameters that control the fit of the EM aglorithm.  With the \dots argument there are several arguments that can be additionally specified. \code{optim.method} controls the method used in \link[stats]{optim} to optimize the penalties.  The \code{optim.start} can be specified to give a starting value for the optimization of the penalties.
}
\value{
a list with the class set to pfdaControl.  Named elements are the same as the named arguments, plus any additional named arguments given in \code{\dots}
}
\author{Andrew Redd}
\seealso{
\code{\link{pfda}}
}
\examples{
pfdaControl()
}
\keyword{ utilities }
