\name{pfdaSingle}
\alias{pfdaSingle}
\alias{pfdaDual}
\alias{PairedFunctionalData}
\alias{FunctionalData}
\title{Fit single or paired functional curves to data.}
\description{
Fits curves to data estimating the mean curve(s), principle components, and the principle component scores.
}
\usage{
pfdaSingle(formula, data , obs, knots=NULL, k=NULL, delta=1e-4, penalties=NULL, tol=1e-2, MaxIter=10000, NCompTol=1/12, cv=.1,...)
pfdaDual(formula, data , obs, knots=NULL, k=NULL, delta=1e-4, penalties=NULL, tol=1e-2, MaxIter=10000, NCompTol=1/12, cv=.1,...)
}
\arguments{
  \item{formula}{A forumla of the form \code{y~t} or \code{cbind(y,z)~t}}
  \item{data}{A data frame object or environment where the variable in formula are defined.}
  \item{obs}{Factor indicating the observational unit for the oberstaions.  Indicates the individual curves to be estimated.}
  \item{knots}{The knots for the }
  \item{k}{The number of principle components. \code{pfdaDual} requires a length 2 vector for the number of principle components for Y and Z respectivly.  A value of \code{NULL} instructs to estimate the number of principle components.}
  \item{delta}{A factor included for ridge type regression in estimating initial values for the algorithm.}
  \item{penalties}{The penalties for adjusting the smoothness of the mean functions and principle components.  \code{pfdaSinlge} requires a length two vector for the mean function and the principle components, respectivly.  \code{pfdaDual} requires a length four vector with the components representing first the mean functions, Y then Z, then the principle componentns, Y then Z, respectivly.  A value of \code{NULL} specifies that the penalties should be estimated.  See details.}
  \item{tol}{The tolerance for the convergence of the algorithm.}
  \item{MaxIter}{The maximum number of itterations before terminiating with an error.}
  \item{NCompTol}{The tolerance for estimating the number of principle components. Only relevant if \code{k=NULL}}
  \item{cv}{The percentage to be withheld on each cross validation run. Only relevent if \code{penalties=NULL}}
  \item{...}{extra arguments SEE DETAILS}
}
\details{
	The \code{pfdaSingle} function fits a separate curve for each observational unit.  The curve fitted is of the form 
	\deqn{\mathbf{y_i}(t)=\mu(t)+f_i(t)}{y_i(t)=mu(t)+fi(t)}
	where \eqn{\mathbf{y_i}}{y_i} is the vector of observations for observational unit \eqn{i}.  
	\eqn{\mu(t)}{mu(t)} is the common mean function.  
	\eqn{f_i(t)}{fi(t)} is the individual deviation for observational unit \eqn{i}.
	\eqn{\mu}{mu} and \eqn{f_i} are estimated non-parametrically by funcional principle component analysis by breaking the function into a product of basis functions and parameters.
	\eqn{\mu(t)=B(t)*\theta_\mu}{mu(t)=B(t)*theta_mu} and \eqn{f_i(t)=f_1(t)*\alpha_{i,1}+\ldots+f_k(t)*\alpha_{i,k}=B(t)\Theta_f\alpha_i} where \eqn{\Theta_f} is a matrix representing the coefficients for forming the principle components from the basis functions \eqn{B(t)}, and \eqn{\alpha_i} represents the principle component scores.
	The basis functions \eqn{B(t)} are orthogonalized spline basis functions using the \code{knots} specified.  If the \code{knots} are not specified quantiles of the time/predictor variable are used. See \code{\link{OrthogonalizedSplineBasis}}.
	
	the \code{pfdaDual} function is used to fit bivariate paired curves. the model is the same as for the single model except correlation of the principle component scores is taken into account.  This covariance is represented by a slot in the output of \code{Lambda}.
	
	If \code{k} is not specified the number of principle components are estimated.  Models are fit with increasing number of principle components stopping when either the variance of the principle component is too small to continue estimation (as the algorithm is unstable for too many princiople components) or the ratio of the variances of the principle component scores for consecutive principle components is below \code{NCompTol}.  In the case of \code{pfdaDual} the number of principle components are estimated using separate (independent) models for each response Y and Z.
	
	If \code{penalties} are not specified, the penalties are estimated using \code{\link{optim}} on the cross validated log likelihood.  The cross validation is performed withholding a \code{cv} percent of the data, fitting model on the remaining data, and computing the log likelihood  of the withheld potion.  This is repeated until all the data has been withheld.  The value used in optimization is the sum of all log liklihood for the different withheld potions.  If neither \code{k} nor \code{penalties} are specified the penalties are estimated for each different number of principle components.  The penalties are also reestimated in the joint model starting from the final penalties in the individual models.  This is accomplished through an nonformal named parameter \code{optimstart}. 
}
\value{An object of class \code{\link{FDModel}} or \code{\link{PFDModel}}}
\seealso{\code{\link{FDModel}}, \code{\link{PFDModel}}, \code{\link{SplineBasis}}}
\references{
\emph{Joint modelling of paired sparse functional data using principal components}
Lan Zhou; Jianhua Z. Huang; Raymond J. Carroll
Biometrika 2008 95: 601-619
}


\keyword{smooth}
\keyword{models}
