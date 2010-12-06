\name{pfda}
\alias{pfda}
\title{Principal Components for Functional Data Analysis}
\description{Functional data analysis using functional principle components.  The principal components are formed using orthogonalized spline basis.
}
\usage{
pfda(model, data = environment(model), ..., driver)
}
\arguments{
  \item{model}{a formula describing the model.}
  \item{data}{optional environment with variables in \code{model}.}
  \item{\dots}{optional arguments, see details.}
  \item{driver}{argument to force a driver function, see details for options.}
}
\details{
	This is a generic function that interprets the formula into a principal component model.  The \code{driver} argument controls which function the data is passed to, and is inferred if omitted.  This will be most useful when forcing a binary response model over a continuous.
	
	The \code{model} argument is a formula object that describes the variable and model that should be estimated.  It uses unique special formula operators \code{\%&\%} and \code{\%|\%}. The following table gives the name of the driver an example formula, and any notes.
\tabular{lll}{
	\strong{Name} \tab \strong{Formula} \tab \strong{Notes} \cr
	single.continuous \tab \code{y ~ t \%|\% id + z1 + z2} \tab A single continuous response variable, \code{x}, as a function of \code{t}. \cr
	single.binary \tab \code{y ~ t \%|\% id} \tab \code{y} should be a logical or factor. \cr
	dual.continuous \tab \code{y \%&\% z ~ t \%|\% id} \tab Paired model where both responses are continuous \cr
	dual.continuous \tab \code{y \%&\% z ~ t \%|\% id} \tab Paired model where \code{y} is a binary (logical or factor) variable \cr
	additive \tab \code{y ~ t \%&\% x \%|\% id + z1 + z2} \tab Paired model where \code{y} is a binary (logical or factor) variable \cr
}
where
	\itemize{
		\item{\code{y} denotes a response variable. If this is a binary variable this must be of class \code{\link[base]{logical}} or a 2 level \code{\link[base]{factor}}.  If the binary drivers are used \code{y} must be coercible into a binary form. }
		\item{\code{z} is the second variable in dual/paired models, and is always continuous.}
		\item{\code{t} the domain variable of the functions.}
		\item{\code{x} A second domain variable for the additive model.}
		\item{\code{id} is a \code{factor} variable that specifies the subject IDs.}
		\item{\code{z1, z2, ...} are extra additive effect variables, only applies to the additive and the single continuous models.}
	}

The \code{...} argument is to hold the optional arguments that control the fit.  All optional arguments have default values that are implied if the argument is omitted.  The possible arguments are
\itemize{
	\item{\code{knots} The knots for the basis functions.  For the additive model the knots argument should be a list of length 2.  The list can be named with the names the same as those in the m\code{model}, but otherwise are taken in order.}
	\item{\code{penalties} The penalties for the model.  Takes precedence over the \code{df} argument.}
	\item{\code{df} The penalties for the curves specified in terms of target degrees.}
	\item{\code{k} the number of principal components.  Should be a vector of length 2 for dual and additive models.}
	\item{\code{control} List of arguments for controling the fit of the EM algorith.  See \code{\link{pfdaControl}} for details.}
}
The \code{penalties} and \code{df} arguments should be vectors of length two for the single models the first being the penalty for the mean curve the second for the principal components.  For the dual and additive models the arguments can be a vector or matrix the first two elements for the mean curves and the sencond for the principal components.  This translates to, in terms of a 2 by 2 matrix, the rows for the variables and the columns for mean and principal component penalties.
	
}
\seealso{
\code{\link{pfda}}
}

\references{
\emph{Joint modelling of paired sparse functional data using principal components}
Lan Zhou; Jianhua Z. Huang; Raymond J. Carroll,
Biometrika, 2008 95: 601-619
}
\author{
	Andrew Redd
	Maintainer: Andrew Redd <aredd at stat.tamu.edu>
}

\keyword{package}
