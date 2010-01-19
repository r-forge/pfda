\name{FunctionalDataModelAnalysis}
\alias{FDModel}
\alias{PFDModel}
\alias{FDAModelParameters}
\alias{plot.FDModel}
\alias{plot.PFDModel}


\title{Using Functional Data Analysis Model Objects}

\description{
	The objects given as results of model fitting in package \link{pfda-package}{pfda}.
}
\details{
	Given as a result of either pfdaSingle or pfdaDual.  Embeds the fitting data as a \code{\link{FunctionalData}} or \code{\link{PairedFunctionalData}} object, under the slot \code{@FittingData}.  Also embeds the basis functions used in fitting the data as a \code{\link{OrthogonalSplineBasis}} object under the slot \code{@Basis}.  The parameters of the fit are stored in a \code{FDAModelParameters} object.
	
	FDModel and PFDModel both inherit from FunctionalDataModelAnalysis class for which most generic functions are defined.  All are S4 classes.  The following methods are defined.
	\begin{enumerate}
		\item \code{show}
		\item \code{logLik}
		\item \code{residuals}
		\item \code{plot}
	\end{enumerate}
	Other data can be accessed directly through the slots.
		
	\code{logLik} gives the log liklihood for the given data or for new data if given in the \code{newdata} parameter.
	
	\code{residuals} gives the residuals for the model.  Returns a vector for \code{FDModel} object and a matrix for \code{PFDModel} object, one column for y and one for z.
	
	The plot method works similar to the method for \code{\link{lm}} which preduces several diagnotic plots which can be specified by the argument \code{which}.
	
	\code{FDAModelParameters} have the following slots with the interpretations.
	\begin{enumerate}
		\item \code{theta_mu} - (numeric) holds the coefficients (weights) of the basis functions that define the mean function.
		\item \code{Theta_f} - (matrix)  each column defined the coefficients for the principle components.
		\item \code{sigma} - (numeric) estimated variance of the error
		\item \code{Alpha} - (matrix) Matrix of principle component scores.  Each row represents the weights for an observational unit.
		\item \code{Da} - (numeric) variance(s) of the principle component scores.  The variance(s) are given as a vector of the diagonal elements of variance matrix, which is defined to be diagonal.
		\item \code{npc} - (integer) the number of the princple components.
		\item \code{penalties} - (numeric) the penalties that were used in the final fit.
		\item \code{Sigma} - (array) An array of the conditional variances of Alpha given the data.
		\item \code{Ry} - (numeric) vector fo the residuals
	\end{enumerate}
	
	\code{FDModel} has a \code{FDModelParameters} object named \code{parameters}.  The \code{PFDModel} has two parameters slots named \code{Analysis_y} and \code{Analysis_z} for the y and z components of the response.  In addition there are the slots in \code{PFDModel} that store the information about the covariance of the principle component scores.  Slot \code{Lambda} stores the covariance between the scores for Y and Z, and \code{Simga_yz} gives the conditional covariance of the scores given the data.
}
\seealso{
	\code{\link{pfdaSingle}},\code{\link{pfdaDual}},\code{\link{SplineBasis}},\code{\link{FunctionalData}},\code{\link{pfda-package}},
}