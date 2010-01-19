#' Formula paring for pfda package that supports the & and | symbols.
#' 
#' Designed to correctly structure a data frame for the pfda functions
#' 
#' @param formula a formula object
#' @param data a data frame or envirnoment
#' @param ... passed to model.frame
#' 
#' @return a data.frame object with three elements, correctly formated for the FunctionalData class.
#' 		The first element should be a data.frame objects. The motivation for this being that we need 
#'		to preserve classes of integer, logical, and factor variables to determine appropriate methods.
#' 		In Second position is the domain variable.  Third should be a factor and defines the subject grouping.
#' 
#' 
#' 
#' @callgraph
#' @seealso subset, formula, model.frame
#' 
pfdaParseFormula<-function(formula, data=environment(formula), ...)
{
	expr2char <- function(x) paste(deparse(x), collapse = "")
	if(length(formula)==3){ 
		op<-formula[[1]]
		if(op=='~'){
			lhs<-if(length(formula)==3) Recall(formula[[2]],data=data,...) else NULL
			rhs<-Recall(formula[[length(formula)]],data=data,...)
			
			mm<-data.frame(NA,rhs)
			mm[[1]]<-lhs
			colnames(mm)[1]<-expr2char(formula[[2]])
			mm
		} else if(op=='&') {
			first  <- Recall(formula[[2]],data=data,...)
			second <- Recall(formula[[3]],data=data,...)
			data.frame(first,second)
		} else if(op=='|') {
			variable  <- Recall(formula[[2]], data=data, ...)
			subjectID <- Recall(formula[[3]], data=data, ...)
			#if((NCOL(variable)>1)||(NCOL(subjectID)>1))stop("only 1 variable allowed on each side of `|` in formula.")
			if(!is.factor(subjectID[[1]]))
				stop(gettextf("The variable `%s` on the right hand side of `|` should be a factor.",expr2char(formula[[3]])))
			vardata<-data.frame(variable,subjectID)
		} else {
			stop(paste("I don't know what to do with this operator ",op))
		} 
	} else { #individual variables
		x = substitute(~m, list(m = formula))
		model.frame(x,data=data,...)
	}
}
