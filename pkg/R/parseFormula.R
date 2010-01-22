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
pfdaParseFormula<-function(formula, data=environment(formula)){
	expr2char <- function(x) paste(deparse(x), collapse = "")
	if(length(formula)==3){ 
		op<-formula[[1]]
		if(op=='~'){
			lhs<-if(length(formula)==3) Recall(formula[[2]],data=data) else NULL
			attr(lhs,'pfda.role')<-'response'
			rhs<-Recall(formula[[length(formula)]],data=data)
			if(is(rhs,"pfda.model.frame"))
				pfda.model.join(addto = rhs, response = lhs)
			else if(attr(rhs,'pfda.role')=='splinegroup')
				pfda.model.frame.new(response = lhs, splinegroup = rhs)
			else stop("invalid model")
		} else if(op=='%&%') {
			first  <- Recall(formula[[2]],data=data)
			second <- Recall(formula[[3]],data=data)
			bound<-data.frame(first,second)
			attr(bound,'pfda.role')<-'bound'
			bound
		} else if(op=='%|%') {
			variable  <- Recall(formula[[2]], data=data)
			attr(variable,'pfda.role')<-'splinevariable'
			subjectID <- Recall(formula[[3]], data=data)
			attr(subjectID,'pfda.role')<-'subjectID'
			if(!is.factor(subjectID[[1]]))
				stop(gettextf("The variable `%s` on the right hand side of `%%|%%` should be a factor.",expr2char(formula[[3]])))
			vardata<-data.frame(variable,subjectID)
			attr(vardata,"pfda.role")<-"splinegroup"
			vardata
		} else if(op=='+') {
			var1<-Recall(formula[[2]],data=data)
			var2<-Recall(formula[[3]],data=data) 
			if((is.null(attr(var1,'pfda.role')) || attr(var1,'pfda.role')!='splinegroup') && (is.null(attr(var2,'pfda.role')) || attr(var2,'pfda.role')!='splinegroup')) {
				if(is(var1,'pfda.model.frame'))
					pfda.model.join(addto = var1, additive = var2)
				else if(is(var2,'pfda.model.frame'))
					pfda.model.join(addto = var2, additive = var1)
				else
					data.frame(var1,var2)
			} else if(attr(var1,'pfda.role')=='splinegroup') {
				pfda.model.frame.new(splinegroup=var1,additive = var2)
			} else if(attr(var2,'pfda.role')=='splinegroup') {
				pfda.model.frame.new(splinegroup=var2, additive = var1)
			} else stop("how the hell did you get here.")
		} else {
			stop(paste("I don't know what to do with this operator ",op))
		} 
	} else { #individual variables
		x = substitute(~m, list(m = formula))
		model.frame(x,data=data)
	}
}
traverseFormula<-function(x,level=0){
	switch(class(x),
		name = { replicate(level, cat(". ")); cat(x,'\n')},
		formula = { traverseFormula(x[[2]],level=level+1); replicate(level, cat(". ")); cat('~\n'); traverseFormula(x[[3]],level=level+1)},
		call = { replicate(level, cat(". ")); cat(x[[1]],'\n'); lapply(x[-1],traverseFormula,level=level+1)}
	)
	invisible(NULL)
}
pfda.model.frame.new<-function(response=NULL,splinegroup=NULL,additive=NULL){
	r<-list(response=response,splinegroup=splinegroup,additive=additive)
	class(r)<- c('list','pfda.model.frame')
	return(r)
}
pfda.model.join<-function(addto = NULL, response = NULL, splinegroup = NULL, additive = NULL){
	if(!is.null(addto))stopifnot(is(addto,'pfda.model.frame'))
	if(!is.null(response))stopifnot(is(response,'data.frame'))
	if(!is.null(splinegroup))stopifnot(is(splinegroup,'data.frame'))
	if(!is.null(additive))stopifnot(is(additive,'data.frame'))
	if(!is.null(addto)){
		if(valid.pfda.model.frame(addto)){
			if(!is.null(response)){
				if(is.null(addto$response)) addto$response = response else stop("this object already has a response variable.")
			}
			if(!is.null(splinegroup)){
				if(is.null(addto$splinegroup)) addto$splinegroup = splinegroup else stop("this object already has a spline group.")
			}
			if(!is.null(additive)){
       	if(is.null(addto$additive)) addto$additive = additive else addto$additive = data.frame(addto$additive, additive)
			}
		} else stop("not a valid object to add to.")
	} else {
		addto <- pfda.model.frame.new(response = response, splinegroup = splinegroup, additive = additive)
	}
	addto
}
valid.pfda.model.frame<-function(object){
	if(!("pfda.model.frame" %in% class(object))) return(FALSE)
	if(length(object) != 3) return(FALSE)
	if(!all(c('response','splinegroup','additive') %in% names(object))) return(FALSE)
	TRUE	
}
valid.pfda.splinegroup<-function(object){
	if(!is(object,'data.frame')) return(FALSE)
	if(length(object)!= 2 && length(object)!=3)return(FALSE)
	if(!is(object[,NCOL(object)],'factor'))return(FALSE)
	TRUE
}
