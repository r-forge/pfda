pfdaParseFormula<-function(formula, data, envir=environment(formula)){
	expr2char <- function(x) paste(deparse(x), collapse = "")
	if(length(formula)==3){ 
		op<-formula[[1]]
		if(op=='~'){
			lhs<-if(length(formula)==3) Recall(formula[[2]],data) else NULL
			attr(lhs,'pfda.role')<-'response'
			rhs<-Recall(formula[[length(formula)]],data)
			if(is(rhs,"pfda.model.frame"))
				pfda.model.join(addto = rhs, response = lhs)
			else if(attr(rhs,'pfda.role')=='splinegroup')
				pfda.model.frame.new(response = lhs, splinegroup = rhs)
			else stop("invalid model")
		} else 
		if(op=='%&%') {
			first  <- Recall(formula[[2]],data)
			second <- Recall(formula[[3]],data)
			bound<-data.frame(first,second)
			attr(bound,'pfda.role')<-'bound'
			bound
		} else 
		if(op=='%|%') {
			variable  <- Recall(formula[[2]], data)
			attr(variable,'pfda.role')<-'splinevariable'
			subjectID <- Recall(formula[[3]], data)
			attr(subjectID,'pfda.role')<-'subjectID'
			if(!is.factor(subjectID[[1]]))
				stop(gettextf("The variable `%s` on the right hand side of `%%|%%` should be a factor.",expr2char(formula[[3]])))
			vardata<-data.frame(variable,subjectID)
			attr(vardata,"pfda.role")<-"splinegroup"
			vardata
		} else 
		if(op=='+') {
			var1<-Recall(formula[[2]],data)
			var2<-Recall(formula[[3]],data) 
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
		structure(eval(formula,data,envir),name=expr2char(formula))
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
infer.driver<-function(mf){
	noadds<-function(c)if(!is.null(mf$additive) && NCOL(mf$additive))stop(gettextf("Inferred class '%s' does not allow additive variables.",c))
	stopifnot(is(mf,'pfda.model.frame'))
	stopifnot(valid.pfda.model.frame(mf))
	if(NCOL(mf$response)==1){ # single or additive
		if(NCOL(mf$splinegroup)==2) { # single
			if(is(mf$response,'factor')) {
				noadds('single.binary')
				if(nlevels(mf$response[[1]])) return('single.binary') else stop("Too many factor levels.")
			} else if(is(mf$response[[1]],'logical') || (is(mf$response[[1]],'integer')&&length(unique(mf$response[[1]]))==2)) {
				noadds('single.binary')
				return('single.binary')
			} else if(is(mf$response[[1]],'numeric')) {
				return('single.continuous')
			} else stop('Unusable response class')
		} else if(NCOL(mf$splinegroup)==3) { #additive
			return('additive')
		} else stop('Bad spline group')
	} else if(NCOL(mf$response)==2) {
		noadds('dual')
		y = mf$response[[1]]
		z = mf$response[[2]]
		if(is(z,'factor')|| is(z,'logical')) stop(gettextf("Second variable in dual models is not allowed to be a '%s'",class(z)))
		else if(is(z,'numeric')) {
			if(is(y,'factor')){
				if(nlevels(y)!=2) stop(gettextf("too many levels in %s, you might need to use drop.levels",names(mf$response)[1]))
				else return('dual.mixed')
			} else if(is(y,'logical') || (is(y,'integer')&&length(unique(y))==2)) return('dual.mixed')
			else if(is(y,'numeric')) return('dual.continuous')
			else stop(gettextf("bad class, '%s', for variable %s",class(y),names(mf$response)[1]))
		} else stop(gettextf("bad class, '%s', for variable %s",class(z),names(mf$response)[2]))
	} else stop("Unusable number of response variables.")
}

