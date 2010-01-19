# FunctionalData_class.R

setClass("FunctionalData",representation('data.frame'))
setValidity("FunctionalData",function(object){
	if(length(object)!=3)return(gettextf('incorrect length. Needs 3 has %d',length(object)))
	if(class(object[[1]])!='data.frame')return(gettextf('incorrect response class, needs `data.frame` given %s',class(object[[1]])))
	if(!is.numeric(object[[2]])) return(gettextf('domain not numeric'))
	if(!is.factor(object[[3]])) return(gettextf("subject ID variable not a factor (it's a %s)",class(object[[3]])))
	TRUE
	})
setMethod("nlevels","FunctionalData",function(x)nlevels(x[[3]]))
setMethod("levels","FunctionalData",function(x)levels(x[[3]]))
setMethod("show","FunctionalData",function(object){
	cat("Functional Data\n")
	cat("Subject IDs:",levels(object),"\n\n")
	callNextMethod()
})
setMethod("plot",signature("FunctionalData","missing"),function(x,y,...){
	plotdata=data.frame(x[[1]],x[2])
	subjects=as.integer(x[[3]])
	pairs(plotdata,pch=subjects,col=subjects)
 } )
