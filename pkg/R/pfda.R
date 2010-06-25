#pfda.R

make_c_debug<-function(...){
	debugs<-na.omit(as.integer(c(...)))
	debugs<-debugs[debugs>0]
	if(length(debugs)){
		Cdebug<-rep(0L,max(debugs))
		Cdebug[debugs]<-TRUE
		return(Cdebug)
	} else return(integer(0))
}
Cdebug<-function(...){
	cd<-make_c_debug(...)
	return(as.integer(c(length(cd),cd)))
}
pfdaControl<-function(...,penalty.method=c('AIC','CV'),minimum.variance = 1e-4, convergence.tolerance = 1e-2, max.iterations=10000, nfolds = 10,	binary.k0=100, binary.kr=10, binary.burnin = 100, nknots=11){
	structure(modifyList(list(
		penalty.method=match.arg(penalty.method),
		minimum.variance=as.double(minimum.variance),
		convergence.tolerance = as.double(convergence.tolerance),
		max.iterations=as.integer(max.iterations),
		nfolds =nfolds,
		useC=TRUE,
		binary.k0=as.integer(binary.k0), binary.kr=as.integer(binary.kr), binary.burnin=as.integer(binary.burnin),
		nknots=nknots,
		optim.method = "Nelder-Mead"
	),list(...)),class=c('pfdaControl','list'))
}
