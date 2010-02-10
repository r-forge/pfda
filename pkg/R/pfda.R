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
pfdaControl<-function(...,penalty.method=c('AIC','CV'),minimum.variance = 1e-4, convergence.tolerance = 1e-2, 
	max.iterations=10000, pc.tolerance=1/12, nfolds = 10, trace=TRUE,
	useC=TRUE,binary.k0=100, binary.kr=10, C.debug=Cdebug()){
	list(
		penalty.method=match.arg(penalty.method),
		minimum.variance=as.double(minimum.variance),
		convergence.tolerance = as.double(convergence.tolerance),
		max.iterations=as.integer(max.iterations),
		pc.tolerance=as.double(pc.tolerance),
		nfolds =nfolds,
		# expected.value.tolerance =expected.value.tolerance,
		trace=trace,useC=useC,
		binary.k0=binary.k0, binary.kr=binary.kr,
		C.debug=as.integer(C.debug),
		...
	)
}
