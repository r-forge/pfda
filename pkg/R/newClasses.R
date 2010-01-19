#  newClasses.R

# FunctionalDataModelAnalysis Class
# setClass("FunctionalDataModelAnalysis", representation("VIRTUAL",Basis = "OrthogonalSplineBasis", ConvergenceCriteria="numeric", iterations ="integer", FittingData = "FunctionalData", CVlogLik="numeric"))
# setMethod("show","FunctionalDataModelAnalysis",function(object){
	# cat("Functional data model analysis, using an orthogonalized spline basis as principle components\n")
	# cat("Number of subjects:",  nlevels(object@FittingData),"\n")
 # } )
 

# pfdaLinearModelParameters
# setClass("pfdaLinearModelParameters",representation(theta_mu="numeric", Theta_f="matrix",sigma = "numeric", Alpha="matrix", Da = "numeric", npc="integer",penalties="numeric", Sigma="array", Ry="numeric"))
# setValidity("pfdaLinearModelParameters",function(object){
	# k=object@npc
	# if( NCOL(object@Alpha)!=k || length(object@Da) !=k || NCOL(object@Theta_f)!=k || any(dim(object@Sigma)[1:2] != k) ) return('Incompatible dimentions found')
	# TRUE
# })
# setMethod("show","pfdaLinearModelParameters",function(object){
	# cat("Penalties", object@penalties,"\n")
	# cat("Residual Variance: ", object@sigma, "\n")
	# cat("Da", object@Da,"\n")
	# cat("Number of principle components:", object@npc,"\n")
	# d<-data.frame(theta_mu=object@theta_mu,theta_f=object@Theta_f)
	# show(d)
# }) 

#pfdaBinaryModelParameters
# setClass("pfdaBinaryModelParameters",representation("pfdaLinearModelParameters"))
# setValidity("pfdaBinaryModelParameters",function(object){
	# cat(object@sigma)
	# if(length(object@sigma))return("value in sigma slot,should have length zero.")
	## if(!(is.na(object@sigma) || is.null(object@sigma)))return("value in place of sigma slot should be NA")
	# if(length(object@Ry))return("cannot assign value to residuals (Ry), for binary models.")
	# TRUE
# })

# setClass("pfdaUnivariateModel", representation("VIRTUAL","FunctionalDataModelAnalysis"))
# setClass("pfdaLinearModel", representation("pfdaUnivariateModel", parameters="pfdaLinearModelParameters"))
# setClass("pfdaBinaryModel", representation("pfdaUnivariateModel", parameters="pfdaBinaryModelParameters"))

# setClass("pfdaBivariateModel", representation("VIRTUAL","FunctionalDataModelAnalysis", Sigma_yz = "array", Lambda="matrix"))
# setClass("pfdaLinearLinearModel",representation("pfdaBivariateModel", Ay='pfdaLinearModelParameters',Az='pfdaLinearModelParameters'))
# setClass("pfdaLinearBinaryModel",representation("pfdaBivariateModel", Ay='pfdaLinearModelParameters',Az='pfdaBinaryModelParameters'))
# setClass("pfdaBinaryBinaryModel",representation("pfdaBivariateModel", Ay='pfdaBinaryModelParameters',Az='pfdaBinaryModelParameters'))

# FDModel
# setGeneric("penalty",function(object,...)standardGeneric("penalty"))
# setMethod("penalty",signature("pfdaUnivariateModel"),function(object,...){
	# p<-object@parameters@penalties
	# names(p)<-c('mean','principle.components')
	# return(p)
# })
# setMethod("show","pfdaUnivariateModel",function(object){
	# callNextMethod()
	# show(object@parameters)
# })
# plot.Univariate<-function(x, y, ...){
	# dots<-list(...)
	# obase<-x@Basis
	# xlab<- if(hasArg(xlab)) dots$xlab else "t"
	# ylab<- if(hasArg(xlab)) dots$xlab else "y"
	# whichplot <-if(hasArg(which)) dots$which else 1:4
	
	# tvalues<-seq(min(x@FittingData@t), max(x@FittingData@t),length.out=100)
	# B<-evaluate(obase, tvalues)
	
	## if(length(whichplot)>1){
		## s<-seq_along(whichplot)
		## layout(matrix(s, ncol=round(sqrt(length(whichplot)))))
	## }
	
	# for(p in whichplot){
		# if(p==1){ # mu and nu plot
			# mu<-B %*% x@parameters@theta_mu
			# plot(tvalues,mu, xlab=xlab, ylab=ylab, main=expression(mu(t)),type='l')
		# } 
		# else if(p==2) { # principle components
			# f<-B %*% x@parameters@Theta_f
			# matplot(tvalues, f, xlab=xlab, ylab=ylab, main=expression(f[i](t)),type='b')
		# } 
		# else if(p==3) { # residual Q-Q plot
			# qqnorm(residuals(x), main="Residual Q-Q Plot")
		# }
		# else if(p==4) { # Y curves
			# levels.subset <- if(hasArg(levels.subset)) dots$levels.subset else levels(x@FittingData)
			# stopifnot(length(levels.subset)>0)
			# mu<-B %*% x@parameters@theta_mu
			# f<-B %*% x@parameters@Theta_f
			# apply(x@parameters@Alpha[levels(x@FittingData)%in% levels.subset,,drop=FALSE], 1, function(Alpha){
				# mu+f%*%Alpha
			# })->yplot
			# cols<-rainbow(nlevels(x@FittingData))
			# index<-x@FittingData@obs %in% levels.subset
			# plot(x@FittingData@t[index], x@FittingData@y[index], col=cols[as(factor(x@FittingData@obs[index]),"integer")],xlab=xlab,ylab=ylab,main="Fitted Curves")
			# matlines(tvalues,yplot,col=cols,type='l')
		# }
	# }
# }
# setMethod("plot",signature("pfdaUnivariateModel", "missing"),plot.Univariate)

# logLik.pfdaLinearModel<-function(object,...){
	# dots<-list(...)
	# data<-as(if(hasArg(newdata)) dots$newdata else object@FittingData,"FunctionalData")
	# t<-data@t
	# y<-data@y
	# o<-data@obs[drop=TRUE]
	# B<-evaluate(object@Basis,t)
	# BTf<-B%*%object@parameters@Theta_f
	# L<-numeric(1)
	# for(i in seq_len(nlevels(o))){
		# ind<-i==as(o,"integer")
		# Sy<-BTf[ind,,drop=FALSE]%*%tcrossprod(diag(x=object@parameters@Da,ncol=object@parameters@npc),BTf[ind,,drop=FALSE])+diag(object@parameters@sigma,nrow=sum(ind))
		# Ry<-y[ind]-B[ind,]%*%object@parameters@theta_mu
		# D<-determinant(Sy,logarithm=TRUE)
		# if(D$sign<0)stop("Negative determinant encountered for variance matrix.")
		# L<-L+as.numeric(D$modulus+crossprod(Ry,solve(Sy,Ry)))
	# }
	## length(object@predicted)*log(sigma)+1/sigma*crossprod(residuals(object))
	# L
# }
# setMethod("logLik","pfdaLinearModel",logLik.pfdaLinearModel)
# setMethod("residuals","pfdaLinearModel",function(object,...)return(object@parameters@Ry))

# setMethod("penalty",signature("pfdaBivariateModel"),function(object,...){
	# py<-object@Analysis_y@penalties
	# pz<-object@Analysis_z@penalties
	# p<-rbind(py,pz)
	# colnames(p) <- c('mean','principle.components')
	# rownames(p) <- colnames(object@FittingData[[1]])
	# return(p)
# })
# plot.Bivariate<-function(x, y, ...){
	# dots<-list(...)
	# opar<-par(no.readonly=TRUE)
	# obase<-x@Basis
	# xlab<- if(hasArg(xlab)) dots$xlab else "t"
	# ylab<- if(hasArg(ylab)) dots$ylab else "y"
	# zlab<- if(hasArg(ylab)) dots$zlab else "z"
	# whichplot <-if(hasArg(which)) dots$which else 1:6
	
	# tvalues<-seq(min(x@FittingData@t), max(x@FittingData@t),length.out=100)
	# B<-evaluate(obase, tvalues)
	
	## if(length(whichplot)>1){
		## s<-seq_along(whichplot)
		## layout(matrix(s, ncol=round(sqrt(length(whichplot)))))
	## }
	# for(p in whichplot){
		# if(p==1){ # mu and nu plot
			# mu<-B %*% x@Analysis_y@theta_mu
			# nu<-B %*% x@Analysis_z@theta_mu
			# ycol<-heat.colors(1)
			# zcol<-topo.colors(1)
			# par(mar=c(5, 4, 4, 4) + 0.1)
			# plot(tvalues,mu, xlab=xlab, ylab="", main="Mean functions",type='l', col=ycol)
			# mtext(expression(mu(t)),side=2,line=2.5, col=ycol)
			# usr<-par("usr")
			# ratio<-diff(usr[3:4])/diff(range(mu))
			# usr[3:4]<-range(nu)+c(-1,1)*(ratio-1)*diff(range(nu))/2
			# par(usr=usr)
			# Axis(x=range(nu),side=4)
			# mtext(expression(nu(x)),side=4,line=2.5, col=zcol)
			# lines(tvalues,nu, col=zcol)
		# } 
		# else if(p==2) { # principle components
			# par(mar=c(5, 4, 4, 4) + 0.1)			
			# ycols<-if(hasArg(ycols)) dots$ycols else heat.colors(x@Analysis_y@npc)
			# zcols<-if(hasArg(zcols)) dots$zcols else topo.colors(x@Analysis_z@npc)
			# f<-B %*% x@Analysis_y@Theta_f
			# g<-B %*% x@Analysis_z@Theta_f
			# matplot(tvalues, f, xlab=xlab, ylab="", main="Principle Components",type='l', col=ycols)
			# usr<-par("usr")
			# ratio<-diff(usr[3:4])/diff(range(f))
			# usr[3:4]<-range(g)+c(-1,1)*(ratio-1)*diff(range(g))/2
			# par(usr=usr)
			# Axis(x=range(g),side=4)
			# matlines(tvalues, g,type='l',col=zcols)
			# mtext(expression(f[i](t)), side=2, line=2.5, col=ycols[1])
			# mtext(expression(g[i](t)), side=4, line=2.5, col=zcols[1])
			##TODO
			##add overplotting points and intervals
		# } 
		# else if(p==3) { # residual Q-Q plot for Y
			# qqnorm(residuals(x)[,1], main="Y Residual Q-Q Plot")
		# }
		# else if(p==4) { # residual Q-Q plot for z
			# qqnorm(residuals(x)[,2], main="Z Residual Q-Q Plot")
		# }
		# else if(p==5) { # Y curves
			# levels.subset <- if(hasArg(levels.subset)) dots$levels.subset else levels(x@FittingData)
			# stopifnot(length(levels.subset)>0)
			# mu<-B %*% x@Analysis_y@theta_mu
			# f<-B %*% x@Analysis_y@Theta_f
			# apply(x@Analysis_y@Alpha[levels(x@FittingData)%in% levels.subset,,drop=FALSE], 1, function(Alpha){
				# mu+f%*%Alpha
			# })->yplot
			# cols<-rainbow(nlevels(x@FittingData))
			# index<-x@FittingData@obs %in% levels.subset
			# plot(x@FittingData@t[index], x@FittingData@y[index], col=cols[as(factor(x@FittingData@obs[index]),"integer")],xlab=xlab,ylab=ylab,main="Fitted Curves")
			# matlines(tvalues,yplot,col=cols,type='l')
		# }		
		# else if(p==6) { # Z curves
			# levels.subset <- if(hasArg(levels.subset)) dots$levels.subset else levels(x@FittingData)
			# stopifnot(length(levels.subset)>0)
			# mu<-B %*% x@Analysis_z@theta_mu
			# f<-B %*% x@Analysis_z@Theta_f
			# apply(x@Analysis_z@Alpha[levels(x@FittingData)%in% levels.subset,,drop=FALSE], 1, function(Alpha){
				# mu+f%*%Alpha
			# })->yplot
			# cols<-rainbow(nlevels(x@FittingData))
			# index<-x@FittingData@obs %in% levels.subset
			# plot(x@FittingData@t[index], x@FittingData@z[index], col=cols[as(factor(x@FittingData@obs[index]),"integer")],xlab=xlab,ylab=zlab,main="Fitted Curves")
			# matlines(tvalues,yplot,col=cols,type='l')
		# }
	# }
	# par(opar)
# }
# setMethod("plot",signature("pfdaBivariateModel","missing"),plot.Bivariate)
# setMethod("show","pfdaBivariateModel", function(object){
	# callNextMethod()
	# cat("Parameters for y:\n")
	# show(object@Analysis_y)
	# cat("Parameters for z:\n")
	# show(object@Analysis_z)
# })


# logLik.pfdaLinearLinearModel<-function(object,...){
	# dots<-list(...)
	# data<-if(hasArg(newdata)) dots$newdata else object@FittingData
	# if(class(data)=="PairedFunctionalData"){
		# t<-data@t
		# y<-data@y
		# z<-data@z
		# o<-data@obs[drop=T]
	# } else if(class(data)=="data.frame"){
		# stopifnot(NCOL(data)==3)
		# stopifnot(NCOL(data[[1]])==2)
		# stopifnot(is.factor(data[[3]]))
		# t<-data[[2]]		
		# y<-data[[1]][,1,drop=T]
		# z<-data[[1]][,2,drop=T]
		# o<-data[[3]][drop=T]
	# } else stop("Unsupported class for new data")
	# B<-evaluate(object@Basis,t)
	# ka<-object@Analysis_y@npc
	# kb<-object@Analysis_z@npc
	# BTf<-B%*%object@Analysis_y@Theta_f
	# BTg<-B%*%object@Analysis_z@Theta_f
	# seps<-object@Analysis_y@sigma
	# sxi<-object@Analysis_z@sigma
	# Da<-object@Analysis_y@Da
	# Db<-object@Analysis_z@Da
	# L<-object@Lambda
	# ll<-0
	# for(i in seq_len(nlevels(o))){
		# ind<-i==as(o,"integer")
		# ni<-sum(ind)
		# Sy<-BTf[ind,,drop=FALSE]%*%diag(x=Da,ka,ka)%*%t(BTf[ind,,drop=FALSE])+diag(seps,sum(ind))
		# Sz<-BTg[ind,,drop=FALSE]%*%diag(x=Db,kb,kb)%*%t(BTg[ind,,drop=FALSE])+diag(sxi,sum(ind))
		# Syz<-BTf[ind,,drop=FALSE]%*%diag(Da,ka,ka)%*%t(L)%*%t(BTg[ind,,drop=FALSE])
		# Ry<-y[ind]-B[ind,]%*%object@Analysis_y@theta_mu
		# Rz<-z[ind]-B[ind,]%*%object@Analysis_z@theta_mu
		
		# Cmat<-rbind(cbind(Sy,Syz),cbind(t(Syz),Sz))
		# Cimat<-solve(Cmat)
		# Sy.inv <-Cimat[seq_len(ni),seq_len(ni)]			#solve(Sy-crossprod(Syz, solve(Sz,t(Syz))))
		# Sz.inv <- Cimat[seq_len(ni)+ni,seq_len(ni)+ni]			#Szi+Szi%*%t(Syz)%*%Sy.inv%*%Syx%*%Szi		#solve(Sz-tcrossprod(Syz,solve(Sy,t(Syz))))
		# Syz.inv<- Cimat[seq_len(ni),seq_len(ni)+ni]			#-solve(Sy, Syz)%*%Sz.inv		
		
		# Dety<-determinant(Sy,logarithm=TRUE)
		# Detzinv<-determinant(Sz.inv,logarithm=TRUE)
		# if((Dety$sign<0)|(Detzinv$sign<0))stop("Negative determinant encountered for variance matrix.")
		# ll<- ll + Dety$modulus-Detzinv$modulus
		# ll<- ll + crossprod(Ry,Sy.inv%*%Ry)
		# ll<- ll + crossprod(Rz,Sz.inv%*%Rz)
		# ll<- ll + 2*crossprod(Ry,Syz.inv%*%Rz)
	# }
	# ll
# }
# setMethod("logLik", "pfdaLinearLinearModel", logLik.pfdaLinearLinearModel)
# setMethod("residuals","pfdaLinearLinearModel", function(object){
	# cbind(resid.y=object@Analysis_y@Ry,resid.z=object@Analysis_z@Ry)
# })

