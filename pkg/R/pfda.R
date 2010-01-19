#pfda.R

# pfdaSingle<-function(formula, data , obs,  knots=NULL, k=NULL, delta=1e-4, penalties=NULL, tol=1e-2, MaxIter=10000,NCompTol=1/12, cv=.1,...){
	## if k is not specified then the number of principle components will be infered.
	## if the penalties are not specified the wil be chosen by crossvalidated log likelihood.
	## recursion for k takes precidence over recursion for penalties.
	# dots<-list(...)
	# if(is.null(dots$Cdebuglevel))Cdebuglevel=0L else Cdebuglevel<-dots$Cdebuglevel
	# if(is.null(dots$Rdebuglevel))Rdebuglevel=0L else Rdebuglevel<-dots$Rdebuglevel
	# if(any(Rdebuglevel<0))cat("entering pfdaSingle(r).\n")
	# mf<-cl<-match.call()
    # m <- match(c("formula", "data", "subset", "obs", "na.action"), names(mf), nomatch=0L)
    # mf <- mf[c(1L, m)]
    # mf$drop.unused.levels <- TRUE
    # mf[[1L]] <- as.name("model.frame")
    # mf <- eval(mf, parent.frame())
	# if(length(attr(mf,"variables"))>2)stop("Only one predictor can be specified")
	# obs=as.factor(mf[[3]])
	# obs=obs[drop=TRUE]
	# oo<-order(obs)
	# mf<-mf[oo,];obs<-obs[oo]  # ensure that the data is sorted properly
	# y<-model.response(mf)
	# if(NCOL(y)>1)stop("Multiple responses specified")
	# t<-mf[,NCOL(y)+1]
	# N<-nlevels(obs)
	# if(is.null(dots$Basis)){ #handles generation of Basis Object
		# splineorder<-if(hasArg(splineorder)) dots$splineorder else 4
		# if(is.null(knots)) { 
			# knots<-quantile(t,seq(0,1,length=11))[c(rep(1,splineorder),seq(from=2,to=10),rep(11,splineorder))]
		# }
		# if(any(Rdebuglevel<0))cat("Computing new Basis...\n")
		# obase<-OBasis(knots,order=splineorder)
	# } else {
		# if(any(Rdebuglevel<0))cat("Using Given Basis...\n")
		# obase<-dots$Basis
	# }
	# if(is.null(dots$Cdebuglevel))Cdebuglevel=0L else Cdebuglevel<-dots$Cdebuglevel
	# if(is.null(dots$Rdebuglevel))Rdebuglevel=0L else Rdebuglevel<-dots$Rdebuglevel
	# if(cv<0||cv>.5)stop("cv parameter must be between 0 and .5")
	
	# if(is.null(k)){ #estimate the number of principle components
		# k=1
		# if(any(Rdebuglevel<0))cat("Finding k=",k,"\n")
		# cl.k<-cl
		# cl.k$k=k
		# results.k <- eval(cl.k, parent.frame())			# pfdaSingle(mf,knots=knots, k=k, delta=delta, penalties=penalties, tol=tol, MaxIter=MaxIter, NCompTol=NCompTol,...)
		# if(results.k@parameters@Da[1]<=delta)stop("Estimated variance for D_alpha is smaller than minimum variance (delta). Try reducing delta.")
		# tryCatch(while(TRUE) {
			# l=k+1
			# if(any(Rdebuglevel<0))cat("Finding l=",l,"\n")
			# cl.l<-cl
			# cl.l$k<-l
			# if(is.null(penalties) && is.null(dots$optimstart))if(is.null(dots$UseCaryoverPenalties) || dots$UseCaryoverPenalties) cl.l$optimstart=results.k@parameters@penalties
			# results.l<-eval(cl.l, parent.frame())			#pfdaSingle(mf,knots=knots, k=l, delta=delta, penalties=penalties, tol=tol, MaxIter=MaxIter, NCompTol=NCompTol,...)
			# Dc<-(results.k@parameters@Da-results.l@parameters@Da[1:k])/results.k@parameters@Da
			# if((all(Dc<NCompTol)&(results.l@parameters@Da[l]<NCompTol*results.l@parameters@Da[k]))||results.l@parameters@Da[l]==delta) break else {
				# results.k <- results.l
				# k <- l
			# }
		# } ,error=function(e)warning(paste("Error encountered computing estimates for k =",l,"number of principle components.  Recovering from error and continuing with last estimable coefficients.\n",e)),
		# finally=results<-results.k)
		# return(results)
	# } else if(is.null(penalties))
	# { # do penalties estimation
		# if(any(Rdebuglevel<0))cat("Finding number of principle components.\n")
		# if(k<1)stop("invalid number of principle components (k>=1).")
		# cvfunction<-function(penalties){
			# if(any(Rdebuglevel<0))cat("Evaluating penalties at:",penalties,"\n")
			# if(any(penalties<0))return(NA)
			# n2ll<-numeric(1)
			# inc=max(round(N*cv),1)
			# tryCatch(for(i in seq(1,N,by=inc)){
				# if(Rdebuglevel==-2)cat("Crossvalidation starting with: ",i,"\n")
				# exclude <- seq(from=i,by=1, length=inc)
				# ind <- !(obs %in% exclude)
				# fitdata<-mf[ind,]
				# testdata<-mf[!ind,]
				# if(Rdebuglevel==-1)assign(paste("fitdata",i,sep=""), fitdata, envir=globalenv())
				# model<-pfdaSingle(fitdata, knots=NULL, k=k, delta=delta, penalties=penalties, tol=tol, MaxIter=MaxIter, NCompTol=NCompTol, Basis=obase,...)
				# n2ll<-n2ll+logLik(model,newdata=fd(testdata[[2]],testdata[[1]],testdata[[3]]))
			# }, 
			# error=function(e){ n2ll<<-NA;if(any(Rdebuglevel<0))print(e) }, 
			# warning=function(w){ n2ll<<-NA;if(any(Rdebuglevel<0))print(w) })
			# if(exists("TrackCVFunction"))TrackCVFunction<<-rbind(TrackCVFunction,c(penalties,n2ll))
			# n2ll
		# }
		# optimmethod<-if(hasArg(optimMethod)) dots$optimmethod else "Nelder-Mead"
		# optimstart<-if(hasArg(optimstart)) dots$optimstart else c(1,1)
		# optimpar<-optim(optimstart,cvfunction, method=optimmethod)
		# if(any(Rdebuglevel<0))cat("Algorithm converged on:\n", optimpar$par, "\nafter ", optimpar$counts[1], " function calls reulsting in approximatly ", optimpar$counts*10, "total calls to pfdaSingle.\n")
		# returnvalue<-pfdaSingle(mf, knots=NULL, k=k, delta=delta, penalties=optimpar$par, tol=tol, MaxIter=MaxIter, NCompTol=NCompTol,Basis=obase,...)
		# returnvalue@CVlogLik<-optimpar$value
		# returnvalue
	# } else 
	# { # Given k and penalties
		# if(any(Rdebuglevel<0))cat("k is given",k,". penalties are given:",penalties,".\n")
		# if(any(Rdebuglevel==-3))browser()
		# if(k<1)stop("invalid number of principle components (k>=1).")
		# nobs<-as.integer(tapply(obs,obs,length))
		# M	<- sum(nobs)
		## if(M!=length(t))stop("Number of observations do not match (t, y, ons, nobs).")
		# B	<- evaluate(obase,t)
		# if(any(Rdebuglevel<0))cat("Computing ISD...\n")
		# ISD <- OuterProdSecondDerivative(obase)
		# p	<- as.integer(dim(obase)[2])
		# { #Compute Memory Requirements
			# k<-k
			# dpl_1 <- M + M*k + 2*k^2 
			# dpl_2 <- p^2 + M + M*k
			# dpl_3 <- M + p^2 + p + k^2
			# dpl_5 <- k^2 + k + 2*p^2 + p +  max( k*p , 8*p)
			# dpl_E <- M + M*k + N*k + 2* k^2 
			# dpl   <- N*p^2 + p + p*k + k +
			          # max(dpl_1, dpl_2, dpl_3, dpl_5, dpl_E ) #+ M + p
			# ipl<-8*p
		# }
		# if(any(Rdebuglevel<0)){
			# cat("garbage collection.\n")
			# gc()}
		# if(any(Rdebuglevel<0))cat("Passing computations to compiled code (.C)...\n")
		# results<-.C("pfdaSingle",
			# t			= as.double(t),
			# y			= as.double(y),
			# nobs		= as.integer(nobs),
			# M			= as.integer(M),
			# N			= as.integer(N),
			# k			= as.integer(k), 
			# B			= as.double(B), 
			# p			= as.integer(p),
			# delta		= as.double(delta),
			# lambda_mu	= as.double(penalties[1]),
			# lambda_f	= as.double(penalties[2]),
			# ISD 		= as.double(ISD),
			## /* State Values */
			# theta_mu	= double(p),
			# Theta_f		= double(p*k),
			# Da			= double(k),
			# seps		= double(1),
			# Alpha		= double(k*N),
			# Sigma_alpha_alpha = double(N*k**2),
			## /* Control Values */
			# tolerance	= as.double(tol),
			# Iterations	= as.integer(MaxIter),
			# debuglevel  = as.integer(Cdebuglevel),
			# double(dpl),integer(ipl)
		# )
		# if(any(Rdebuglevel<0))cat("returned from compiled code(",results$errorcode,").\n")
		# if(any(Rdebuglevel<0)){ cat("garbage collection.\n"); gc()}

		# dim(results$Theta_f)<-c(p,k)
		# dim(results$Alpha)<-c(N,k)
		# dim(results$Sigma_alpha_alpha)<-c(k,k,N)
		# if(results$Iterations==MaxIter) warning("Method did not converge. Results are returned for the last previous state.")
		# if(any(Rdebuglevel<0))cat("Making new FDAModelParameters object .\n")
		# parameters=new("FDAModelParameters",theta_mu=results$theta_mu, Theta_f=results$Theta_f,sigma = results$seps, Alpha=results$Alpha, Da = results$Da, npc=results$k, penalties=penalties, Sigma=results$Sigma_alpha_alpha, Ry=results$y)
		# if(any(Rdebuglevel<0))cat("Making new FDModel object.\n")
		# newFDModel = new("FDModel",Basis=obase, ConvergenceCriteria=results$tolerance, iterations = results$Iterations, FittingData = FunctionalData(t,y,obs),parameters=parameters)
		# if(any(Rdebuglevel<0))cat("exiting pfdaSingle(R).\n")
		# return(newFDModel)
	# }
# }

# setGeneric("fitModel",function(object,...)standardGeneric("fitModel"),where=globalenv())
# setMethod("fitModel","FunctionalData",function(object,...){
	# pfdaSingle(y~t,data=as(object,"data.frame"),obs=obs,...)
# })

# predict.FDModel<-function(object, ...){
	# dots<-list(...)
	# data<-as(if(hasArg(newdata)) dots$newdata else object@FittingData,"FunctionalData")
	# t<-as.numeric(data@t)
	# y<-as.numeric(data@y)
	# obs<-data@obs[drop=TRUE]
	
	# B<-evaluate(object@Basis, t)
	# BTf<-B%*%object@parameters@Theta_f
	# Btm<-B%*%object@parameters@theta_mu
	# Alpha<-matrix(nrow=nlevels(obs),ncol=object@parameters@npc)
	# yhat<-numeric(length(t))
	# for(i in nlevels(obs)){
		# ind <- i==as.integer(obs)
		# Alpha[i,]<-1/object@parameters@sigma*solve(diag(x=1/object@parameters@Da,nrow=object@parameters@npc)+1/object@parameters@sigma*crossprod(BTf[ind,]),crossprod(BTf[ind,],y[ind]-Btm[ind]))
		# yhat[ind]<-Btm[ind]+BTf[ind,,drop=FALSE]%*%Alpha[i,,drop=FALSE]
	# }
	# new("FDPredict",Alpha=Alpha,predicted=yhat,model=object, data = data)
# }

# setGeneric("predict",function(object,...)standardGeneric("predict"),where=globalenv())
# setMethod("predict", "FDModel", function(object,...)pfdaSinglePredict(obejct,...),where=globalenv())

# pfdaDual<-function(formula, data , obs,  knots=NULL, k=NULL, delta=1e-4, penalties=NULL, tol=1e-2, MaxIter=10000,NCompTol=1/12,cv=.1,...){
	 ## Hidden Named Arguments
	 ## useoptimstart tells whether to use the penalties from k step as starting points in penalty optimization
	 ## optimstart (c(1,1,1,1)) gives starting point for optim
	# dots<-list(...)
	# splineorder<-if(hasArg(splineorder)) dots$splineorder else 4
	# mf<-cl<-match.call()
    # m <- match(c("formula", "data", "subset", "obs", "na.action"), names(mf), nomatch=0L)
    # mf <- mf[c(1L, m)]
    # mf$drop.unused.levels <- TRUE
    # mf[[1L]] <- as.name("model.frame")
    # mf <- eval(mf, parent.frame())
	# if(length(attr(mf,"variables"))>3)stop("Only one predictor can be specified.")
	# obs<-as.factor(mf[[3]])[drop=TRUE]
	# oo<-order(obs)
	# mf<-mf[oo,];
	# mf[[3]]<-obs<-obs[oo]  # ensure that the data is sorted properly
	# mr<-model.response(mf)
	# if(NCOL(mr)>2)stop("Too many responses specified")
	# y<-mr[,1]
	# z<-mr[,2]
	# t<-mf[,2]
	# N<-nlevels(obs)
	# if(is.null(knots)){
		# knots<-quantile(t,seq(0,1,length=11))[c(rep(1,splineorder),seq(from=2,to=10),rep(11,splineorder))]
	# } 
	# else if(!all((t>=knots[splineorder]) && (t<=knots[length(knots)-splineorder+1])))stop("Knots do not sufficiently cover the domain.")
	# if(cv<0||cv>.5)stop("cv parameter must be between 0 and .5")
	# if(is.null(dots$Rdebuglevel))Rdebuglevel=0L else Rdebuglevel<-dots$Rdebuglevel
	# if(hasArg(Cdebug)) Cdebug=dots$Cdebug else if(is.null(dots$Cdebuglevel)) Cdebug=integer(0) else {
		# Cdebug<-make_c_debug(dots$Cdebuglevel)
	# }
	# if(is.null(k)){
		# if(any(Rdebuglevel<0))cat("Finding number of principle components.\n")
		# if(any(Rdebuglevel<0))cat("Finding single model for y.\n")
		# rs.y<-suppressWarnings(pfdaSingle(data.frame(y,t,obs),knots=knots,k=NULL,delta=delta,
			  # penalties=if(is.null(penalties)) NULL else penalties[c(1,3)], tol=tol, MaxIter=MaxIter, NCompTol=NCompTol,...))
		# if(any(Rdebuglevel<0))cat("Finding single model for z.\n")
		# rs.z<-suppressWarnings(pfdaSingle(data.frame(z,t,obs),knots=knots,k=NULL,delta=delta,
			  # penalties=if(is.null(penalties)) NULL else penalties[c(1,3)], tol=tol, MaxIter=MaxIter, NCompTol=NCompTol,...))
		# ka<-rs.y@parameters@npc
		# kb<-rs.z@parameters@npc
		# py<-rs.y@parameters@penalties
		# pz<-rs.z@parameters@penalties
		## StartFrom<-if(hasArg(StartFrom)) dots$StartFrom else if((!is.null(dots$useStartFrom) && dots$useStartFrom))new("PFDModel",Basis=rs.y@Basis,Analysis_y=rs.y@parameters,Analysis_z=rs.z@parameters,Lambda=solve(crossprod(rs.y@parameters@Alpha),crossprod(rs.y@parameters@Alpha,rs.z@parameters@Alpha))) else NULL
		# optimstart<-if(!hasArg(optimstart) && (is.null(dots$useoptimstart) || useoptimstart)) c(py[1],pz[1],py[2],pz[2]) else dots$optimstart
		# pfdaDual(mf, k=c(ka,kb), penalties=penalties, delta=delta,tol=tol, MaxIter=MaxIter,NCompTol=NCompTol,optimstart=optimstart,...)
	# } 
	# else if(is.null(penalties)) {
		# if(any(Rdebuglevel<0))print("Finding optimal penalties.")
		# cvfunction<-function(penalties){
			# if(any(Rdebuglevel<0))print(paste("Evaluating penalties at:",paste(penalties,collapse=" ")))
			# if(Rdebuglevel==-3)browser()
			# if(any(penalties<0))return(NA)
			# n2ll<-numeric(1)
			# inc=max(round(N*cv),1)
			# tryCatch(for(i in seq(1,N,by=inc)){
				# if(Rdebuglevel==-2)cat("Crossvalidation starting with: ",i,"\n")
				# exclude <- seq(from=i,by=1, length=inc)
				# ind <- !(obs %in% exclude)
				# fitdata<-mf[ind,]
				# testdata<-mf[!ind,]
				# if(Rdebuglevel==-1)assign(paste("fitdata",i,sep=""), fitdata, envir=globalenv())
				# model<-pfdaDual(fitdata, knots=knots, k=k, delta=delta, penalties=penalties, tol=tol, MaxIter=MaxIter, NCompTol=NCompTol,...)
				# n2ll<-n2ll+logLik(model,newdata=testdata)
			# }, 
			# error=function(e){ n2ll<-NA; if(any(Rdebuglevel<0))print(e) }, 
			# warning=function(w){ n2ll<-NA; if(any(Rdebuglevel<0))print(w) }
			# )
			# if(exists("pfdaDebug_TrackCVFunction"))pfdaDebug_TrackCVFunction<<-rbind(pfdaDebug_TrackCVFunction,c(penalties,n2ll))
			# n2ll
		# }
		# optimmethod<-if(hasArg(optimMethod)) dots$optimmethod else "Nelder-Mead"
		# optimstart<-if(hasArg(optimstart)) dots$optimstart else c(1,1,1,1)
		# optimpar<-optim(optimstart,cvfunction, method=optimmethod)
		# if(any(Rdebuglevel<0))cat("Algorithm converged on:\n", optimpar$par, "\nafter ", optimpar$counts[1], " function calls reulsting in approximatly", optimpar$counts*10, "total calls to pfdaSingle.\n")
		# returnvalue<-pfdaDual(mf, knots=knots, k=k, delta=delta, penalties=optimpar$par, tol=tol, MaxIter=MaxIter, NCompTol=NCompTol,...)
		# returnvalue@CVlogLik<-optimpar$value
		# returnvalue
	# }
	# else {
		# if(any(Rdebuglevel<0))cat("k is given",k,". penalties are given:",penalties,".\n")
		# if(Rdebuglevel==-3)browser()
		# if(any(k<1))stop("invalid number of principle components (k>=1).")
		# M	<- NROW(mf)
		# nobs <- tapply(obs,obs,length)
		# ka <- k_alpha <- k[1]
		# kb <- k_beta  <- k[2]
		# {
			# obase<-OrthogonalizeBasis(SplineBasis(knots))
			# p<-dim(obase)[2]
			# theta_mu			<-	double(p)
			# theta_nu			<-	double(p)
			# Theta_f 			<-	double(p*k_alpha)
			# Theta_g 			<-	double(p*k_beta)
			# D_alpha  			<-	double(k_alpha)
			# D_beta  			<-	double(k_beta)
			# Lambda  			<-	double(k_alpha*k_beta)
			# sigma_epsilon 		<-	double(1)
			# sigma_xi  			<-	double(1)
			# Alpha 				<- 	double(N*k_alpha)
			# Beta  				<-	double(N*k_beta)
			# Sigma_alpha_alpha 	<-	double((k_alpha**2)*N)
			# Sigma_alpha_beta  	<-	double(k_alpha*k_beta*N)
			# Sigma_beta_beta   	<-	double((k_beta**2)*N)
			# incInits<-FALSE
		# }
		# B	<- evaluate(obase,t)
		# ISD <- OuterProdSecondDerivative(obase)
		# { #Compute Memory Requirements
			# k<-max(ka,kb)
			# dpl_1     <- M + M*ka + 2*k^2 
			# dpl_2     <- p^2 + M + M*k
			# dpl_3     <- M + p^2 + p + k^2
			# dpl_4     <- k^2 + ka*kb
			# dpl_5_1   <- k + 2*p^2 + p +  
				# max(  outer_qf = k*p , eigens = 8*p)
			# dpl_5_2   <- ka^2 + kb*ka
			# dpl_5     <- ka^2 + kb^2 + 
				# max(dpl_5_1  ,dpl_5_2 )
			# dpl_E_1   <- kb^2 + ka^2 + 
				# max(outer_qf = ka*kb, 
					# sym_inv = 2* kb^2, 
					# inner_qf = ka*kb)
			# dpl_E_2_1 <- 2*k^2
			# dpl_E_2_2 <- 3*k^2
			# dpl_E_2   <- kb^2 + max(dpl_E_2_1,dpl_E_2_2,ka*kb)
			# dpl_E_3_1 <- 2*k*max(nobs)
			# dpl_E_3   <- dpl_E_3_1
			# dpl_E     <- 2*M + M*ka + M*kb + ka^2 + kb^2 + ka*kb + 
				# max(dpl_E_1, dpl_E_2, dpl_E_3)
			# dpl<- N*(p**2) + p*2 + p*ka+p*kb + ka*kb + ka + kb + 
				# max(dpl_1, dpl_2, dpl_3, dpl_4, dpl_5, dpl_E )
			# ipl<-8*p
		# }
		# if(any(Rdebuglevel<= -2)){
			# cat("dpl=",dpl,"\n")
			# cat("dpl_1=",dpl_1,"\n")
			# cat("dpl_2=",dpl_2,"\n")
			# cat("dpl_3=",dpl_3,"\n")
			# cat("dpl_4=",dpl_4,"\n")
			# cat("dpl_5=",dpl_5,"\n")
			# cat("dpl_E=",dpl_E,"\n")
		# }
		# results<-{ .C("pfdaDual",
			# t=as.double(t), y=as.double(y), z=as.double(z),
			# nobs=as.integer(nobs),
			# M=as.integer(M), 
			# N=N, 
			# k_alpha=as.integer(k_alpha),
			# k_beta=as.integer(k_beta),
			# B=as.double(B), 
			# p=as.integer(p),
			# delta=as.double(delta),
			# penalties = as.double(penalties),
			# ISD = as.double(OuterProdSecondDerivative(obase)),
			## end input values
			# theta_mu=as.double(theta_mu),
			# theta_nu=as.double(theta_nu),
			# Theta_f =as.double(Theta_f),
			# Theta_g =as.double(Theta_g),
			# D_alpha =as.double(D_alpha),
			# D_beta =as.double(D_beta),
			# Lambda =as.double(Lambda),
			# sigma_epsilon =	as.double(sigma_epsilon),
			# sigma_xi =as.double(sigma_xi),
			# Alpha =as.double(Alpha),
			# Beta =as.double(Beta),
			# Sigma_alpha_alpha =as.double(Sigma_alpha_alpha),
			# Sigma_alpha_beta  =as.double(Sigma_alpha_beta),
			# Sigma_beta_beta   =as.double(Sigma_beta_beta),
			## end Receiving values
			# tolerance=as.double(tol),
			# Iterations=as.integer(MaxIter),
			# incInits=as.integer(incInits),
			# debuglevel=as.integer(c(length(Cdebug),Cdebug)),
			# double(dpl),
			# integer(ipl)
		# ) }
		# if(any(Rdebuglevel<0))cat("returned from compiled code(",results$errorcode,").\n")
		# if(any(Rdebuglevel<0)){ cat("garbage collection.\n");gc() }
		# dim(results$Theta_f)<-c(results$p,k_alpha)
		# dim(results$Theta_g)<-c(results$p,k_beta)
		# dim(results$Alpha)<-c(results$N,k_alpha)
		# dim(results$Beta)<-c(results$N,k_beta)
		# dim(results$Lambda)<-c(k_beta,k_alpha)
		# dim(results$Alpha)<-c(results$N,k_alpha)
		# dim(results$Beta)<-c(results$N,k_beta)
		# dim(results$Sigma_alpha_alpha)<-c(k_alpha,k_alpha,results$N)
		# dim(results$Sigma_alpha_beta)<-c(k_alpha,k_beta,results$N)
		# dim(results$Sigma_beta_beta)<-c(k_beta,k_beta,results$N)
		# colnames(results$Lambda)<-paste("Alpha",seq_len(k_alpha))
		# rownames(results$Lambda)<-paste("Beta", seq_len(k_beta))
		# if(results$Iterations==MaxIter) warning("Method did not converge. Results are returned for the last previous state.")
		# new("PFDModel",Basis=obase, ConvergenceCriteria=results$tolerance, iterations = results$Iterations, 
			# FittingData = PairedFunctionalData(t=t,y=y,z=z,obs=obs),
			# Analysis_y=new("FDAModelParameters",theta_mu=results$theta_mu, Theta_f=results$Theta_f,sigma = results$sigma_epsilon, 
							# Alpha=results$Alpha, Da = results$D_alpha, npc=results$k_alpha, penalties=penalties[c(1,3)], 
							# Sigma=results$Sigma_alpha_alpha, Ry=results$y),
			# Analysis_z=new("FDAModelParameters",theta_mu=results$theta_nu, Theta_f=results$Theta_g,sigma = results$sigma_xi, 
							# Alpha=results$Beta, Da = results$D_beta, npc=results$k_beta, penalties=penalties[c(2,4)], 
							# Sigma=results$Sigma_beta_beta, Ry=results$z),
			# Sigma_yz = results$Sigma_alpha_beta, 
			# Lambda=results$Lambda
		# )
	# }
# }

# setMethod("fitModel","PFDModel",function(object,...){pfdaDual(cbind(y,z)~t,data=as(object,"data.frame"),obs=obs,...)})
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


#########################################################################################
#'  	Control values for the pfda function
#' 	
#' 	this is just a conveinience function to show the possible control values and gather them in one place.
#' 	
#' 	@param minimum.variance the minimun variance allowed in the algorhtim.  Used to prevent singularity problems.
#' 	@param convergence.tolerance  the tolerance to determine the convergence of the algorithm.
#' 	@param max.iterations the maximum number of iterations allowed in the algorithm
#' 	@param pc.tolerance the relative ratio of variance for adding a principle component to the model
#' 	@param nfolds the number of folds to use in cross validataion
#' 	@param expected.values.tolerance tolerance for determining convergence of estimating expected values in binary W step
#' 	@param turn on tracing messages
#' 	
#' 	
pfdaControl<-function(...,minimum.variance = 1e-4, convergence.tolerance = 1e-2, 
	max.iterations=10000, pc.tolerance=1/12, nfolds = 10, trace=TRUE,
	useC=TRUE,C.debug=Cdebug()){
	list(
		minimum.variance=as.double(minimum.variance),
		convergence.tolerance = as.double(convergence.tolerance),
		max.iterations=as.integer(max.iterations),
		pc.tolerance=as.double(pc.tolerance),
		nfolds =nfolds,
		# expected.value.tolerance =expected.value.tolerance,
		trace=trace,useC=useC,C.debug=as.integer(C.debug),
		...
	)
}


#########################################################################################
#'	Handles the dual principle componet for single response model (aka Calcium model)
#'	
#'	
#'	
#'	
#'	
#'	
#'	
#'	



#########################################################################################
#'	General Model Fitting For Principle Component Functional Data Analsysis
#'	
#'	@param formula
#'	@param data
#'	@param knots
#'	@param k
#'	@param penalties
#'	@param response.class
#'	@param control
#'	
#'	
#'	
#'	
#'	
#'	
# pfda<-function(formula, data=environment(formula), knots=NULL, k=NULL, penalties=NULL, response.class=NULL, control=pfdaControl(),...)
# { 	
	# if(class(formula)=='formula')          model <- as(pfdaParseFormula(formula,data=data),'FunctionalData')
	# else if(is(formula, 'FunctionalData')) model <- formula
	# else if(is(formula, 'data.frame'))     model <- as(formula,'FunctionalData')
	# else stop(gettextf('Unusable class (%s) passed for formula parameter.  Try specifying a formula.',class(formula)))
	# validObject(model)
	
	# response	<- model[[1]]
	# domain		<- model[[2]]
	# subject		<- model[[3]]
	# if(!(is.factor(subject)||is.integer(subject)))stop("bad class for subject")
	
	# NR<-NCOL(response)
	# if(NR>2)stop("too many responses found.")
	# if(NR==1&&control$trace)message("found single response")
	# if(NR==2&&control$trace)message("found two responses")
	
	# if(is.null(k)||any(is.na(k))){
	
	# } 
	# else if(is.null(penalties)||any(is.na(penalties))) {
		# if(is.null(penalties))penalties<-rep(NA,2*NR)
		## ' @TODO program fitting for specific penalties given others.
	# }
	# else {
		# if(length(penalties)!=NR*2){ if(!is.matrix(penalties))penalties<-matrix(penalties,NR,2)
			# message('penalties have been duplicated for both responses.')}
		# if(length(k)!=NR){ k<-rep(k,length.out=NR)
			# message("k has been replicated for responses.")}
		# if(is.null(response.class)){
			# pfdaInferClass<-function(x){
				# c<-class(x)
				# if(c=='numeric') 'linear' else
				# if(c=='logical') 'binary' else
				# if(c=='factor') 'binary' else
				# if(c=='integer'){
					# if(length(unique(x))==2) 'binary' else 'linear'
				# }
			# }
			# response.class=character(NR)
			# for(i in seq(NR)){
				# response.class[i]<-pfdaInferClass(response[[i]])
			# }
			# message(paste("infering response classes from data classes.",
					# paste("\tresponse class for ",names(response),"=",response.class,collapse="\n"),sep="\n"))
		# } else if(length(response.class)!=NR) { 
			# response.class=rep(response.class,length.out=NR)
			# message("response.class replicated for all responses.")
		# }
		
		# { #compute basis and related quantities
			# order<-if(hasArg(order))order else 4
			# if(is.null(knots))knots<-expand.knots(unique(quantile(domain,seq(0,1,length=11))),order=order) #' knots are determined by quantiles of the domain using linear interpolation
			# obase<-OBasis(knots,order=order)
			# Bmatrix<-evaluate(obase,domain)
			# Kmatrix <- OuterProdSecondDerivative(obase)
			# p	<- as.integer(dim(obase)[2])
		# }
		# { Cdebug<-if(is.null(control$Cdebug)) 0L else control$Cdebug
		  # if(length(Cdebug)!=Cdebug[1]&&all(Cdebug==0L))Cdebug<-c(length(Cdebug),Cdebug)}
		# M<- NROW(response)
		# if(NR==1){
			# nobs<-as.integer(tapply(subject,subject,length))
			# N<-length(nobs)
			# if(response.class=='linear'){
				# { #Compute Memory Requirements
					# k<-k
					# dpl_1 <- M + M*k + 2*k^2 
					# dpl_2 <- p^2 + M + M*k
					# dpl_3 <- M + p^2 + p + k^2
					# dpl_5 <- k^2 + k + 2*p^2 + p +  max( k*p , 8*p)
					# dpl_E <- M + M*k + N*k + 2* k^2 
					# dpl   <- N*p^2 + p + p*k + k +
							# max(dpl_1, dpl_2, dpl_3, dpl_5, dpl_E ) #+ M + p
					# ipl<-8*p
				# }
				# results <- { .C("pfdaSingle",
					# t              = as.double(domain),
					# y              = as.double(response[[1]]),
					# nobs           = as.integer(nobs),
					# M              = as.integer(M),
					# N              = as.integer(N),
					# k              = as.integer(k), 
					# B              = as.double(Bmatrix), 
					# p              = as.integer(p),
					# delta          = as.double(control$minimum.variance),
					# lambda_mu      = as.double(penalties[1]),
					# lambda_f       = as.double(penalties[2]),
					# Kmatrix        = as.double(Kmatrix),
					## /* State Values */
					# theta_mu       = double(p),
					# Theta_f        = double(p*k),
					# Da             = double(k),
					# sigma          = double(1),
					# Alpha          = double(k*N),
					# Sigma          = double(N*k**2),
					## /* Control Values */
					# tolerance	     = as.double(control$convergence.tolerance),
					# iterations	= as.integer(control$max.iterations),
					# debuglevel     = as.integer(control$C.debug),
					# double(dpl),integer(ipl)
				# )}	
				# { #post C processing
					# dim(results$Theta_f)<-c(p,k)
					# dim(results$Alpha)<-c(N,k)
					# dim(results$Sigma)<-c(k,k,N)
					# if(results$iterations==control$max.iterations) 
						# warning("Method did not converge. Results are returned for the last previous state.")
					# parameters=new("FDAModelParameters",theta_mu=results$theta_mu, Theta_f=results$Theta_f,sigma = results$seps, Alpha=results$Alpha, Da = results$Da, npc=results$k, penalties=penalties, Sigma=results$Sigma_alpha_alpha, Ry=results$y)
					# newFDModel = new("FDModel",Basis=obase, ConvergenceCriteria=results$tolerance, iterations = results$Iterations, FittingData = as('FunctionalData'),parameters=parameters)
					# return(newFDModel)
				# }
			# } 
			# else if(response.class=='binary'){
				# { #compute memory
					# ni=max(nobs)
					# kr=100
					# dpl = M + sum(nobs^2) + N*k^2 + N*p^2 + p + p*k + k + max(
						# ni*kr + kr + (ni * k + ni + k^2 + k + p + p^2 + k*p + 3*k),# step W
						# M + M*k + N*k + 2* k^2 + 2 * k + k *ni        ,# step 1/E
						# p^2+ M+ M*k                                   ,# step 2  
						# M + p^2 + p + k^2                             ,# step 3  
						# k + 2*k^2 + 2*p^2 + p + p*max(k,8)            ,# step 4  
						# 2*p^2 + M + p*N + 8*p)# inits
					# ipl = 8*p
				# }
				# results <- { .C("pfda_bin_single",
					# y              = as.integer(response[[1]]),
					# nobs           = as.integer(nobs),
					# M              = as.integer(M),
					# N              = as.integer(N),
					# k              = as.integer(k), 
					# B              = as.double(Bmatrix), 
					# p              = as.integer(p),
					# lambda_mu      = as.double(penalties[1]),
					# lambda_f       = as.double(penalties[2]),
					# Kmatrix        = as.double(Kmatrix),
					## /* State Values */
					# theta_mu	     = double(p),
					# Theta_f        = double(p*k),
					# Da             = double(k),
					# Alpha          = double(k*N),
					# Sigma          = double(N*k**2),
					## /* Control Values */
					# minimun.variance = as.double(control$minimum.variance),
					# c.tol            = as.double(control$convergence.tolerance),
					# iterations	     = as.integer(control$max.iterations),
					# debuglevel       = as.integer(Cdebug),
					# double(dpl),integer(ipl)
				# )}
				# { #post C processing
					# if(control$trace){
						# cat('Testing garbage collection\n')
						# gc()
					# }
					# if(!is.null(control$debug.capture.w)){
						# if(control$trace)cat('capturing the current values of w and ww variables\n(this breaks several laws of programming, \ndo not attempt to do this unless you really know what you are doing.\n')
						# .pfda.debug.w<<-results[[22]][1:M]
						
						# .ww.raw<-results[[22]][M+1:sum(nobs^2)]
						# .ww<-vector('list',N)
						# for(i in 1:N){
							# .ww[[i]] <- matrix(.ww.raw[seq(nobs[i]^2)],nobs[i],nobs[i])
							# .ww[[i]][lower.tri(.ww[[i]])]<-t(.ww[[i]])[lower.tri(.ww[[i]])]
							# .ww.raw<-.ww.raw[-seq(nobs[i]^2)]
						# }
						# .pfda.debug.ww<<-.ww
						
						# .pfda.debug.aa<<-results[[22]][M+sum(nobs^2)+1:(N*k^2)]
						# dim(.pfda.debug.aa)<-c(k,k,N)
					# }
					# Theta_f<-matrix(results$Theta_f,p,k)
					# Alpha<-matrix(results$Alpha,N,k)
					# Sigma<-array(results$Sigma,dim=c(k,k,N))
					# if(results$iterations==control$max.iterations) 
						# warning("Method did not converge. Results are returned for the last previous state.")
					# parameters={ new("pfdaBinaryModelParameters",
						# theta_mu=results$theta_mu, Theta_f=Theta_f, 
						# Alpha=Alpha, Da = results$Da, npc=results$k, 
						# penalties=penalties, Sigma=Sigma)}
					# newFDModel = new("pfdaBinaryModel",Basis=obase, ConvergenceCriteria=results$c.tol, iterations = results$iterations, FittingData = as(model,'FunctionalData'),parameters=parameters)
					# return(newFDModel)
				# }
			# } 
			# else stop("how did you get here? are you trying to break my code?")
		# } 
		# else if(NR==2) {
			# y<-response[[1]]
			# z<-response[[2]]
		

		# } 
		# else stop("there must be a bug because you should not have reached here :(")
	# }
# }
