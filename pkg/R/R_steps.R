{ # global scope variables
# these variables are defined for the sole purpose of confusing the checker
# Each of these should be masked or ignored by a different frame path.
# This is necessary due to the eval/structure programming method employed here.
subject<-
y<-
Z<-
B<-
Bt<-
Bx<-
fname<-
s.xi<-
s.eps<-
ka<-
kb<-
kz<-
nobs<-
M<-
N<-
p<-
K<-
Kt<-
Kx<-
ln<-
lf<-
lg<-
min.v<-
max.I<-
tol<-
name.t<-
name.x<-
tbase<-
xbase<-
NULL
}
{ # general steps
.gen.tm<-function(y,B,subject,tf,alpha,sigma,l,K){
	left<-sigma*l*K
	right<-matrix(0,NCOL(B),1)
	for(i in seq(nlevels(subject))){
		ix <- i == as.integer(subject)
		Bi<-B[ix,]
		left  = left + crossprod(Bi)
		right = right + crossprod(Bi,y[ix]-Bi%*%tf%*%alpha[i,])
	}
	solve(left,right)
}
.gen.tf<-function(y,B,subject,tm,tf,alpha,sigma,aa,lf,K){
	n<-nlevels(subject)
	k<-NCOL(tf)
	for(est in seq(k)){
		left  <- sigma*lf*K
		right <- matrix(0,NCOL(B),1)
		for(i in seq(n)){
			ix <- i == as.integer(subject)
			Bi <- B[ix,]
 			left  = left  + crossprod(Bi)*(aa[est,est,i])
			sne <- matrix(0,NROW(Bi),1)
			for(l in setdiff(seq(k),est))sne = sne + Bi%*%tf[,l]*aa[est,l,i]
			right = right + crossprod(Bi,	(y[ix]-Bi%*%tm)*alpha[i,est] - sne	)
		}
		tf[,est]<-solve(left,right)
	}
	tf*rep(sign(tf[1,]),each=NROW(tf))
}
.gen.dual.lambda<-function(sum.aa,sum.ab){
	t(sum.ab)%*%solve(sum.aa)
}
.gen.orthog<-function(tf,alpha,sum.aa){
	n<-NROW(alpha)
	k<-NCOL(tf)
	Sa<-sum.aa/n
	e<-eigen(tf%*%tcrossprod(Sa,tf))
	qf = e$vectors[,seq_len(k),drop=FALSE]
	qf = qf*rep(sign(qf[1,]),each=NROW(qf))
	Da = x=e$values[seq_len(k)]
	R = crossprod(qf,tf)
	list(tf=qf,D=Da,d=tcrossprod(alpha,R),transformation=R)
}
.gen.dual.sigmas<-function(By,Bz,tf,tg,lambda,Da,Db,s.eps,s.xi){
	Seta <- diag(Db,length(Db))-lambda%*%tcrossprod(diag(Da,length(Da)),lambda)
	Setai<- solve(Seta)
	phi  <- By %*% tf
	psi  <- Bz %*% tg
	Sa   <- diag(1/Da,length(Da)) + crossprod(lambda,Setai)%*%lambda + crossprod(phi)/s.eps
	Sab1 <- -crossprod(lambda,Setai)
	Sb   <- Setai+crossprod(psi)/s.xi
	Saa  <- solve(Sa - Sab1%*%solve(Sb,t(Sab1)))
	Sab  <- -Saa%*%Sab1%*%solve(Sb)
	Sbb  <- solve(Sb - crossprod(Sab1,solve(Sa,Sab1)))
	list(Saa=Saa,Sab=Sab,Sbb=Sbb)
}
.gen.symblock.solve<-function(A,B,C){
#!	Solves matrixof form
#!  [[  A   B]
#!	 [t(B) C]]
#!	Returning the block inverses
#!

	stopifnot(NROW(A)==NROW(B))
	stopifnot(NCOL(C)==NCOL(B))

	svAB = solve(A,B)
	Ci = solve(C-crossprod(B,svAB))
	Bi = -svAB%*%Ci
	Ai = solve(A)+ svAB%*%tcrossprod(Ci,svAB)

	list(Ai,Bi,Ci)
}
.X.handle.z<-expression({ # handle Z
	stopifnot(exists('Z',inherits=FALSE))
	if(missing(Z)||is.null(Z))
		Z = matrix(nrow=length(y),ncol=0)
	else {
		if(is(Z,"formula"))
			Z = model.matrix(Z)
	  else Z = as.matrix(Z)
	}
})
.X.subset<-expression({
	stopifnot(exists('subset',inherits=FALSE))
	if(!missing(subset) && !is.null(subset)){
		if(exists('Z',inherits=FALSE))Z<-subset(Z,subset)
		y<-subset(y,subset)
		if(exists('z',inherits=FALSE))z<-subset(z,subset)
		t<-subset(t,subset)
		if(exists('x',inherits=FALSE))x<-subset(x,subset)
		subject<-subject[subject,drop=T]
	}
})
.X.single.knots<-expression({ # knots identification
	if(!exists("tbase",inherits=FALSE) || is.null(tbase)){
		if(is.null(knots)){
			kt<-expand.knots(unique(quantile(t,seq(0,1,length.out=control$nknots))))
		} else kt<-knots
		tbase = OBasis(kt)
	}
	Bt = evaluate(tbase,t)
	Kt = OuterProdSecondDerivative(tbase)
})
.X.single.penalties<-expression({ # penalties
	if(is.null(penalties)){
		penalties<- if(is.null(df)) rep(NA,2) else c(l.from.df(df[1],Bt,Kt),l.from.df(df[2],Bt,Kt))
	}
})
.X.single.optimize.npc<-expression({
	# funcall$k=1
	# if(exists('y'))funcall$y=y
	# if(exists('Z'))funcall$y=Z
	# if(exists('t'))funcall$y=t
	# if(exists('subject'))funcall$y=subject
	# model.k0<-eval(funcall,env=attr(funcall,'envir'))
	message("optimizing number of principal components using AIC")
	k=1
	model.k0 <- RecallWith(k=k)
	while(TRUE){
		# funcall$k<-funcall$k+1
		# model.k1<-try(eval(funcall,env=attr(funcall,'envir')),silent=TRUE)
		k=k+1
		model.k1 <- RecallWith(k=k)
		if(class(model.k1)[1]=="try-error")return(model.k0)
		else if(AIC(model.k0)<AIC(model.k1))return(model.k0)
		else model.k0<-model.k1
	}
})
.F.single.optimize.npc<-function(){
	message("optimizing number of principal components using AIC")
	# localfuncs("RecallWith")
	k=1
	model.k0 <- tryCatch(RecallWith(k=k,fname=fname),warning=function(w)stop(w$message))
	tryCatch(while(TRUE){ 
		k=k+1
		model.k1 <- RecallWith(k=k,fname=fname)
		if(AIC(model.k0)<AIC(model.k1))return(model.k0)
		else model.k0<-model.k1
	}, warning=function(w){
		if(w$message=="EM-algorithm did not converge")model.k0 else stop(w$message)
	}, error=function(e)model.k0)
}
.X.dual.k<-expression({ # resloves k input for dual pc (cc/bc/add)
	if(is.null(k)) k<-rep(NA,2)
	if(length(k)!=2) k<-rep(as.vector(k),length.out=2)
	ka = as.integer(k[1])
	kb = as.integer(k[2])
	# if(names(k)==NULL) names(k)<-c(name.t,name.x)
})
.X.dual.penalties <-expression({ #  resolve penalties input
	if(is.null(penalties)){
		penalties<- if(exists('Bx',inherits=FALSE))
			if(is.null(df)) matrix(NA,2,2) else c(l.from.df(df[1],Bt,Kt),l.from.df(df[2],Bx,Kx),l.from.df(df[3],Bt,Kt),l.from.df(df[4],Bx,Kx))
	  else
			if(is.null(df)) matrix(NA,2,2) else c(l.from.df(df[1],Bt,Kt),l.from.df(df[2],Bt,Kt),l.from.df(df[3],Bt,Kt),l.from.df(df[4],Bt,Kt))
	}
	if(!class(penalties)=='matrix') penalties<-matrix(penalties,2,2)
	# if(is.null(colnames(penalties))) colnames(penalties)<-c(name.t,name.x)
})
.X.read.parameters <- expression({
		if(exists('Z',inherits=FALSE)){ kz  = if(is.null(Z))0L else NCOL(Z)}
		k  = as.integer(k)
		N   = nlevels(subject)
		nobs= table(subject)
		M   = length(y)
		p  = ncol(Bt)
})
.X.binary.y<-expression({ # enforces the conditions on a binary y
	y<-if(is(y,'factor')){
		stopifnot(nlevels(y)==2)
		as.integer(y)-1
	} else if(is(y,'logical')){
		as.integer(y)
	} else {
		stopifnot(identical(as.integer(sort(unique(y))), as.integer(0:1)))
		as.integer(y)
	}
})
.X.funcall.replacevars<-expression({
	funcall$y=y
	if(exists("z",inherits = FALSE))funcall$z=z
	if(exists("Z",inherits = FALSE))funcall$Z=Z
	funcall$t = t
	if(exists("x",inherits = FALSE))funcall$x = x
	funcall$subject = subject
	funcall$k=k
	funcall$penalties = penalties
	funcall$df = df
	if(exists("bases",inherits=FALSE))funcall$bases=bases
	funcall$knots = knots
	funcall$control = control
})
.X.optimize.penalties<-expression({
	eval(.X.funcall.replacevars)
	pix<-which(is.na(penalties))
	if(control$penalty.method=='CV'){
		message("optimizing penalties using cross-validation")
		if(is.null(control$folds))control$folds<-cv.folds(nlevels(subject),10)
		ews<-function(s,p){
			ix = !(subject %in% s)
			fc<-funcall
			fc$penalties=p
			fc$subset=ix
			model<-eval(fc,env=attr(fc,'envir'))
			logLik(model,newdata=list(y=y[!ix],Z=Z[!ix,,drop=FALSE],Bt=Bt[!ix,,drop=FALSE],Bx=Bx[!ix,,drop=FALSE],subject=subject[!ix,drop=TRUE]),n2L=TRUE)
		}
		cvf<-function(pen,...){
			p<-penalties
			p[pix]<-exp(pen)
			if(any(is.infinite(p)))return(Inf)
			cvl<-try(lapply(control$folds,ews,p=p),TRUE)
			if(class(cvl)=="try-error")return(NA)
			(sum(unlist(cvl)))
		}
		if(is.null(control$optim.method))control$optim.method<-"Nelder-Mead"
		if(is.null(control$optim.start))control$optim.start<-rep(1,length(pix))
		optimpar<-optim(log(control$optim.start),cvf,method=control$optim.method)
		penalties[pix]<-exp(optimpar$par)
		funcall$penalties=penalties
		eval(funcall,env=attr(funcall,'envir'))
	} else
	if(control$penalty.method=='AIC') {
		message("optimizing penalties using AIC")
		aicf<-function(pen){
			fc<-funcall
			p<-penalties
	  	p[pix]<-exp(pen)
			if(any(is.infinite(p)))return(Inf)
			fc$penalties<-p
			m<-try(eval(fc,env=attr(fc,'envir')),silent=TRUE)
			if(class(m)[1]=="try-error") NA else AIC(m)
		}
		if(is.null(control$optim.method))control$optim.method<-"Nelder-Mead"
		if(is.null(control$optim.start))control$optim.start<-rep(l.from.df(2.1,Bt,Kt),length(pix))
		optimpar<-optim(log(control$optim.start),aicf,method=control$optim.method)
		penalties[pix] <- exp(optimpar$par)
		funcall$penalties=penalties
		eval(funcall,env=attr(funcall,'envir'))
	}
})
.F.optimize.penalties<-function(){
#	eval(.X.funcall.replacevars)
	pix<-which(is.na(penalties))
	if(control$penalty.method=='CV'){
		message("optimizing penalties using cross-validation")
		if(is.null(control$folds))control$folds<-cv.folds(nlevels(subject),10)
		ews<-function(s,p){
			ix = !(subject %in% s)
			fc<-funcall
			fc$penalties=p
			fc$subset=ix
			model<-eval(fc,env=attr(fc,'envir'))
			logLik(model,newdata=list(y=y[!ix],Z=Z[!ix,,drop=FALSE],Bt=Bt[!ix,,drop=FALSE],Bx=Bx[!ix,,drop=FALSE],subject=subject[!ix,drop=TRUE]),n2L=TRUE)
		}
		cvf<-function(pen,...){
			p<-penalties
			p[pix]<-exp(pen)
			if(any(is.infinite(p)))return(Inf)
			cvl<-try(lapply(control$folds,ews,p=p),TRUE)
			if(class(cvl)=="try-error")return(NA)
			(sum(unlist(cvl)))
		}
		if(is.null(control$optim.method))control$optim.method<-"Nelder-Mead"
		if(is.null(control$optim.start))control$optim.start<-rep(1,length(pix))
		optimpar<-optim(log(control$optim.start),cvf,method=control$optim.method)
		penalties[pix]<-exp(optimpar$par)
		funcall$penalties=penalties
		eval(funcall,env=attr(funcall,'envir'))
	} else
	if(control$penalty.method=='AIC') {
		message("optimizing penalties using AIC")
		aicf<-function(pen){
			# fc<-funcall
			p<-penalties
	  	p[pix]<-exp(pen)
			if(any(is.infinite(p)))return(Inf)
			# fc$penalties<-p
			tryCatch({
				m<-RecallWith(penalties=p,fname=fname)
				if(class(m)[1]=="try-error") NA else AIC(m)
				}
				,error=function(e)NA
				,warning=function(w){
					if(w$message=="EM-algorithm did not converge") return(NA)
					else stop(w$message)
				}
			)
		}
		if(is.null(control$optim.method))control$optim.method<-"Nelder-Mead"
		if(is.null(control$optim.start))control$optim.start<-rep(l.from.df(2.1,Bt,Kt),length(pix))
		optimpar<-optim(log(control$optim.start),aicf,method=control$optim.method)
		penalties[pix] <- exp(optimpar$par)
		# funcall$penalties=penalties
		# eval(funcall,env=attr(funcall,'envir'))
		RecallWith(penalties=penalties,fname=fname)
	}
}
}
{ # utilities and convenience functions
.u.single.resid<-function(y,B,subject,tm,tf,alpha){
	for(i in seq_len(nlevels(subject))){
		ix = i == as.integer(subject)
		y[ix]<-y[ix]-B[ix,,drop=FALSE]%*%(tm+tf%*%alpha[i,])
	}
	y
}
.u.orthogonalize<-function(U,v){
	# U is assumed to be orthogonal already
	u = as.vector(if(ncol(U)>0) v-U%*%as.vector(crossprod(v,U)) else v)
	u/sqrt(crossprod(u))
}
cv.folds<-function(n, folds = 10){
    split(sample(1:n), rep(1:folds, length = n))
}
positive.first.row<-function(X){
	stopifnot(is.matrix(X))
	X*rep(sign(X[1,]),each=nrow(X))
}
.pfda.df<-Vectorize(function(B,l,K,sigma=1){
	tr<-function(x)sum(diag(x))
	while(TRUE){
		L <- crossprod(B)+sigma*l*K
		if(kappa(L)<1e14)break
		else if(l>1) l <- l*.90 else NA #l <- l*1.1
	}
	tr(solve(L,crossprod(B)))
},'l')
l.from.df<-function(df,B,K,sigma=1)if(is.na(df)) NA else uniroot(function(l).pfda.df(B,l,K,sigma)-df,c(1e-4,1e10))$root
RecallWith<-function(...,fname){
	calls <- unlist(lapply(sys.calls(),function(x)deparse(x[[1]])[1]))
	if(missing(fname)){
		iRW<-max(which(calls=="RecallWith"))
		if(iRW==1) stop("RecallWith must be called within another function")
		unacceptable.calls<-c(which(calls=="{"),unlist(sapply(c("^eval","^function","^try","TryCatch","do.in.envir","^debug:::","^get","assign","mlocal"),grep,calls)))
		acceptable.calls<- setdiff(seq_len(iRW-1),unacceptable.calls)
		w = max(acceptable.calls)
		fname<-calls[w]# deparse(match.call(call = sys.call(w))[[1]])
	} else {
		w = max(which(calls==fname))
	}
	oldArgs<-formals(sys.function(w))
	if(is.null(oldArgs))oldArgs<-list()
	for(a in names(oldArgs)){
		try(oldArgs[[a]]<-get(a,sys.frame(w)),silent=TRUE)
	}
	do.call(fname,as.list(modifyList(oldArgs,list(...))))
}
}
{ # single version
.single.c.i<-function(y,B,subject,k,min.v){
	tm = solve(crossprod(B)+diag(min.v,NCOL(B)),crossprod(B,y))
	Ry = y-B%*%tm
	N = nlevels(subject)
	tfa<-matrix(0,NCOL(B),N)
	T=matrix(0,NCOL(B),NCOL(B))
	for(i in seq_len(nlevels(subject))){
		ix <- i==as.integer(subject)
		tfa[,i] = solve(crossprod(B[ix,])+min.v*diag(NCOL(B)),crossprod(B[ix,],Ry[ix]))
		Ry[ix]<-Ry[ix]-B[ix,]%*%tfa[,i]
		T<-T+tcrossprod(tfa[,i])
	}
	e<-eigen(T)
	tf<-e$vectors[,seq_len(k),drop=FALSE]
	tf<-tf*rep(sign(tf[1,]),each=NROW(tf))
	alpha<-matrix(0,N,k)
	aa<-list()
	for(i in seq_len(nlevels(subject))){
		alpha[i,]<-solve(crossprod(tf),crossprod(tf,tfa[,i]))
		aa[[i]] <- tcrossprod(alpha[i,])
	}
	sigma=crossprod(.u.single.resid(y,B,subject,tm,tf,alpha))/length(y)
	return(list(tm=tm,tf=tf,alpha=alpha,aa=aa,Da=e$values[seq_len(k)],sigma=sigma))
}
.single.c.1.1<-function(yi,Bi,tm,tf,Da,sigma){
	Sa = solve(crossprod(Bi%*%tf)+diag(1/Da,length(Da)))
	mu = Sa%*%crossprod(Bi%*%tf,yi-Bi%*%tm)/sigma
	aa = tcrossprod(mu)+Sa
	list(alpha=mu,Sa=Sa,aa=aa)
}
.single.c.1<-function(y,B,subject,tm,tf,Da,sigma){
	alpha<-matrix(nrow=nlevels(subject),ncol=NCOL(tf))
	aa<-array(0,dim=c(NCOL(tf),NCOL(tf),nlevels(subject)))
	Sa<-array(0,dim=c(NCOL(tf),NCOL(tf),nlevels(subject)))
	for(i in seq_len(nlevels(subject))){
		ix = i==as.integer(subject)
		r<-.single.c.1.1(y[ix],B[ix,],tm,tf,Da,sigma)
		alpha[i,]<-r$alpha
		Sa[,,i]<-r$Sa
		aa[,,i]<-r$aa
	}
	list(alpha=alpha,Saa=Sa,aa=aa)
}
.single.c.2<-function(y,B,subject,tm,tf,alpha,Sa){
	Ry <- y-B%*%tm
	s<-0
	for(i in seq_len(nlevels(subject))){
		ix<- i==as.integer(subject)
		phi=B[ix,,drop=FALSE]%*%tf
		s<-s+crossprod(Ry[ix]-phi%*%alpha[i,])+sum(diag(phi%*%Sa[,,i]%*%t(phi)))
	}
	list(sigma=as.vector(s/length(subject)))
}
.single.c.3<-function(y,B,subject,tf,alpha,sigma,lm,K){
	list(tm=.gen.tm(y,B,subject,tf,alpha,sigma,lm,K))
}
.single.c.4<-function(y,B,subject,tm,tf,alpha,sigma,aa,lf,K){
	list(tf=.gen.tf(y,B,subject,tf,tm,alpha,sigma,aa,lf,K))
}
.single.c.5<-function(tf,alpha,aa){
	saa<-matrix(0,NCOL(tf),NCOL(tf))
	for(sub in seq(NROW(alpha))) saa = saa + aa[[sub]]
	a<-.gen.orthog(tf,alpha,saa)
	list(tf=a$tf,Da=a$D, alpha = a$d)
}
single.c.core<-function(y,B,subject,k,lm,lf,K,min.v,max.I,tol){
	{ #Initial values
		r<-.single.c.i(y,B,subject,k,min.v)
		tm=r$tm
		tf=r$tf
		alpha=r$alpha
		aa=r$aa
		Da=r$Da
		sigma=r$sigma
		I=0
		.cc<-numeric(0)
	}
	while(I<-I+1){
		s1<-.single.c.1(y,B,subject,tm,tf,Da,sigma)
		s2<-.single.c.2(y,B,subject,tm,tf,s1$alpha,s1$Saa)
		s3<-.single.c.3(y,B,subject,tf,s1$alpha,s2$sigma,lm,K)
		s4<-.single.c.4(y,B,subject,s3$tm,tf,s1$alpha,s2$sigma,s1$aa,lf,K)
		s5<-.single.c.5(s4$tf,s1$alpha,s1$aa)
		{ # convergence checks
			ccl<-numeric(0)
			ccl['tm'] <- sum(abs((tm-s3$tm)/s3$tm))
			ccl['tf'] <- sum(abs((tf-s5$tf)/s5$tf))
			ccl['Da'] <- sum(abs((Da-s5$Da)/s5$Da))
			cc <- sum(ccl)
			.cc<-c(.cc,cc)
			if(cc<tol)break
			if(I>=max.I){
				warning('Maximum number of iterations exceeded, convergence not obtained.')
				break
			}
			{ # reassign values to
				tm<-s3$tm
				tf<-s5$tf
				alpha<-s5$alpha
				Da<-s5$Da
				sigma<-s2$sigma
				aa<-s5$aa
				Saa<-s1$Saa
			}
		}
	}
	list(y=y,Bt=B,subject=subject,tm=tm,tf=tf,alpha=alpha, Da=Da,sigma=sigma,aa=aa,Saa=Saa,lm=lm,lf=lf,K=K,k=k,I=I,cc=cc)
}
.single.c.n2L<-function(y, subject, B, tm, tf, Da, sigma){
	phi<-B%*%tf
	L <- 0.0
	for(i in seq_len(nlevels(subject))){
		ix<-i==as(subject,"integer")
		Sy<-phi[ix,,drop=FALSE]%*%tcrossprod(diag(x=Da,ncol=length(Da)),phi[ix,,drop=FALSE])+diag(sigma,nrow=sum(ix))
		Ry<-y[ix]-B[ix,]%*%tm
		D<-determinant(Sy,logarithm=TRUE)
		if(D$sign<0)stop("Negative determinant encountered for variance matrix.")
		L<-L+as.numeric(D$modulus+crossprod(Ry,solve(Sy,Ry)))
	}
	L
}
AIC.pfda.single.c<-function(object,...){
	with(object,{
		as.vector(.single.c.n2L(y, subject, Bt, tm, tf, Da, sigma))+
		2*(.pfda.df(Bt,lm,K,sigma)+k*.pfda.df(Bt,lf,K,sigma))
	})
}
# AIC.pfda.single.c.R<-function(object,...)
logLik.pfda.single.c.R<-logLik.pfda.single.c.rawC<-function(object,...,newdata=NULL, n2L=TRUE){
	r<-with(object,with(newdata,.single.c.n2L(y,subject,B,tm,tf,Da,sigma)))
	if(n2L) r else exp(-r/2)
}
single.c<-function(y,Z,t,subject,knots=NULL,penalties=NULL,df=NULL,k=NULL,control=pfdaControl(),subset=NULL){
	{ # setup
	fname = deparse(match.call()[[1L]])
	localfuncs(c('.F.single.optimize.npc','.F.optimize.penalties'))
	name.t = deparse(substitute(t))
	eval(.X.handle.z)
	eval(.X.subset)
	eval(.X.single.knots)
	eval(.X.single.penalties)
	}
	if(is.null(k)||any(is.na(k))){
		# funcall <- match.call()
		#stop("identification of number of principle components is not done yet.")
		# eval(.X.single.optimize.npc)
		.F.single.optimize.npc()
	} else
	if (any(is.na(penalties))) {
		funcall <- match.call()
		.F.optimize.penalties()
	}
	else {
		rtn<-if(control$useC){
			{ # setup for passing to Compiled code
				eval(.X.read.parameters)
				{ #Compute Memory Requirements
					pfda_computeResid =           M*k
					pfda_s_i =                    p*p + M + N*p + p*p + max(
					                                pfda_eigens =8*p,
					                                pfda_computeResid )
					single_c_resid =              pfda_computeResid
					pfda_m1 =                     M + M*k + 2*k^2
					pfda_m2 =                     M + M*k + p^2
					pfda_m3_for_subjectnum =      p
					pfda_m3_for_estcol =          p^2 + pfda_m3_for_subjectnum
					pfda_m3_core =                M + pfda_m3_for_estcol
					pfda_m3 =                     k^2 + pfda_m3_core
					pfda_m5_0 =                   p^2 + p^2 + p + k + max(
																					pfda_eigens = 8*p,
																					pfda_matrix_outer_quadratic_form =k*p)
					pfda_m5_1 =                   k^2 + pfda_m5_0
					pfdaSingle_m5 =               k^2 + pfda_m5_1
					pfdaSingle_e_1 =              k + 2*k^2
					pfdaSingle_e =                M + M*k + k*N + max(
																					pfda_computeResid ,
																					pfdaSingle_e_1 )
					single_c_E =                  M + pfdaSingle_e
					single_c_unpenalized =        M + kz*kz + max(
																					single_c_resid ,
																					dsysv_=10*kz)
					single_c_penalized =          M + pfda_m2
					single_c_princcomp =          M + p*p*N + pfda_m3_core
					single_c_variances =          M + max(
																					pfda_m1 ,
																					pfdaSingle_m5 )
					single_c_core =               p*p*N + pfda_s_i + p + p*k + k + max(
																					single_c_E ,
																					single_c_unpenalized ,
																					single_c_penalized ,
																					single_c_princcomp ,
																					single_c_variances ,
																					single_c_resid )
					ipl<-8*p
					dpl<- single_c_core
				}
			}
			structure(.C('single_c_core', residuals=y, Z=Z, Bt=Bt, tz=double(kz), tm=double(p), tf=matrix(0,p,k), alpha=matrix(0,N,k), Da=double(k), sigma=0, aa=array(0,dim=c(k,k,N)), Saa=array(0,dim=c(k,k,N)), nobs=nobs, N=N, M=M, kz=kz, k=k, p=p, lm=penalties[1], lf=penalties[2], K=Kt, minV=control$minimum.variance, Iterations=control$max.iterations, tol=control$convergence.tolerance, dl=control$C.debug, dp=double(dpl), ip=integer(max(6*p,kz)))
				,class=c('pfda.single.c.rawC','pfda.single.c','list'))
		} else {
			structure(single.c.core(y,Bt,subject,k,lm=penalties[1],lf=penalties[2],Kt,control$minimum.variance,control$max.iterations,control$convergence.tolerance)
				,class=c('pfda.single.c.R','pfda.single.c','list'))
		}
		rtn$tbase<-tbase
		rtn$y<-y
		rtn$subject<-subject
		return(rtn)
	}
}
print.pfda.single.c<-function(x,...){
	cat('Univariate Functional Principal Component Model\n')
	cat('Formula: ', deparse(attr(x,'formula')),'\n')
	cat(NCOL(x$tf),' principal components\n')
	cat('penalties are \n');print( penalty.pfda.single.c(x))
}
penalty.pfda.single.c<-function(x){
	with(x,c(mean=lm,pc=lf))
}
plot.pfda.single.c<-function(x,...){
	with(x,{
		layout(matrix(1:2,nrow=2,ncol=1,byrow=T))
		plot(tbase,tm, main=paste("Plot of mean curve for",attr(x,'name.t')),xlab=attr(x,'name.t'),ylab=attr(x,'name.y'))
		plot(tbase,tf, main=paste("Principle components for",attr(x,'name.t')),xlab=attr(x,'name.t'),ylab='')
	})
}
}
{ # single binaries
.roberts1<-function(c){
	a<-(c+sqrt(c^2+4))/2
	while(TRUE){
		u=runif(1)
		z=c-log(u)/a
		p=exp(-(z-a)^2/2)
		v=runif(1)
		if(v<=p)return(z)
	 }
}
.reject1<-function(c){
	while(TRUE){
		z=rnorm(1)
		if(z>c)return(z)
	}
}
.rtruncnormlower<-function(m,mean,sd,lower){
	c = (lower-mean)/sd
	if(c>0){
		replicate(m,.roberts1(c))
	} else {
		replicate(m,.reject1(c))
	}
}
.single.b.i<-function(y,B,subject,k,min.v){
	.single.c.i(y,B,subject,k,min.v)[1:5]
}
.single.b.1a.Sa.i<-function(B,tf,Da)solve(crossprod(B%*%tf)+solve(Da))
.single.b.1a.mu.i<-function(B,w, tm, tf, Sa)Sa%*%crossprod(B%*%tf, w-B%*%tm)
.single.b.1b.aa.i<-function(B,w,ww,tm,tf,Sa){
	phi=B%*%tf
	s1=-Sa%*%crossprod(phi,B%*%tm)
	s2=Sa%*%crossprod(phi,tcrossprod(as.matrix(w),s1))


	Sa%*%t(phi)%*%ww%*%phi%*%t(Sa) +
	s2 + t(s2) +
	tcrossprod(s1) +
	Sa
}
.single.b.1.i<-function(B,w,ww,tm,tf,Da){
	Sa = .single.b.1a.Sa.i(B,tf,Da)
	alpha = .single.b.1a.mu.i(B, w, tm, tf, Sa)
	aa = .single.b.1b.aa.i(B,w,ww,tm,tf,Sa)
	list(alpha=alpha,aa=aa,Sa=Sa)
}
.single.b.1<-function(B,subject,w,ww,tm,tf,Da){
	subject<-as.integer(as.factor(subject))
	n<-max(subject)
	k<-NCOL(tf)
	aa<-vector('list',n)
	alpha<-matrix(nrow=n,ncol=ncol(tf))
	Sa<-array(0,dim=c(k,k,n))
	for(i in seq(n)){
		ix<-i == as.integer(subject)
		r<-.single.b.1.i(B[ix,],w[ix],ww[[i]],tm,tf,Da)
		Sa[,,i]<-r$Sa
		alpha[i,]<-r$alpha
		aa[[i]]<-r$aa
	}
	list(alpha=alpha,aa=structure(unlist(aa),dim=c(1,1,10)),Sa=Sa)
}
.single.b.2<-function(B,subject,w,tf,alpha,lm,K){
	list(tm=.gen.tm(w,B,subject,tf,alpha,1,lm,K))
}
.single.b.3<-function(B,subject,w,tm,tf,alpha,aa,lf,K){
	list(tf=.gen.tf(w,B,subject,tm,tf,alpha,1,aa,lf,K))
}
.single.b.4<-function(tf,alpha,aa){
	saa<-matrix(0,NCOL(tf),NCOL(tf))
	for(sub in seq(NROW(alpha))) saa = saa + aa[[sub]]
	a<-.gen.orthog(tf,alpha,saa)
	list(tf=a$tf,Da=a$D, alpha = a$d)
}
.single.b.w.genw<-function(yi,wi,Bi,tm,tf,Da,kr,j){
	# R version of pfda_bin_s_gen_w in C code
	Bij <- Bi[-j,,drop=F]
	wij <- wi[-j]

	# pfda_bin_single_generate_w_parms1
	Ss <- solve(crossprod(Bij%*%tf)+solve(diag(Da,length(Da))))
	mu <- Ss%*%t(Bij%*%tf)%*%(wij-Bij%*%tm)

	# pfda_bin_single_generate_w_parms2
	a <- t(Bi[j,])%*%tm + t(Bi[j,])%*%tf%*%mu
	s <- 1+t(Bi[j,])%*%tf%*%Ss%*%t(tf)%*%Bi[j,]

	a+if(yi[j]) s*.rtruncnormlower(kr,0,1,-a/s) else -s*.rtruncnormlower(kr,0,1,a/s)
}
.single.b.w.1<-function(yi,Bi,wi,tm,tf,Da,kr){
	# R version of pfda_bin_single_approximate_moments_forobs
	n<-NROW(Bi)
	w_sim<-matrix(sapply(seq_len(n),.single.b.w.genw,kr=kr,yi=yi,Bi=Bi,wi=wi,tm=tm,tf=tf,Da=Da),nrow=kr)
	wi<-apply(w_sim,2,mean)
	# wwa<-array(apply(w_sim,2,tcrossprod),dim=c(n,n,kr))
	# wws<-matrix(0,n,n)
	# for(l in seq_len(kr))wws<-wws+wwa[,,l]
	wwi<-crossprod(w_sim)/kr
	return(list(wi=wi,wwi=wwi))
}
.single.b.w<-function(y,B,subject,w,ww,tm,tf,Da,weight,kr){
	for(i in seq_len(nlevels(subject))){
		ix    <- i==as.integer(subject)
		rtn     <- .single.b.w.1(y[ix],B[ix,,drop=F],w[ix],tm,tf,Da,kr)
		w[ix] <- (1-weight)*w[ix]+weight*rtn[[1]]
		ww[[i]] <- (1-weight)*ww[[i]]+weight*rtn[[2]]
	}
	return(list(w=w,ww=ww))
}
.single.b.core<-function(y,B,subject,k,lm,lf,K,minimum.variance, max.iterations,convergence.tolerance){
	{ # setup variables & initial values
		M<-NROW(y)
		N<-nlevels(subject)
		nobs<-table(subject)
		w <- y
		ww<-vector('list',N)
			for(i in seq_len(N)) ww[[i]]<-matrix(0,nobs[i],nobs[i])
		Saa<-array(0,dim=c(k,k,N))
		rtn<-.single.b.i(y,B,subject,k,minimum.variance)
		tm<-rtn$tm
		tf<-rtn$tf
		Da<-rtn$Da
		alpha<-rtn$alpha
		old<-vector('list',0)
		aa <- rtn$aa	#vector('list',N); for(i in seq_len(N)) aa<-matrix(0,k,k)
		I=0
		.cc<-numeric(0)
	}
	while(I<-I+1){
		{ #step w
			r0=k0=100
			kr=10
			if(I<r0){
				weight=1
				rtn<-.single.b.w(y,B,subject,w,ww,tm,tf,Da,k0,weight)
			} else {
				weight=10/(10+I)
				rtn<-.single.b.w(y,B,subject,w,ww,tm,tf,Da,kr,weight)
			}
			w<-rtn$w
			ww<-rtn$ww
		}
		s1 <- .single.b.1(B,subject,w,ww,tm,tf,Da)
		s2 <- .single.b.2(B,subject,w,tf,s1$alpha,lm,K)
		s3 <- .single.b.3(B,subject,w,s2$tm,tf,s1$alpha,s1$aa,lf,K)
		s4 <- .single.b.4(tf,alpha,aa)
		{ #check convergence
			ccl<-numeric(0)
			ccl['tm'] <- sum(abs((tm-s2$tm)/s2$tm))
			ccl['tf'] <- sum(abs((tf-s4$tf)/s4$tf))
			ccl['Da'] <- sum(abs((Da-s4$Da)/s4$Da))
			cc <- sum(ccl)
			if(cc<convergence.tolerance)break
			if(I>=max.iterations){
				warning('Maximum number of iterations exceeded, convergence not obtained.')
				break
			}
			.cc<-c(.cc,cc)
			{ #reassign variable
				tm    <- s2$tm
				tf    <- s4$tf
				Da    <- s4$Da
				alpha <- s4$alpha
				aa    <- s1$aa
				Saa   <- s1$Sa
			}
		}
	}
	list(tm=tm,tf=tf,Da=Da,aa=aa,Saa=Saa,I=I,cc=cc)
}
.single.b.updatePCS<-function(Yi, Da, phii, Y.rhoi, alpha){
	ni = length(Yi)

	muy = Y.rhoi+phii%*%alpha
	mua = pnorm(muy)
	Y.tilda = muy + (Yi-mua)/dnorm(mua)

	A = tcrossprod(diag(Da,length(Da)), phii)
	Q = solve(diag(1, ni) + phii%*%A)

	list(alpha= as.vector((A %*% Q) %*% Y.tilda))
}
.single.b.estimatePCS<-function(Y, B, subject, tm, tf, Da){
	Y.rho = B%*%tm
	phi = B%*%tf
	ofun<-function(i){
		ix = as.integer(subject)==i
		alpha = rep(0,NCOL(tf))
		dif = 1
		while(dif>1e-3){
			ab<-.single.b.updatePCS(Y[ix], Da, phi[ix,,drop=F], Y.rho[ix,,drop=F], alpha)
			dif<-mean(abs(alpha-ab$alpha))
			alpha = ab$alpha
		}
		alpha
	}
	alpha<-matrix(unlist(lapply(seq_len(nlevels(subject)), ofun)),nlevels(subject),NCOL(tf),byrow=T)
	alpha
}
.single.b.n2L.1<-function(yi,Bi,tm,tf,alpha,Da){
	sum(.single.b.n2L.parts(yi,Bi,tm,tf,alpha,Da))
	}
.single.b.n2L.parts<-function(yi,Bi,tm,tf,alpha,Da){
	phi = Bi%*%tf
	mui = pnorm(Bi%*%tm+phi%*%alpha)
	nui = mui*(1-mui)
	Wi = as.vector(dnorm(mui)^2/nui)
	c(sum(log(Da)),
		determinant(diag(NROW(Da))+crossprod(phi,diag(Wi)%*%phi%*%diag(Da,length(Da))))$modulus, 
		crossprod(alpha,solve(diag(Da,length(Da)))%*%alpha), 
		-2*sum(yi*log(mui)+ (1-yi)*log(1-mui)-Wi-mui)
	)
}
.single.b.n2L.bySubject<-function(y, subject, B, tm, tf, Da, alpha){
	sapply(seq_len(nlevels(subject)),function(si){
			ix=as.integer(subject)==si
			.single.b.n2L.1(y[ix], B[ix,], tm, tf, alpha[si,], Da)
		})
}
.single.b.n2L<-function(y, subject, B, tm, tf, Da){
	alpha<-.single.b.estimatePCS(y, B, subject, tm, tf, Da)
	ll<-sum(sapply(seq_len(nlevels(subject)),function(si){
		ix=as.integer(subject)==si
		.single.b.n2L.1(y[ix], B[ix,], tm, tf, alpha[si,], Da)
	}))
}
loglik.pfda.single.b<-function(object,...,newdata=NULL,n2L=TRUE){
	y<-object@FittingData[[1]][[1]]
	t<-object@FittingData[[2]]
	subject<-object@FittingData[[3]]
	B<-evaluate(object@Basis,t)


}
AIC.pfda.single.b<-function(object,...){
	n2L<-with(object,.single.b.n2L(y, subject, Bt, tm, tf, Da))
	n2L + with(object,2*(.pfda.df(Bt,lm,K)+k*.pfda.df(Bt,lf,K)))
}
single.b<-function(y,t,subject, knots=NULL, penalties=NULL, df=NULL, k=NULL, control=pfdaControl(),subset=NULL){
	{ # setup
	fname = deparse(match.call()[[1L]])
	localfuncs('.F.single.optimize.npc')
	eval(.X.subset)
	eval(.X.single.knots)
	eval(.X.single.penalties)
	}
	if(is.null(k)||any(is.na(k))){
		# stop('number of principal components optimization not finished yet')
		.F.single.optimize.npc()
	} else
	if(any(is.na(penalties))) {
		funcall <- match.call()
		eval(.X.optimize.penalties)
	}
	else {
		eval(.X.binary.y)
		rtn<-if(control$useC){
			{ # setup for passing to Compiled code
				k  = as.integer(k[1])
				N   = nlevels(subject)
				nobs= table(subject)
				M   = length(y)
				kr = with(control,max(binary.k0,binary.kr))
				p  = ncol(Bt)
				{ #Compute Memory Requirements
					ni = max(nobs)
					.inits                      = 2*p^2 + M + p*N + 8*p
					.step.W                	    = ni*kr + kr + (ni * k + ni + k^2 + k + p + p^2 + k*p + 3*k)
					.step.1.E              	    = M + M*k + N*k + 2* k^2 + 2 * k + k *ni
					.step.2                 	= p^2+ M+ M*k
					.step.3                 	= M + p^2 + p + k^2
					.step.4                  	= k + 2*k^2 + 2*p^2 + p + p*max(k,8)
					dpl <- M + sum(nobs^2) + N*k^2 + N*p^2 + p + p*k + k + N + max(.inits, .step.W, .step.1.E, .step.2, .step.3, .step.4)
					ipl <- 8*p
				}
			}
			structure(.C('pfda_bin_single', y=y, nobs=nobs, M=M, N=N, k=k, Bt=Bt, p=p, lm=penalties[1], lf=penalties[2], K=Kt, tm=double(p), tf=matrix(0,p,k), Da=double(k), alpha=matrix(0,N,k), Saa=array(0,dim=c(k,k,N)),  minV=control$minimum.variance,  tol=control$convergence.tolerance,  Iterations=control$max.iterations, burninlength=as.integer(control$binary.burnin), burningenerate=as.integer(control$binary.k0), weightedgenerate=as.integer(control$binary.kr), dl=control$C.debug, dp=double(dpl) , p=integer(ipl))
				,class=c('pfda.single.b.rawC','pfda.single.b','list'))
		} else {
			structure(.single.b.core(y,Bt,subject,k,penalties[1],penalties[2],K=Kt,control$minimum.variance, control$max.iterations,control$convergence.tolerance)
				,class=c('pfda.single.b.R','pfda.single.b','list'))
		}
		rtn$tbase<-tbase
		rtn$y<-y
		rtn$subject<-subject
		return(rtn)
	}
}
print.pfda.single.b<-function(x,...){
	cat('Univariate Binary Functional Principal Component Model\n')
	cat('Formula: ', deparse(attr(x,'formula')),'\n')
	cat(NCOL(x$tf),' principal components\n')
	cat('penalties are \n');print( penalty.pfda.single.b(x))
}
penalty.pfda.single.b<-function(x){
	with(x,c(mean=lm,pc=lf))
}
}
{ # Dual (Continuous/Continuous) case
.dual.cc.i<-function(y,z,B,subject,ka,kb,min.v){
	a<-.single.c.i(y,B,subject,ka,min.v)
	b<-.single.c.i(z,B,subject,kb,min.v)
	lambda=crossprod(b$alpha,a$alpha)%*%solve(crossprod(a$alpha))
	ab<-array(0,dim=c(ka,kb,nlevels(subject)))
	for(i in seq_len(nlevels(subject)))ab[,,i]<-tcrossprod(a$alpha[i,],b$alpha[i,])

	{ list(tm=a$tm,tn=b$tm,
		tf=a$tf,tg=b$tf,
		alpha=a$alpha,beta=b$alpha,
		Da=a$Da,Db=b$Da,
		lambda=lambda,
		s.eps=a$sigma,s.xi=b$sigma,
		aa=a$aa,ab=ab,bb=b$aa
	)}
}
.dual.cc.1.1<-function(yi,zi,Bi,tm,tn,tf,tg,lambda,Da,Db,s.eps,s.xi,Setai){
	phi = Bi%*%tf
	psi = Bi%*%tg

	Saa = diag(1/Da,length(Da)) + t(lambda)%*%Setai%*%lambda+crossprod(phi)/s.eps
	Sab = crossprod(-lambda,Setai)
	Sbb = Setai+crossprod(psi)/s.xi

	Sbbi = solve(Sbb)

	Zaa = solve(Saa-Sab%*%Sbbi%*%t(Sab))
	Zab = -Zaa%*%Sab%*%Sbbi
	Zbb = solve(Sbb-t(Sab)%*%solve(Saa)%*%Sab)

	mu.a = Zaa%*%crossprod(phi,yi-Bi%*%tm)/s.eps + Zab%*%crossprod(psi,zi-Bi%*%tn)
	mu.b = t(Zab)%*%crossprod(phi,yi-Bi%*%tm)/s.eps + Zbb%*%crossprod(psi,zi-Bi%*%tn)

	aa = tcrossprod(mu.a)+Zaa
	ab = tcrossprod(mu.a,mu.b)+Zab
	bb = tcrossprod(mu.b)+Zbb
	list(alpha=mu.a, beta=mu.b, Saa=Zaa, Sab=Zab, Sbb = Zbb,aa=aa,ab=ab,bb=bb)
}
.dual.cc.1<-function(y,z,B,subject,tm,tn,tf,tg,lambda,Da,Db,s.eps,s.xi){
	Setai = solve(Db-lambda%*%Da%*%t(lambda))
	#TODO
	{
		Saa = array(0,dim=c(NCOL(tf),NCOL(tf),nlevels(subject)))
		Sab = array(0,dim=c(NCOL(tf),NCOL(tg),nlevels(subject)))
		Sbb = array(0,dim=c(NCOL(tg),NCOL(tg),nlevels(subject)))
		alpha = matrix(nrow=nlevels(subject),ncol=NCOL(tf))
		beta  = matrix(nrow=nlevels(subject),ncol=NCOL(tg))
		aa = array(0,dim=c(NCOL(tf),NCOL(tf),nlevels(subject)))
		ab = array(0,dim=c(NCOL(tf),NCOL(tg),nlevels(subject)))
		bb = array(0,dim=c(NCOL(tg),NCOL(tg),nlevels(subject)))
	}
	for(i in seq_len(nlevels(subject))){
		ix = i == as.integer(subject)
		r<-.dual.cc.1.1(y[ix],z[ix],B[ix,],tm,tn,tf,tg,lambda,Da,Db,s.eps,s.xi,Setai)
		{  # assign values
			Saa[,,i] <-r$Saa
			Sab[,,i] <-r$Sab
			Sbb[,,i] <-r$Sbb
			alpha[i,]<-r$alpha
			beta[i,] <-r$beta
			aa[,,i]  <-r$aa
			ab[,,i]  <-r$ab
			bb[,,i]  <-r$bb
		}
	}
	list(alpha=alpha,beta=beta,Saa=Saa,Sab=Sab,Sbb=Sbb,aa=aa,ab=ab,bb=bb)
}
.dual.cc.2<-function(y,z,B,subject,tm,tn,tf,tg,alpha,beta,Saa,Sbb){
	list(
	s.eps = .single.c.2(y,B,subject,tm,tf,alpha,Saa)$sigma,
	s.xi  = .single.c.2(z,B,subject,tn,tg,beta ,Sbb)$sigma)
}
.dual.cc.3<-function(y,z,B,subject,tf,tg,alpha,beta,s.eps,s.xi,lm,ln,K){
	list(
	tm = .gen.tm(y,B,subject,tf,alpha,s.eps,lm,K),
	tn = .gen.tm(y,B,subject,tg,beta ,s.xi ,ln,K))
}
.dual.cc.4<-function(y,z,B,subject,tm,tn,tf,tg,alpha,beta,s.eps,s.xi,aa,bb,lf,lg,K){
	list(
	tf = .gen.tf(y,B,subject,tm,tf,alpha,s.eps,aa,lf,K),
	tg = .gen.tf(z,B,subject,tn,tg,beta ,s.xi ,bb,lg,K))
}
.dual.cc.5<-function(aa,ab){
	sum.aa<-matrix(0,dim(aa)[1],dim(aa)[2])
	sum.ab<-matrix(0,dim(ab)[1],dim(ab)[2])
	for(i in seq_len(dim(aa)[3])){
		sum.aa = sum.aa + aa[,,i]
		sum.ab = sum.ab + ab[,,i]
	}
	list(lambda = .gen.dual.lambda(sum.aa,sum.ab))
}
.dual.cc.6<-function(tf,tg,alpha,beta,lambda,aa,bb){
	sum.aa<-matrix(0,dim(aa)[1],dim(aa)[2])
	sum.bb<-matrix(0,dim(bb)[1],dim(bb)[2])
	for(sub in seq(NROW(alpha))){
		sum.aa = sum.aa + aa[,,sub]
		sum.bb = sum.bb + bb[,,sub]
	}
	a<-.gen.orthog(tf,alpha,sum.aa)
	b<-.gen.orthog(tg,beta ,sum.bb)

	list(tf = a$tf, tg = b$tf, Da=a$D, Db=b$D, alpha = a$d, beta = b$d,
		lambda = b$transformation%*%lambda%*%solve(a$transformation))
}
.dual.cc.n2L<-function(t, y, z, B, subject, tm, tn, tf, tg, lambda, Da, Db, seps, sxi){
	ka<-length(Da)
	kb<-length(Db)
	phi<-B%*%tf
	psi<-B%*%tg
	ll<-0
	for(i in seq_len(nlevels(subject))){
		ind<-i==as(subject,"integer")
		ni<-sum(ind)
		Sy<-phi[ind,,drop=FALSE]%*%diag(x=Da,ka,ka)%*%t(phi[ind,,drop=FALSE])+diag(seps,sum(ind))
		Sz<-psi[ind,,drop=FALSE]%*%diag(x=Db,kb,kb)%*%t(psi[ind,,drop=FALSE])+diag(sxi,sum(ind))
		Syz<-phi[ind,,drop=FALSE]%*%diag(Da,ka,ka)%*%t(lambda)%*%t(psi[ind,,drop=FALSE])
		Ry<-y[ind]-B[ind,]%*%tm
		Rz<-z[ind]-B[ind,]%*%tn

		Cmat<-rbind(cbind(Sy,Syz),cbind(t(Syz),Sz))
		Cimat<-solve(Cmat)
		Sy.inv <-Cimat[seq_len(ni),seq_len(ni)]			#solve(Sy-crossprod(Syz, solve(Sz,t(Syz))))
		Sz.inv <- Cimat[seq_len(ni)+ni,seq_len(ni)+ni]			#Szi+Szi%*%t(Syz)%*%Sy.inv%*%Syx%*%Szi		#solve(Sz-tcrossprod(Syz,solve(Sy,t(Syz))))
		Syz.inv<- Cimat[seq_len(ni),seq_len(ni)+ni]			#-solve(Sy, Syz)%*%Sz.inv

		Dety<-determinant(Sy,logarithm=TRUE)
		Detzinv<-determinant(Sz.inv,logarithm=TRUE)
		if((Dety$sign<0)|(Detzinv$sign<0))stop("Negative determinant encountered for variance matrix.")
		ll<- ll + Dety$modulus-Detzinv$modulus
		ll<- ll + crossprod(Ry,Sy.inv%*%Ry)
		ll<- ll + crossprod(Rz,Sz.inv%*%Rz)
		ll<- ll + 2*crossprod(Ry,Syz.inv%*%Rz)
	}
	ll
}
AIC.pfda.dual.cc<-function(object,...){
	with(object,.dual.ge.AIC(
		.dual.cc.n2L(t, y, z, B, subject, tm, tn, tf, tg, lambda, Da, Db, seps, sxi),0,
		B   , B   , B  , B ,
		1   , ka  , 1  , kb,
		penalties[1]  , penalties[2]  , penalties[3] , penalties[4],
		K   , K   , K  , K ,
		seps, seps, sxi, sxi)
	)
	# .dual.cc.AIC.rawC(object, t, y, z, B, subject, lm, ln, lf, lg, K))
}
logLik.pfda.dual.cc<-function(object,...,newdata=NULL,n2L=TRUE){
	with(object,with(newdata,.dual.cc.n2L(t,y,z,B,subject,tm,tn,tf,tg,lambda,Da,Db,seps,sxi)))
}
dual.cc.core<-function(y,z,B,subject,ka,kb,lm,ln,lf,lg,K,min.v,max.I,tol){
	{ #initial values
		init<-.dual.cc.i(y,z,B,subject,ka,kb,min.v)
		tm<-init$tm
		tn<-init$tn
		tf<-init$tf
		tg<-init$tg
		alpha<-init$alpha
		beta<-init$beta
		aa<-init$aa
		ab<-init$ab
		bb<-init$bb
		Da<-init$Da
		Db<-init$Db
		lambda<-init$lambda
		n<-nlevels(subject)
		Saa<-vector('list',n); for(i in seq_len(n))Saa[[i]]<-matrix(0,ka,ka)
		Sab<-vector('list',n); for(i in seq_len(n))Sab[[i]]<-matrix(0,ka,kb)
		Sbb<-vector('list',n); for(i in seq_len(n))Sbb[[i]]<-matrix(0,kb,kb)
		s.eps<-as.vector(init$s.eps)
		s.xi<-as.vector(init$s.xi)
		I<-0
		nobs<-table(subject)
		.cc<-numeric(0)
	}
	while(I<-I+1){
		s1<-.dual.cc.1(y,z,B,subject,tm,tn,tf,tg,lambda,Da,Db,s.eps,s.xi)
		s2<-.dual.cc.2(y,z,B,subject,tm,tn,tf,tg,s1$alpha,s1$beta,s1$Saa,s1$Sbb)
		s3<-.dual.cc.3(y,z,B,subject,tf,tg,alpha,beta,s2$s.eps,s2$s.xi,lm,ln,K)
		s4<-.dual.cc.4(y,z,B,subject,s3$tm,s3$tn,tf,tg,s1$alpha,s1$beta,s2$s.eps,s2$s.xi,s1$aa,s1$bb,lf,lg,K)
		s5<-.dual.cc.5(s1$aa,s1$ab)
		s6<-.dual.cc.6(s4$tf,s4$tg,s1$alpha,s1$beta,s5$lambda,s1$aa,s1$bb)
		{ # convergence
			ccl<-numeric(0)
			ccl['tm']<-sum(abs(tm-s3$tm))
			ccl['tn']<-sum(abs(tn-s3$tn))
			ccl['tf']<-sum(abs(tf-s6$tf))
			ccl['tg']<-sum(abs(tf-s6$tg))
			ccl['Da']<-sum(abs(Da-s6$Da))
			ccl['Db']<-sum(abs(Da-s6$Db))
			ccl['lambda']<-sum(abs(lambda-s6$lambda))
			(cc<-sum(ccl))
			{ # reassign values
				aa<-s1$aa
				ab<-s1$ab
				bb<-s1$bb
				Saa<-s1$Saa
				Sab<-s1$Sab
				Sbb<-s1$Sab
				s.eps<-s2$s.eps
				s.xi<-s2$s.xi
				tm<-s3$tm
				tn<-s3$tn
				tf<-s6$tf
				tg<-s6$tg
				Da<-s6$Da
				Db<-s6$Db
				alpha<-s6$alpha
				beta<-s6$beta
				lambda<-s6$lambda
			}
			if(cc<tol)break
			if(I>max.I){
				warning('Maximum number of iterations exceeded, convergence not obtained.')
				break
			}
			.cc<-c(.cc,cc)
		}
	}
	list(tm=tm,tn=tn,tf=tf,tg=tg,alpha=alpha,beta=beta,lambda=lambda,Da=Da,Db=Db,s.eps=s.eps,s.xi=s.xi)
}
dual.cc<-function(y,z,t,subject, knots=NULL, penalties=NULL,df=NULL, k=NULL, control=pfdaControl(),tbase=NULL, subset=NULL){
	fname = deparse(match.call()[[1L]])
	eval(.X.subset)
	eval(.X.dual.k)
	eval(.X.single.knots)
	eval(.X.dual.penalties)
	if(any(is.na(k))){
		# stop("not finished with number of principal component optimization.")
		if(is.na(ka)){
			model.y.single <- single.c(y,NULL,t,subject, knots=knots, penalties=penalties[,1],k=NULL, control=control,subset=subset)
			ka <- model.y.single$k
		}
		if(is.na(kb)){
			model.z.single <- single.c(z,NULL,t,subject, knots=knots, penalties=penalties[,2],k=NULL, control=control,subset=subset)
			kb<-model.z.single$k
		}
		Recall(y,z,t,subject, knots=knots,penalties=penalties,k=c(ka,kb), control=control)
	}
	else if(any(is.na(penalties))){
		# funcall <- match.call()
		# eval(.X.optimize.penalties)}
		localfuncs('.F.optimize.penalties')
		.F.optimize.penalties()}
	else { # Pass to core function
		rtn<-{ modifyList(
			if(control$useC){
				eval(.X.read.parameters)
				{ #Compute Memory Requirements
					k<-max(ka,kb)
					dpl_1     <- M + M*ka + 2*k^2
					dpl_2     <- p^2 + M + M*k
					dpl_3     <- M + p^2 + p + k^2
					dpl_4     <- k^2 + ka*kb
					dpl_5_1   <- k + 2*p^2 + p +
						max(  outer_qf = k*p , eigens = 8*p)
					dpl_5_2   <- ka^2 + kb*ka
					dpl_5     <- ka^2 + kb^2 +
						max(dpl_5_1  ,dpl_5_2 )
					dpl_E_1   <- kb^2 + ka^2 +
						max( outer_qf = ka*kb,
								 sym_inv = 2* kb^2,
								 inner_qf = ka*kb )
					dpl_E_2_1 <- 2*k^2
					dpl_E_2_2 <- 3*k^2
					dpl_E_2   <- kb^2 + max(dpl_E_2_1,dpl_E_2_2,ka*kb)
					dpl_E_3_1 <- 2*k*max(nobs)
					dpl_E_3   <- dpl_E_3_1
					dpl_E     <- 2*M + M*ka + M*kb + ka^2 + kb^2 + ka*kb +
						max(dpl_E_1, dpl_E_2, dpl_E_3)
					dpl<- N*(p**2) + p*2 + p*ka+p*kb + ka*kb + ka + kb +
						max(dpl_1, dpl_2, dpl_3, dpl_4, dpl_5, dpl_E )
					ipl<-8*p
				}
				{ structure(.C("pfdaDual",
					residuals.y=as.double(y), residuals.z=as.double(z),
					B=(Bt),
					tm=double(p),			tn=double(p),
					tf=matrix(0,p,ka),			tg=matrix(0,p,kb),
					alpha =matrix(0,N,ka),
					beta =matrix(0,N,kb),
					lambda =matrix(0,kb,ka),
					Da =double(ka),
					Db =double(kb),
					seps =	double(1),
					sxi =double(1),
					Saa =array(0,c(ka,ka,N)),
					Sab =array(0,c(ka,kb,N)),
					Sbb =array(0,c(kb,kb,N)),
					nobs=as.integer(nobs),
					M=as.integer(M),
					N=N,
					ka=as.integer(ka),
					kb=as.integer(kb),
					p=as.integer(p),
					penalties = penalties,
					K = Kt,
					minV=control$minimum.variance,
					Iterations=control$max.iterations,
					tol=control$convergence.tolerance,
					dl=control$C.debug,
					dp=double(dpl),
					ip=integer(ipl)
				),class=c('pfda.dual.cc.rawC','pfda.dual.cc','list'))}
			} else {
				structure(dual.cc.core(y,z,B,subject,ka,kb,lm,ln,lf,lg,K,min.v,max.I,tol)
					,class=c('pfda.dual.cc.R','pfda.dual.cc','list'))
			},
		list(tbase=tbase,subject=subject,y=y,z=z))}
	}
}
print.pfda.dual.cc<-function(x,...){
	cat('Paired Functional Principal Component Model (Continuous/Continuous)\n')
	cat('Formula: ', deparse(attr(x,'formula')),'\n')
	cat(attr(x,'name.y'),' has ', NCOL(x$tf),' principal components\n')
	cat(attr(x,'name.z'),' has ', NCOL(x$tg),' principal components\n')
	cat('penalties are \n');print( penalty.pfda.dual.cc(x))
}
penalty.pfda.dual.cc<-function(object,..)with(object,structure(matrix(penalties,2,2),dimnames=list(c(attr(object,'name.y'),attr(object,'name.z')),c('mean','pc'))))
}
{ # Dual(Binary/Continuous) case
.dual.bc.i<-function(y,z,B,subject,ka,kb,min.v){
	a<-.single.b.i(y,B,subject,ka,min.v)
	b<-.single.c.i(z,B,subject,kb,min.v)
	lambda=crossprod(b$alpha,a$alpha)%*%solve(crossprod(a$alpha))

	list(tm=a$tm,tn=b$tm,tf=a$tf,tg=b$tf,alpha=a$alpha,beta=b$alpha,aa=a$aa,bb=b$aa,Da=a$Da,Db=b$Da,lambda=lambda,s.xi=b$sigma)
}
.dual.bc.1a<-function(zi,Bi,wi,tm,tn,tf,tg,s.xi,Saa,Sab){
	crossprod(Bi%*%tf%*%Saa,wi-Bi%*%tm)+Sab%*%crossprod(Bi%*%tg,zi-Bi%*%tn)/s.xi
}
.dual.bc.1b<-function(zi,Bi,wi,tm,tn,tf,tg,s.xi,Sab,Sbb){
	crossprod(Bi%*%tf%*%Sab,wi-Bi%*%tm)+Sbb%*%crossprod(Bi%*%tg,zi-Bi%*%tn)/s.xi
}
.dual.bc.1cde<-function(zi,Bi,wi,wwi,tm,tn,tf,tg,s.xi,Saa,Sab,Sbb){
	S1<-.dual.bc.1a(zi,Bi,0,tm,tn,tf,tg,s.xi,Saa,Sab)
	S2<-.dual.bc.1b(zi,Bi,0,tm,tn,tf,tg,s.xi,Sab,Sbb)
	S3<-tcrossprod(Saa,Bi%*%tf)
	S4<-S3%*%wi
	S5<-t(Bi%*%tf%*%Sab)
	S6<-S5%*%wi

	list(
		aa= Saa+tcrossprod(S1)+tcrossprod(S1,S4)+tcrossprod(S4,S1) + S3%*%wwi%*%t(S3)
		,
		ab= Sab+tcrossprod(S1,S2)+tcrossprod(S4,S2)+tcrossprod(S1,S6)+S3%*%tcrossprod(wwi,S5)
		,
		bb= Sbb+tcrossprod(S2)+tcrossprod(S6,S2)+tcrossprod(S2,S6)+S5%*%tcrossprod(wwi,S5)
	)
}
.dual.bc.1<-function(z,B,subject,w,ww,tm,tn,tf,tg,lambda,Da,Db,s.xi){
	m<-NROW(B)
	n<-nlevels(subject)
	ka<-NCOL(tf)
	kb<-NCOL(tg)

	alpha<-matrix(nrow=n,ncol=ka)
	beta <-matrix(nrow=n,ncol=kb)
	aa <- array(dim=c(ka,ka,n))
	ab <- array(dim=c(ka,kb,n))
	bb <- array(dim=c(kb,kb,n))
	Saa<-vector('list',n)
	Sab<-vector('list',n)
	Sbb<-vector('list',n)
	for(sn in seq_len(n)){
		ix<-sn==as.integer(subject)
		Bi<-B[ix,]
		s<-.gen.dual.sigmas(Bi,Bi,tf,tg,lambda,Da,Db,1,s.xi)
		alpha[sn,]<-.dual.bc.1a(z[ix],Bi,w[ix],tm,tn,tf,tg,s.xi,s$Saa,s$Sab)
		beta [sn,]<-.dual.bc.1b(z[ix],Bi,w[ix],tm,tn,tf,tg,s.xi,s$Sab,s$Sbb)
		cde<-.dual.bc.1cde(z[ix],Bi,w[ix],ww[[sn]],tm,tn,tf,tg,s.xi,s$Saa,s$Sab,s$Sbb)
		aa[,,sn]<-cde$aa
		ab[,,sn]<-cde$ab
		bb[,,sn]<-cde$bb
		Saa[[sn]]<-s$Saa
		Sab[[sn]]<-s$Sab
		Sbb[[sn]]<-s$Sbb
	}
	list(alpha=alpha,beta=beta,aa=aa,ab=ab,bb=bb,Saa=Saa,Sab=Sab,Sbb=Sbb)
}
.dual.bc.2<-function(z,B,subject,tn,tg,beta,Sbb){
	m<-NROW(B)
	Rz<-z-B%*%tn
	sum.xi<-0
	for(i in seq_len(nlevels(subject))){
		ix<- i == as.integer(subject)
		sum.xi<-sum.xi+sum(diag(B[ix,]%*%tg%*%tcrossprod(Sbb[[i]],B[ix,]%*%tg)))
		Rz[ix]<-Rz[ix]-B[ix,]%*%tg%*%beta[i,]
	}
	list(s.xi=(sum.xi+crossprod(Rz))/m)
}
.dual.bc.3<-function(z,B,subject,w,tf,tg,alpha,beta,s.xi,lm,ln,K){
	list(
	tm=.gen.tm(w,B,subject,tf,alpha,1   ,lm,K),
	tn=.gen.tm(z,B,subject,tg,beta ,s.xi,ln,K)
	)
}
.dual.bc.4<-function(z,B,subject,w,tm,tn,tf,tg,alpha,beta,s.xi,aa,bb,lf,lg,K){
	list(
	tf=	.gen.tf(w,B,subject,tm,tf,alpha,1   ,aa,lf,K),
	tg= .gen.tf(z,B,subject,tn,tg,beta ,s.xi,bb,lg,K))
}
.dual.bc.5<-function(aa,ab){
	.dual.cc.5(aa,ab)
}
.dual.bc.6<-function(tf,tg,alpha,beta,lambda,aa,bb){
	.dual.cc.6(tf,tg,alpha,beta,lambda,aa,bb)
}
.dual.bc.genw<-function(j,kr,yi,wi,zi,Bi,tm,tn,tf,tg,Da,Db,lambda){
	#
	Bij <- Bi[-j,,drop=F]
	wij <- wi[-j]

	#
	r<-.gen.dual.sigmas(Bij,Bi,tf,tg,lambda,Da,Db,1,1)
	mu <- r$Saa %*%crossprod(Bij%*%tf,wij-Bij%*%tm) + r$Sab %*% crossprod(Bi%*%tg, zi-Bi%*%tn)
	# mu = .dual.bc.1a(zi,Bij,wij,tm,tn,tf,tg,1.0,r$Saa,r$Sab)

	# Ss <- solve(crossprod(Bij%*%tf)+solve(diag(Da,length(Da))))

	# pfda_bin_single_generate_w_parms2
	a <- crossprod(Bi[j,],tm) + t(Bi[j,])%*%tf%*%mu
	s <- 1+ t(Bi[j,])%*%tf%*%r$Saa%*%t(tf)%*%Bi[j,]

	a+if(yi[j]){
		s*.rtruncnormlower(kr,0,1,-a/s)
	} else {
		-s*.rtruncnormlower(kr,0,1,a/s)
	}
}
.dual.bc.w.1<-function(kr,yi,wi,zi,Bi,tm,tn,tf,tg,Da,Db,lambda){
	# R version of
	n<-NROW(Bi)
	w_sim<-sapply(seq_len(n),.dual.bc.genw,kr,yi,wi,zi,Bi,tm,tn,tf,tg,Da,Db,lambda)
	wi<-apply(w_sim,2,mean)
	wwi<-crossprod(w_sim)/kr
	list(wi=wi,wwi=wwi)
}
.dual.bc.w<-function(y,z,B,subject,w,ww,tm,tn,tf,tg,lambda,Da,Db,weight,kr){
	for(i in seq_len(nlevels(subject))){
		ix      <- i==as.integer(subject)
		rtn     <- .dual.bc.w.1(kr,y[ix],w[ix],z[ix],B[ix,,drop=F],tm,tn,tf,tg,Da,Db,lambda)
		w[ix]   <- (1-weight)*w[ix]+weight*rtn[[1]]
		ww[[i]] <- (1-weight)*ww[[i]]+weight*rtn[[2]]
	}
	return(list(w=w,ww=ww))
}
dual.bc.core<-function(y,z,B,subject,ka,kb,lm,ln,lf,lg,K,min.v,max.I,tol){
	{ #initial values
		init<-.dual.bc.i(y,z,B,subject,ka,kb,min.v)
		tm<-init$tm
		tn<-init$tn
		tf<-init$tf
		tg<-init$tg
		alpha<-init$alpha
		beta<-init$beta
		aa<-init$aa
		ab<-init$ab
		bb<-init$bb
		Da<-init$Da
		Db<-init$Db
		lambda<-init$lambda
		n<-nlevels(subject)
		Saa<-vector('list',n); for(i in seq_len(n))Saa[[i]]<-matrix(0,ka,ka)
		Sab<-vector('list',n); for(i in seq_len(n))Sab[[i]]<-matrix(0,ka,kb)
		Sbb<-vector('list',n); for(i in seq_len(n))Sbb[[i]]<-matrix(0,kb,kb)
		s.xi<-as.vector(init$s.xi)
		I<-0
		w<-as.double(y)
		nobs<-table(subject)
		ww<-vector('list',n); for(i in seq_len(n))ww[[i]]<-matrix(0,nobs[i],nobs[i])
		r0=k0=100 ##TODO fix this to read from control values.
		kr=10
		.cc<-numeric(0)
	}
	while(I<-I+1){
		{ # step w
			if(I<r0){
				weight=1
				rtn<-.dual.bc.w(y,z,B,subject,w,ww,tm,tn,tf,tg,lambda,Da,Db,weight,k0)
			} else {
				weight=10/(10+I)
				rtn<-.dual.bc.w(y,z,B,subject,w,ww,tm,tn,tf,tg,lambda,Da,Db,weight,kr)
			}
			w<-rtn$w
			ww<-rtn$ww
		}
		s1<-.dual.bc.1(z,B,subject,w,ww,tm,tn,tf,tg,lambda,Da,Db,s.xi)
		s2<-.dual.bc.2(z,B,subject,tn,tg,s1$Saa)
		s3<-.dual.bc.3(z,B,subject,w,tf,tg,s1$alpha,s2$s.xi,s1$beta,lm,ln,K)
		s4<-.dual.bc.4(z,B,w,s3$tm,s3$tn,tf,tg,s1$alpha,s1$beta,s2$s.xi,s1$aa,s1$bb,lf,lg,K)
		s5<-.dual.bc.5(s1$aa,s1$ab)
		s6<-.dual.bc.6(s4$tf,s4$tg,s1$alpha,s1$beta,s5$lambda,s1$aa,s1$bb)
		{ # convergence
			ccl<-numeric(0)
			ccl['tm']<-sum(abs(tm-s3$tm))
			ccl['tn']<-sum(abs(tn-s3$tn))
			ccl['tf']<-sum(abs(tf-s6$tf))
			ccl['tg']<-sum(abs(tf-s6$tg))
			ccl['Da']<-sum(abs(Da-s6$Da))
			ccl['Db']<-sum(abs(Da-s6$Db))
			ccl['lambda']<-sum(abs(lambda-s6$lambda))
			(cc<-sum(ccl))
			{ # reassign values
				aa<-s1$aa
				ab<-s1$ab
				bb<-s1$bb
				Saa<-s1$Saa
				Sab<-s1$Sab
				Sbb<-s1$Sab
				s.xi<-s2$s.xi
				tm<-s3$tm
				tn<-s3$tn
				tf<-s6$tf
				tg<-s6$tg
				Da<-s6$Da
				Db<-s6$Db
				alpha<-s6$alpha
				beta<-s6$beta
				lambda<-s6$lambda
			}
			if(cc<tol)break
			if(I>max.I){
				warning('Maximum number of iterations exceeded, convergence not obtained.')
				break
			}
			.cc<-c(.cc,cc)
		}
	}
	list(tm=tm,tn=tn,tf=tf,tg=tg,alpha=alpha,beta=beta,lambda=lambda,Da=Da,Db=Db,s.xi=s.xi)
}
dual.bc<-function(y,z,t,subject, knots=NULL, penalties=NULL,df=NULL, k=NULL, control=pfdaControl(),tbase=NULL,subset=NULL){
	{ #Startup
	eval(.X.subset)
	eval(.X.dual.k)
	eval(.X.single.knots)
	eval(.X.dual.penalties)
	}
	if(any(is.na(k))){
		# stop("not finished with number of principal component optimization.")
		if(is.na(k[1])) model.y.single <- single.b(y,     t,subject, knots=knots, penalties=penalties[,1],k=NULL, control=control,subset=subset)
		if(is.na(k[2])) model.z.single <- single.c(z,NULL,t,subject, knots=knots, penalties=penalties[,2],k=NULL, control=control,subset=subset)
		Recall(y,z,t,subject, knots=knots,penalties=penalties,k=c(model.y.single$k,model.z.single$k), control=control)
	}
	else if(any(is.na(penalties))){
		funcall <- match.call()
		eval(.X.optimize.penalties) }
	else { # Pass to core function
		eval(.X.binary.y)
		rtn<-modifyList(if(control$useC){
			eval(.X.read.parameters)
			{ #compute memory requirements
				k=max(ka,kb)
				ni=max(nobs)
				kr=max(control$binary.k0, control$binary.kr)
				dpl = M + sum(nobs^2) + N*k^2 + N*p^2 + p + p*k + k + max(
					M*(3+ka+kb) + (ni*kr+kr+2+5*ka+kb+10*k^2) ,
					#ni*kr + kr + (ni * k + ni + k^2 + k + p + p^2 + k*p + 3*k)+1000,# step W
					M + M*k + N*k + 2* k^2 + 2 * k + k *ni                    ,# step 1/E
					p^2+ M+ M*k                                               ,# step 2
					M + p^2 + p + k^2                                         ,# step 3
					k + 2*k^2 + 2*p^2 + p + p*max(k,8)                        ,# step 4
					2*p^2 + M + p*N + 8*p)                                     # inits
				ipl = 8*p
			}
			{ structure(.C("dual_bc_core",
				y         =as.integer(y), 
				z         =as.double(z),
				B         =Bt,
				tm        =double(p),
				tn         =double(p),
				tf        =matrix(0.0,p,ka),
				tg        =matrix(0.0,p,kb),
				alpha     =matrix(0.0,N,ka),
				beta      =matrix(0.0,N,kb),
				lambda    =matrix(0.0,kb,ka),
				Da        =double(ka),
				Db        =double(kb),
				sxi       =double(1),
				aa        =array(0.0,dim=c(ka,ka,N)),
				ab        =array(0.0,dim=c(ka,kb,N)),
				bb        =array(0.0,dim=c(kb,kb,N)),
				Saa        =array(0.0,dim=c(ka,ka,N)),
				Sab        =array(0.0,dim=c(ka,kb,N)),
				Sbb        =array(0.0,dim=c(kb,kb,N)),
				nobs      =as.integer(nobs),
				N         =N,
				M         =as.integer(M),
				ka        =as.integer(ka),
				kb        =as.integer(kb),
				p         =as.integer(p),
				lm        =penalties[1], 
				ln        =penalties[2], 
				lf        =penalties[3], 
				lg        =penalties[4],
				K         = Kt,
				minV      =control$minimum.variance,
				k0        =control$binary.k0, 
				kr        =control$binary.kr,
				Iterations=control$max.iterations,
				tol       =control$convergence.tolerance,
				dl        =control$C.debug,
				dp        =double(dpl),
				ip        =integer(ipl)
			),class=c('pfda.dual.bc.rawC','pfda.dual.bc','list'))}
		} else {
			structure(dual.cc.core(y,z,B,subject,ka,kb,lm,ln,lf,lg,K,min.v,max.I,tol)
				,class=c('pfda.dual.bc.R','pfda.dual.bc','list'))
		},list(subject=subject))
	}
}
.dual.bc.updatePCS<-function(Yi, Zi, Da, Db, lambda, phii, psii, Y.rhoi, Z.rhoi, alpha, beta){
	ni = NROW(phii)
	# mean functions
	muy = Y.rhoi+phii%*%alpha
	mua = pnorm(muy)
	Y.tilda = muy + (Yi-mua)/dnorm(mua)
	Z.tilda = Zi
  # blocks for means
	A = tcrossprod(Da,phii)                               #' $ A = D_\alpha \phi_i^T $
	B = Da %*% t(lambda) %*% t(psii)                      #' $ B = D_\alpha \Lambda^T \psi_i^T $
	C = lambda %*% A                                      #' $ C = \Lambda D_\alpha \psi_i^T = \Lambda A$
	D = Db %*% t(psii)                                    #' $ D = D_\beta \psi_i^T $
	#blocks associated with variances
	E = diag(1, ni) + (phii %*% A)                        #' $ E = I + \phi_i A = I + \phi_i D_\alpha \phi_i^T $
	F = phii %*% B                                        #' $ F = \phi_i D_\alpha \Lambda^T \psi_i  $ 
	G = psii %*% C                                        #' $ G = \psi_i \Lambda D_\alpha \phi_i $
	H = diag(1, ni) + psii %*% D                          #' $ H = I + \psi_i D_\beta \psi_i $
	# inversion by block
	E1 = solve(E)                                         
	T = solve(H-G%*%E1%*%F)                               #' $ T = (H - G E^{-1} F)^{-1} $
	Q = E1+E1%*%F%*%T%*%G%*%E1                            #' $ Q = E^{-1} + E^{-1} F T G E^{-1} $
	R = -E1%*%F%*%T                                       #' $ R = -E^{-1} F T $
	S = -T%*%G%*%E1                                       #' $ S = -T G E^{-1} $
	# results
	list(
	alpha= (A %*% Q + B %*% S) %*% Y.tilda + (A %*% R + B %*% T) %*% Z.tilda,
	beta = (C %*% Q + D %*% S) %*% Y.tilda + (C %*% R + D %*% T) %*% Z.tilda)
}
.dual.bc.estimatePCS<-function(Y, Z, B, subject, tm, tn, tf, tg, Da, Db, lambda, sxi){
	Y.rho = B%*%tm
	Z.rho = B%*%tn                                                                                  
	
	phi = B%*%tf
	psi = B%*%tg
	ofun<-function(i){
		ix = as.integer(subject)==i
		alpha = rep(0,NCOL(tf))
		beta  = rep(0,NCOL(tg))
		dif = 1
		while(dif>1e-3){
			ab<-.dual.bc.updatePCS(Y[ix], Z[ix], Da, Db, lambda, phi[ix,,drop=F], psi[ix,,drop=F], Y.rho[ix,], Z.rho[ix,], alpha, beta)
			dif<-mean(alpha-ab$alpha)+mean(beta - ab$beta)
			alpha = ab$alpha
			beta  = ab$beta
		}
		list(alpha=alpha,beta=beta)
	}
	ab<-lapply(seq_len(nlevels(subject)), ofun)
	alpha<-matrix(0, nlevels(subject), NCOL(tf))
	beta <-matrix(0, nlevels(subject), NCOL(tg))
	for(i in nlevels(subject)){
		alpha[i,]<-ab[[i]]$alpha
		beta[i,] <-ab[[i]]$beta
	}
	list(alpha=alpha,beta=beta)
}
.dual.bc.n2L<-function(y, z, B, subject, tm, tn, tf, tg, lambda, Da, Db, sxi){
  ka = length(Da)
	kb = length(Db)
	N = nlevels(subject)
	M = length(y)
	phi = B %*% tf
	psi = B %*% tg
	pi  = B %*% tm
	rho = B %*% tn
	ab = .dual.bc.estimatePCS(y, z, B, subject, tm, tn, tf, tg, Da, Db, lambda)
	C = rbind(cbind(diag(Da,ka), diag(Da,ka)%*%t(lambda)),cbind(t(lambda%*%diag(Da,ka)), diag(Db,kb)))
	N*(sum(log(Da))+determinant(diag(Db,kb)-lambda %*% tcrossprod(diag(Da,ka), lambda))$modulus) +	 
	M*log(sxi) +
	with(ab,sum(sapply(seq_len(N),function(i){
		ix = subject==i
		mu_alpha=pnorm( pi[ix]+phi[ix,,drop=F]%*%alpha[i,])
		mu_beta =pnorm(rho[ix]+psi[ix,,drop=F]%*% beta[i,])
		W_alpha.i = (dnorm(qnorm(mu_alpha))^2)/mu_alpha*(1-mu_alpha)
		gamma<-c(alpha[i,],beta[i,])

		# $ \log|I+\Omega W_i \Omega C| $
		A<-diag(1,ka)+crossprod(phi[ix,,drop=F], (W_alpha.i * phi[ix,,drop=F]) %*% diag(Da,ka) )
		determinant(A)$modulus + 
		determinant(diag(1,kb) + crossprod(psi[ix,,drop=F]) - lambda %*% (Da * solve(A)) %*% tcrossprod(diag(Da,ka),lambda) )$modulus + 
		# determinant(diag(1,ka)+ crossprod(phi[ix,],W_alpha.i*phi[ix,,drop=F])%*%C,log=T)$modulus+
		#  $\sumj^{n_i} Y_i $
		sum(-2*(y[ix]*log(mu_alpha) + (1-y[ix])*log(1-mu_alpha)) + (z[ix]-mu_beta)^2) + 
		{ # $ \gamma_i^T C^{-1} \gamma_i $
	  Sib <- solve(diag(Db,kb) - lambda %*% Da %*% lambda)
		Siab <- -crossprod(lambda ,Sib)
		Sia <- diag(1/Da,ka) - Siab %*% lambda %*% diag(1/Da,ka)
		crossprod(gamma,solve(C,gamma))
    crossprod(alpha[i,],Sia %*% alpha[i,]) +
		crossprod(alpha[i,],Siab %*% beta[i,])*2 + 
		crossprod(beta[i,] ,Sib %*% beta[i,])
	  }
		})))
}
AIC.pfda.dual.bc<-function(object,...){
	with(object,.dual.ge.AIC(
		.dual.bc.n2L(y, z, B, subject, tm, tn, tf, tg, lambda, Da, Db, sxi),0,
		B   , B   , B  , B ,
		1   , ka  , 1  , kb,
		lm  , ln  , lf , lg,
		K   , K   , K  , K ,
		1, 1, sxi, sxi)
	)
}
print.pfda.dual.bc<-function(x,...){
	cat('Paired Functional Principal Component Model (Binary/Continuous)\n')
	cat('Formula: ', deparse(attr(x,'formula')),'\n')
	cat(attr(x,'name.y'),' has ', NCOL(x$tf),' principal components\n')
	cat(attr(x,'name.z'),' has ', NCOL(x$tg),' principal components\n')
	cat('penalties are \n');print( penalty.pfda.dual.bc(x))
}
penalty.pfda.dual.bc<-function(object,..)with(object,structure(matrix(c(lm,ln,lf,lg),2,2),dimnames=list(c(attr(object,'name.y'),attr(object,'name.z')),c('mean','pc'))))
}
{ # Additive aka Calcium Model
.dual.ca.i<-function(y,Z,Bt,Bx,subject,kg,kd,min.v){
	tz<- if(!is.null(Z) && NCOL(Z)>=1) solve(crossprod(Z),crossprod(Z,y)) else numeric(0)
	a<-.single.c.i(y,Bt,subject,kg,min.v)
	b<-.single.c.i(y,Bx,subject,kd,min.v)
	lambda <- crossprod(b$alpha,a$alpha)%*%solve(crossprod(a$alpha))
	gd<-array(0,dim=c(kg,kd,nlevels(subject)))

	R=y - Z%*%tz - Bt%*%a$tm - Bx%*%b$tm
	for(i in seq_len(nlevels(subject))){
		ix = i==as.integer(subject)
		R[ix]=R[ix] - Bt[ix,]%*%a$tf%*%a$alpha[i,] - Bx[ix,]%*%b$tf%*%b$alpha[i,]
	}
	sigma  = as.vector(crossprod(R)/length(R))

	for(i in seq_len(nlevels(subject)))gd[,,i]<-tcrossprod(a$alpha[i,],b$alpha[i,])
	{ list(
		tz=tz,
		tt=a$tm,tx=b$tm,
		tf=a$tf,tg=b$tf,
		gamma=a$alpha,delta=b$alpha,
		Dg=a$Da,Dd=b$Da,
		lambda=lambda,
		sigma=sigma,
		gg=array(unlist(a$aa),c(kg,kg,nlevels(subject))),
		gd=gd,
		dd=array(unlist(b$aa),c(kd,kd,nlevels(subject)))
	)}
}
.dual.ca.E.1<-function(Ri,tt,tx,sigma, phi, psi, O){
	# phi = Bti %*% tf
	# psi = Bxi %*% tg
	# O is  list of block inverted Dg,C,Dd from .gen.symblock.solve

	S = .gen.symblock.solve(
		O[[1]]+crossprod(phi)/sigma,
		O[[2]]+crossprod(phi,psi)/sigma,
		O[[3]]+crossprod(psi)/sigma)

	a=crossprod(phi,Ri/sigma)
	b=crossprod(psi,Ri/sigma)

	g=S[[1]]%*%a+S[[2]]%*%b
	d=crossprod(S[[2]],a)+S[[3]]%*%b
	list(
	mu.gamma = as.vector(g),
	mu.delta = as.vector(d),
	Sigma=S,
	gg = tcrossprod(g)+S[[1]],
	gd = tcrossprod(g,d)+S[[2]],
	dd = tcrossprod(d)+S[[3]]
	)
}
.dual.ca.E<-function(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,lambda,Dg,Dd,sigma){
	O<-.gen.symblock.solve(
		diag(Dg,length(Dg)),
		tcrossprod(diag(Dg,length(Dg)),lambda),
		diag(Dd,length(Dd)))
	{
		kg<-NCOL(tf)
		kd<-NCOL(tg)
		n<-nlevels(subject)
		Sgg<-array(0,dim=c(kg,kg,n))
		Sgd<-array(0,dim=c(kg,kd,n))
		Sdd<-array(0,dim=c(kd,kd,n))
		gg <-array(0,dim=c(kg,kg,n))
		gd <-array(0,dim=c(kg,kd,n))
		dd <-array(0,dim=c(kd,kd,n))
		gamma<-matrix(nrow=n,ncol=kg)
		delta<-matrix(nrow=n,ncol=kd)
		R<-y-Z%*%tz-Bt%*%tt-Bx%*%tx
		phi = Bt%*%tf
		psi = Bx%*%tg
	}
	for(i in seq_len(n)){
		ix = i==as.integer(subject)
		r<-.dual.ca.E.1(R[ix],tt,tx,sigma,phi[ix,,drop=F],psi[ix,,drop=F],O)
		{
			Sgg[,,i]<-r$Sigma[[1]]
			Sgd[,,i]<-r$Sigma[[2]]
			Sdd[,,i]<-r$Sigma[[3]]
			gg[,,i]<-r$gg
			gd[,,i]<-r$gd
			dd[,,i]<-r$dd
			gamma[i,]<-r$mu.gamma
			delta[i,]<-r$mu.delta
		}
	}
	list(Sgg=Sgg,Sgd=Sgd,Sdd=Sdd,gamma=gamma,delta=delta,gg=gg,gd=gd,dd=dd)
}
.dual.ca.unpenalized<-function(y,Z,Bt,Bx,subject,tt,tx,tf,tg,gamma,delta,sigma){
	if(is.null(Z)||NCOL(Z)==0)return(list(tz=numeric(0)))
	R=y-Bt%*%tt-Bx%*%tx
	for(i in seq_len(nlevels(subject))){
		ix = i==as.integer(subject)
		R[ix]=R[ix] - Bt[ix,]%*%tf%*%gamma[i,] - Bx[ix,]%*%tg%*%delta[i,]
	}
	list(tz=solve(crossprod(Z),crossprod(Z,R)))
}
.dual.ca.penalized.1<-function(R, B, sigma, l, K){
	solve(crossprod(B)+sigma*l*K,crossprod(B,R))
}
.dual.ca.penalized<-function(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,gamma,delta,sigma,lt,lx,Kt,Kx){
	R=y-Z%*%tz
	for(i in seq_len(nlevels(subject))){
		ix = i==as.integer(subject)
		R[ix]=R[ix] - Bt[ix,]%*%tf%*%gamma[i,] - Bx[ix,]%*%tg%*%delta[i,]
	}

	list( tt= .dual.ca.penalized.1(R-Bx%*%tx,Bt,sigma,lt,Kt),
	      tx= .dual.ca.penalized.1(R-Bt%*%tt,Bx,sigma,lx,Kx))
}
.dual.ca.princcomp<-function(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,gamma,delta,sigma,gg,gd,dd, lf,lg,Kt,Kx){
	R = y - Z%*%tz -Bt%*%tt - Bx%*%tx
	# tf
	for(j in seq_len(NCOL(tf))){
		left=sigma*lf*Kt
		right=numeric(NROW(tf))
		for(i in seq_len(nlevels(subject))){
			ix = i==as.integer(subject)
			left = left + gg[j,j,i]*crossprod(Bt[ix,])
			right= right+ crossprod(
				Bt[ix,,drop=F],
				R[ix]*gamma[i,j]-
					Bx[ix,,drop=F]%*%tg%*%gd[j,,i]-
					(if(NCOL(tf)>1)Bt[ix,,drop=F]%*%tf[,-j]%*%gg[j,-j,i] else 0))
		}
		tf[,j]= .u.orthogonalize(tf[,seq_len(j-1),drop=FALSE],solve(left,right))
	}

	# tg
	for(j in seq_len(NCOL(tg))){
		left=sigma*lg*Kx
		right=numeric(NROW(tg))
		for(i in seq_len(nlevels(subject))){
			ix = i==as.integer(subject)
			left = left + dd[j,j,i]*crossprod(Bx[ix,])
			right= right+ crossprod(
				Bx[ix,,drop=F],
				R[ix]*delta[i,j]-
					Bt[ix,,drop=F]%*%tf%*%gd[,j,i] -
					(if(NCOL(tg)>1)Bx[ix,,drop=F]%*%tg[,-j]%*%dd[j,-j,i] else 0))
		}
		tg[,j]= .u.orthogonalize(tg[,seq_len(j-1),drop=FALSE],solve(left,right))
	}

	# return value
	list(tf = tf, tg = tg)
}
.dual.ca.variances<-function(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,gamma,delta,sigma,gg,gd,dd){
	R=y - Z%*%tz - Bt%*%tt - Bx%*%tx
	sum.gg<-matrix(0,NCOL(tf),NCOL(tf))
	sum.gd<-matrix(0,NCOL(tf),NCOL(tg))
	sum.dd<-matrix(0,NCOL(tg),NCOL(tg))
	for(i in seq_len(nlevels(subject))){
		ix = i==as.integer(subject)
		R[ix]=R[ix] - Bt[ix,]%*%tf%*%gamma[i,] - Bx[ix,]%*%tg%*%delta[i,]
		sum.gg<-sum.gg+gg[,,i]
		sum.gd<-sum.gd+gd[,,i]
		sum.dd<-sum.dd+dd[,,i]
	}
	sigma  = as.vector(crossprod(R)/length(R))

	eg<-eigen(sum.gg/nlevels(subject))
	ed<-eigen(sum.dd/nlevels(subject))

	Tg<-eg$vectors*rep(sign(tf[1,]%*%eg$vectors),each=nrow(eg$vectors))
	Td<-ed$vectors*rep(sign(tg[1,]%*%ed$vectors),each=nrow(ed$vectors))

	lambda = t(solve(sum.gg,sum.gd))

	list(sigma=sigma,
		Dg = eg$values,
		Dd = ed$values,
		lambda = t(Td)%*%lambda%*%Tg,
		tf = tf%*%Tg,
		tg = tg%*%Td,
		gamma = gamma%*%Tg,
		delta = delta%*%Td
	)
}
.dual.ca.core<-function(y,Z,Bt,Bx,subject,kg,kd,lt,lx,lf,lg,Kt,Kx,min.v,max.I,tol){
	{ # initial values
		if(is.null(Z)) Z<-matrix(nrow=NROW(y),ncol=0)
		init<-.dual.ca.i(y,Z,Bt,Bx,subject,kg,kd,min.v)
		init$I<-0
		init$cctrace<-numeric(0)
	}
	within(init,while(I<-I+1){
		E <- .dual.ca.E(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,lambda,Dg,Dd,sigma)
		U <- .dual.ca.unpenalized(y,Z,Bt,Bx,subject,       tt,  tx,tf,tg ,E$gamma,E$delta  ,sigma)
		P <- .dual.ca.penalized(y  ,Z,Bt,Bx,subject,U$tz,  tt,  tx,tf,tg ,E$gamma,E$delta  ,sigma,lt,lx,Kt,Kx)
		C <- .dual.ca.princcomp(y  ,Z,Bt,Bx,subject,U$tz,P$tt,P$tx,tf,tg ,E$gamma,E$delta  ,sigma,E$gg,E$gd,E$dd, lf,lg,Kt,Kx)
		V <- .dual.ca.variances(y ,Z,Bt,Bx,subject,U$tz,P$tt,P$tx,C$tf,C$tg,E$gamma,E$delta,sigma,E$gg,E$gd,E$dd)
		{ #Convergence
			ccl<-numeric(0)
			ccl['tz']<-sum(abs((tz-U$tz)/tz))
			ccl['tt']<-sum(abs((tt-P$tt)/tt))
			ccl['tx']<-sum(abs((tx-P$tx)/tx))
			ccl['tf']<-sum(abs((tf-V$tf)/tf))
			ccl['tg']<-sum(abs((tg-V$tg)/tg))
			ccl['Dg']<-sum(abs((Dg-V$Dg)/Dg))
			ccl['Dd']<-sum(abs((Dd-V$Dd)/Dd))
			(cc<-sum(ccl,na.rm=TRUE))
			{ # reassign values
				gg    <-E$gg
				gd    <-E$gd
				dd    <-E$dd
				tz    <-U$tz
				tt    <-P$tt
				tx    <-P$tx
				tf    <-V$tf
				tg    <-V$tg
				gamma <-V$gamma
				delta <-V$delta
				sigma <-V$sigma
				Dg    <-V$Dg
				Dd    <-V$Dd
				lambda<-V$lambda
			}
			if(cc<tol)break
			if(I>max.I){
				warning('Maximum number of iterations exceeded, convergence not obtained.')
				break
			}
			cctrace<-c(cctrace,cc)
		}
	})
}
.dual.ca.n2L<-function(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,lambda,Dg,Dd,sigma,...){
	ll<-0
	phi = Bt%*%tf
	psi = Bx%*%tg
	Dg<-diag(Dg,length(Dg))
	Dd<-diag(Dd,length(Dd))
	Ry<-y-Bt%*%tt-Bx%*%tx-(if(is.null(Z))0 else Z%*%tz)
	for(i in nlevels(subject)){
		ix = i==as.integer(subject)
		U  = phi[ix,,drop=F]%*%tcrossprod(Dg,psi[ix,]%*%lambda)
		S  = phi[ix,,drop=F]%*%tcrossprod(Dg,phi[ix,,drop=F])+ psi[ix,,drop=F]%*%tcrossprod(Dd,psi[ix,,drop=F]) + U + t(U) + diag(sigma,sum(ix))
		ll = ll+determinant(S)$modulus+crossprod(Ry[ix],solve(S)%*%Ry[ix])
	}
	attr(ll,"log")<-TRUE
	ll
}
AIC.pfda.additive<-function(object,...){
	with(object,.dual.ge.AIC(.dual.ca.n2L(y,Z,Bt,Bx,subject, tz, tt, tx, tf, tg, lambda, Dg, Dd, sigma),
		kz,
		Bt, Bt, Bx, Bx,
		1 , kg, 1 , kd,
		lt, lf, lx, lg,
		Kt, Kt, Kx, Kx,
		sigma, sigma, sigma, sigma))
 # with(object,.dual.ca.AIC.rawC(object, y, Z, Bt, Bx, subject, Kt, Kx))
}
logLik.pfda.additive<-function(object,...,newdata=NULL,n2L=TRUE){
	with(object,with(newdata,.dual.ca.n2L(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,lambda,Dg,Dd,sigma,...)))
}
dual.ca<-function(y,Z,t,x,subject,knots=NULL,penalties=NULL,df=NULL,k=NULL,bases=NULL,control=pfdaControl()){
	# name.t = deparse(substitute(t))
	# name.x = deparse(substitute(x))
	{ #setup material
	Z = if(is.null(Z))matrix(nrow=length(y),ncol=0) else as.matrix(Z)
	if(is.null(bases))bases=new.env()
	stopifnot(is.environment(bases))
	if(!(exists("tbase",envir=bases,inherits=FALSE) && exists("xbase",envir=bases, inherits=FALSE))){ # knot and base identification
		if(is.null(knots)){
			if(length(control$nknots)==2) { nkt <- control$nknots[1]; nkx <- control$nknots[2] } else nkt <- nkx <- control$nknots[1]
			kt<-expand.knots(unique(quantile(t,seq(0,1,length.out=nkt))))
			kx<-expand.knots(unique(quantile(x,seq(0,1,length.out=nkt))))
		} else if(length(knots)==2) {
			if(name.t %in% names(knots)) kt <-knots[[name.t]]
			else if('t' %in% names(knots)) kt<-knots[['t']]
			else kt <- knots[[1]]
			if(name.x %in% names(knots)) kx <-knots[[name.x]]
			else if('x' %in% names(knots)) kx<-knots[['x']]
			else kx <- knots[[2]]
		} else stop("knots needs to be a list of length 2")
		assign('tbase',OBasis(kt),envir=bases)
		assign('xbase',OBasis(kx),envir=bases)
		knots<-list(kt,kx)
		# names(knots)<-c(name.t,name.x)
	}
	{ # evaluate Bases
		Bt = evaluate(bases$tbase,t)
		Bx = evaluate(bases$xbase,x)
		Kt = OuterProdSecondDerivative(bases$tbase)
		Kx = OuterProdSecondDerivative(bases$xbase)
		Bx<-Bx[,-1]  #Removing the first column of Bx is to improve stability and remove interdependence between t and x
		Kx<-Kx[-1,-1]
	}
	eval(.X.dual.k)
	eval(.X.dual.penalties)
	}
	if(any(is.na(k))){
		# stop("identification of number of principle components is not done yet.")
		k<-c(1,1)
		model0<-RecallWith(k=k)
		AIC0<-AIC(model0)
		while(TRUE){
			modelA<-try(RecallWith(k=c(model0$kg+1,model0$kd)),silent=TRUE)
			modelB<-try(RecallWith(k=c(model0$kg,model0$kd+1)),silent=TRUE)
			AICA<- if(is(modelA,'pfda.additive')) AIC(modelA) else Inf
			AICB<- if(is(modelB,'pfda.additive')) AIC(modelB) else Inf
			if(AIC0<AICA && AIC0<AICB) break
			else if(AICA<AICB){
				model0<-modelA
				AIC0<-AICA
			} else {
				model0<-modelB
				AIC0<-AICB
			}
		}
		return(model0)
	} else
	if (any(is.na(penalties))) {
		funcall <- structure(match.call(),envir=parent.frame(),expand.dots = FALSE)
		if(is.null(control$optim.start)){
			lt <- l.from.df(2.1,Bt,Kt)
			lx <- l.from.df(2.1,Bx,Kx)
			control$optim.start<-c(lt,lx,lt,lx)
		}
		eval(.X.optimize.penalties) }
	else {
		rtn<-if(control$useC){
			{ # compute memory and parameters
				kz  = if(is.null(Z))0L else NCOL(Z)
				kg  = as.integer(k[1])
				kd  = as.integer(k[2])
				N   = nlevels(subject)
				nobs= table(subject)
				M   = length(y)
				pt  = ncol(Bt)
				px  = ncol(Bx)
				p   = max(pt, px)
				k   = max(kg, kd)
				dpl = kz + kg + kd + pt +px + pt*kg + px*kd + kg*kd + max(
					pt*pt*N + px*px*N + 2*p^2 + M + M*k+ p*N + 8*p,
					2^kg^2 + 2*kd^2 + kg*kd+ max(kg,kd)*10 + 2*kg*kd+(1+kg+kd)*M,
					M + kz^2 + 10*kz + M*max(kg, kd),
					10*p + p^2 +M,
					M + p^2 + p + max( 10*pt ,  M + N*(kg+kd) + M*max(kg,kd) ),
					2*kg^2 + 2*kd^2 + M + 10*k + 9*p + 3*p^2 + N*k + M*k
					)
			}
			structure(.C('dual_ca_core', residuals=y, Z=Z, Bt=Bt, Bx=Bx,
				tz=double(kz), tt=double(pt), tx=double(px),
				tf=matrix(0,pt,kg), tg=matrix(0,px,kd),
				gamma=matrix(0,N,kg), delta=matrix(0,N,kd),
				lambda=matrix(0,kd,kg), Dg=double(kg), Dd=double(kd),
				sigma=0.0,
				gg=array(0,c(kg,kg,N)),   gd=array(0,c(kg,kd,N)),  dd=array(0,c(kd,kd,N)),
				Sgg=array(0,c(kg,kg,N)), Sgd=array(0,c(kg,kd,N)), Sdd=array(0,c(kd,kd,N)),
				nobs=nobs, N=N, M=M, kz=kz, kg=kg, kd=kd, pt=pt, px=px,
				lt=penalties[1], lx=penalties[2], lf=penalties[3], lg=penalties[4], Kt=Kt, Kx=Kx,
				minV=control$minimum.variance, Iterations=control$max.iterations, tol=control$convergence.tolerance,
				dl=control$C.debug, dp=double(dpl), ip=integer(max(6*p,kz)))
			,class=c('pfda.additive.rawC','pfda.additive','list'))
		} else {
			structure(.dual.ca.core(y,Z,Bt,Bx,subject,k[1],k[2],penalties[1],penalties[2],penalties[3],penalties[4],Kt,Kx,control$minimum.variance,control$max.iterations,control$convergence.tolerance)
			,class=c('pfda.additive.R','pfda.additive','list'))
		}
		rtn$tbase<-bases$tbase
		rtn$xbase<-bases$xbase
		rtn$y<-y
		rtn$subject<-subject
		rtn
	}
}
plot.pfda.additive<-function(x,...){
	with(x,{
		layout(matrix(1:4,nrow=2,ncol=2,byrow=T))
		plot(tbase,tt, main=paste("Plot of mean curve for",attr(x,'name.t')),xlab=attr(x,'name.t'),ylab=attr(x,'name.y'))
		plot(tbase,tf, main=paste("Principle components for",attr(x,'name.t')),xlab=attr(x,'name.t'),ylab='')
		plot(xbase,c(0,tx), main=paste("Plot of mean curve for ",attr(x,'name.x')),xlab=attr(x,'name.x'),ylab=attr(x,'name.y'))
		plot(xbase,rbind(0,tg), main=paste("Principle components for ",attr(x,'name.x')),xlab=attr(x,'name.x'),ylab='')
	})
}
print.pfda.additive<-function(x,...){
	cat('Additive Principal Component Model\n')
	cat('Formula: ', deparse(attr(x,'formula')),'\n')
	cat(attr(x,'name.t'),' has ', NCOL(x$tf),' principal components\n')
	cat(attr(x,'name.x'),' has ', NCOL(x$tg),' principal components\n')
	cat('penalties are \n');print( penalty.pfda.additive(x))
}
penalty.pfda.additive<-function(object,..)with(object,structure(matrix(c(lt,lx,lf,lg),2,2),dimnames=list(c(attr(object,'name.t'),attr(object,'name.x')),c('mean','pc'))))
persp.pfda.additive<-function(x,col.fun=NULL,nt=11,nx=11,...){
	if(is.null(col.fun))col.fun<-heat.colors
	with(x,{
	tx.grid<-expand.grid(
		t=seq(tbase@knots[4],tbase@knots[length(tbase@knots)-3],length=nt),
		x=seq(xbase@knots[4],xbase@knots[length(xbase@knots)-3],length=nx))
	Btg<-evaluate(tbase,tx.grid$t)
	Bxg<-evaluate(xbase,tx.grid$x)[,-1]
	yg<-matrix(Btg%*%(tt+tf%*%apply(gamma,2,mean)) + Bxg%*%(tx+tg%*%apply(delta,2,mean)),nx,nt,byrow=T)
	{
		# this code was copied and modified from the persp example page
		z <- yg
		nrz <- nrow(yg)
		ncz <- ncol(yg)
		# Create a function interpolating colors in the range of specified colors
		nbcol <- 200
		color <- col.fun(nbcol)
		zfacet <- yg[-1, -1] + yg[-1, -ncz] + yg[-nrz, -1] + yg[-nrz, -ncz]
		facetcol <- cut(zfacet, nbcol)
	}
	persp(yg,zlab=attr(x,"name.y"),xlab=attr(x,"name.x"),ylab=attr(x,"name.t"), 
		col=color[facetcol], ...)
	})
}
}
{ # general
	.dual.ge.AIC<-function(n2L,kz,
		B1, B2, B3, B4,
		k1, k2, k3, k4,
		l1, l2, l3, l4,
		K1, K2, K3, K4,
		s1, s2, s3, s4){
		as.vector(n2L)+2*(kz +
			k1*.pfda.df(B1,l1,K1,s1) +
			k2*.pfda.df(B2,l2,K2,s2) +
			k3*.pfda.df(B3,l3,K3,s3) +
			k4*.pfda.df(B4,l4,K4,s4)
		)
	}
	pfda<-function(model, data=environment(model), ..., driver){
		mf <- pfdaParseFormula(model,data)
		if(missing(driver))driver = infer.driver(mf)
		structure(with(mf,switch(driver,
			single.continuous = structure(single.c(response,additive,splinegroup[[1]],splinegroup[[2]],...),name.y=attr(response,'name'),name.t=attr(splinegroup[[1]],'name')),
			single.binary = structure(single.b(response,splinegroup[[1]],splinegroup[[2]],...),name.y=attr(response,'name'),name.t=attr(splinegroup[[1]],'name')),
			dual.continuous =  structure(dual.cc(response[[1]],response[[2]],splinegroup[[1]],splinegroup[[2]],...),name.y=attr(response[[1]],'name'),name.z=attr(response[[2]],'name'),name.t=attr(splinegroup[[1]],'name')),
			dual.mixed = structure(dual.bc(response[[1]],response[[2]],splinegroup[[1]],splinegroup[[2]],...),name.y=attr(response[[1]],'name'),name.z=attr(response[[2]],'name'),name.t=attr(splinegroup[[1]],'name'),name.x=attr(splinegroup[[1]],'name')),
			additive = structure(dual.ca(response,additive,splinegroup[[1]],splinegroup[[2]],splinegroup[[3]],...),name.y=attr(response,'name'),name.t=attr(splinegroup[[1]],'name'),name.x=attr(splinegroup[[2]],'name'))
		)), formula = model)
	}
}
