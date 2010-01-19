IT_bin_s_1<-function(){ for(attempt in seq(100)){
	# totally random input
	N<-trunc(runif(1,10,50))
	nobs<-trunc(runif(N,3,8))
	id<-factor(rep(seq(N),nobs))
	M<-sum(nobs)
	y<-sample(c(T,F),M,replace=T)
	x<-rnorm(M)
	
	lm=lf=100
	.pfda.debug.w <- NULL
	.pfda.debug.ww <- NULL
	.pfda.debug.aa <- NULL
	
	suppressWarnings(capture.output(c0<-pfda(y~x|id,k=1,penalties=c(100,100),response.class='binary',
		control=pfdaControl(Cdebug=pfda:::Cdebug(426),debug.capture.w=TRUE,trace=FALSE))))
	w<-.pfda.debug.w
	ww<-.pfda.debug.ww
	(tm<-c0@parameters@theta_mu)
	(tf<-c0@parameters@Theta_f)
	(Alpha<-c0@parameters@Alpha)
	(Da<-c0@parameters@Da)
	Sigma_aa<-c0@parameters@Sigma
	obase<-c0@Basis
	B<-evaluate(obase,x)
	
	r1<-.single.b.1(B,id,w,ww,tm,tf,Da)
	
	if(!isTRUE(all.equal(c0@parameters@Alpha,r1$alpha)))return("failed on alpha")
	if(!isTRUE(all.equal(unlist(r1$Sa), as.numeric(c0@parameters@Sigma))))return("failed on Sigma_aa")
	if(!isTRUE(all.equal(unlist(r1$aa),as.numeric(.pfda.debug.aa))))return("failed on aa_hat")
	TRUE
} }
IT_bin_s_2<-function(){ for(attempt in seq(100)){
	# totally random input
	N<-trunc(runif(1,10,50))
	nobs<-trunc(runif(N,3,8))
	id<-factor(rep(seq(N),nobs))
	M<-sum(nobs)
	y<-sample(c(T,F),M,replace=T)
	x<-rnorm(M)

	lm=lf=100

	.pfda.debug.w <- NULL
	.pfda.debug.ww <- NULL
	.pfda.debug.aa <- NULL

	suppressWarnings(invisible(capture.output(
	c0<-pfda(y~x|id,k=1,penalties=c(100,100),response.class='binary',
		control=pfdaControl(Cdebug=pfda:::Cdebug(427),debug.capture.w=TRUE,trace=F))
	)))
	w<-.pfda.debug.w
	ww<-.pfda.debug.ww
	(tm<-c0@parameters@theta_mu)
	(tf<-c0@parameters@Theta_f)
	(Alpha<-c0@parameters@Alpha)
	(Da<-c0@parameters@Da)
	Sigma_aa<-c0@parameters@Sigma
	obase<-c0@Basis
	B<-evaluate(obase,x)
	K<-OuterProdSecondDerivative(obase)

	# debug(.bin_s_2)
	r1<-.single.b.2(B,id,w,tf,Alpha,lm,K)
		
	if(!isTRUE(all.equal(tm,as.double(r1))))return("failed to find tm correctly")
	TRUE
} }
IT_bin_s_3<-function(){ for(attempt in seq(100)){
	# totally random input
	N<-trunc(runif(1,10,50))
	nobs<-trunc(runif(N,3,8))
	id<-factor(rep(seq(N),nobs))
	M<-sum(nobs)
	y<-sample(c(T,F),M,replace=T)
	x<-rnorm(M)
	k<-trunc(runif(1,1,4))

	lm=lf=100

	.pfda.debug.w <- NULL
	.pfda.debug.ww <- NULL
	.pfda.debug.aa <- NULL

	C.seed<-.Random.seed
	suppressWarnings(invisible(capture.output(
	c0<-pfda(y~x|id,k=k,penalties=c(100,100),response.class='binary',
		control=pfdaControl(Cdebug=pfda:::Cdebug(427),debug.capture.w=TRUE,trace=F))
	)))
	w<-.pfda.debug.w
	ww<-.pfda.debug.ww
	.aa.raw<-.pfda.debug.aa
	dim(.aa.raw)<-c(k,k,N)
	aa<-vector('list',N);for(i in seq(N))aa[[i]]<-matrix(.aa.raw[,,i],k,k)

	(tm<-c0@parameters@theta_mu)
	(tf<-c0@parameters@Theta_f)
	(Alpha<-c0@parameters@Alpha)
	(Da<-c0@parameters@Da)
	Sigma_aa<-c0@parameters@Sigma
	obase<-c0@Basis
	B<-evaluate(obase,x)
	K<-OuterProdSecondDerivative(obase)

	# debug(.bin_s_2)
	r1<-.single.b.3(B,id,w,tm,tf,Alpha,aa,lf,K)

	C.seed->.Random.seed
	# suppressWarnings(invisible(capture.output(
	c1<-pfda(y~x|id,k=k,penalties=c(100,100),response.class='binary',
		control=pfdaControl(Cdebug=pfda:::Cdebug(428, 134:137),debug.capture.w=TRUE,trace=F))
	# )))
	(tf2<-c1@parameters@Theta_f*rep(sign(c1@parameters@Theta_f[1,]),each=NROW(tf)))
	
	if(!isTRUE(all.equal(tf2,r1)))return("failed to find tf correctly")
	TRUE
} }
IT_bin_s_4<-function(){ for(attempt in seq(100)){
	# totally random input
	N<-trunc(runif(1,10,50))
	nobs<-trunc(runif(N,3,8))
	id<-factor(rep(seq(N),nobs))
	M<-sum(nobs)
	y<-sample(c(T,F),M,replace=T)
	x<-rnorm(M)
	k<-trunc(runif(1,1,4))

	lm=lf=100

	.pfda.debug.w <- NULL
	.pfda.debug.ww <- NULL
	.pfda.debug.aa <- NULL
	
	C.seed<-.Random.seed
	suppressWarnings(invisible(capture.output(
	c0<-pfda(y~x|id,k=k,penalties=c(100,100),response.class='binary',
		control=pfdaControl(Cdebug=pfda:::Cdebug(428),debug.capture.w=TRUE,trace=F))
	)))
	w<-.pfda.debug.w
	ww<-.pfda.debug.ww
	.aa.raw<-.pfda.debug.aa
	dim(.aa.raw)<-c(k,k,N)
	aa<-vector('list',N);for(i in seq(N))aa[[i]]<-matrix(.aa.raw[,,i],k,k)

	(tm<-c0@parameters@theta_mu)
	(tf<-c0@parameters@Theta_f)
	(Alpha<-c0@parameters@Alpha)
	(Da<-c0@parameters@Da)
	Sigma_aa<-c0@parameters@Sigma
	obase<-c0@Basis
	B<-evaluate(obase,x)
	K<-OuterProdSecondDerivative(obase)

	# debug(.bin_s_2)
	r1<-.single.b.4(tf,Alpha,aa)

	C.seed->.Random.seed
	# suppressWarnings(invisible(capture.output(
	c1<-pfda(y~x|id,k=k,penalties=c(100,100),response.class='binary',
		control=pfdaControl(Cdebug=pfda:::Cdebug(429),debug.capture.w=TRUE,trace=F))
	# )))
	tf2<-c1@parameters@Theta_f
	alpha<-c1@parameters@Alpha
	Da2<-c1@parameters@Da

	if(!isTRUE(all.equal(tf2,r1$tf)))return("tf does not match")
	if(!isTRUE(all.equal(Da2,r1$Da)))return("Da does not match")
	if(!isTRUE(all.equal(alpha,r1$alpha)))return("alpha does not match")
	TRUE
} }
IT_bin_s_w<-function(){
	# totally random input
	N<-trunc(runif(1,10,50))
	nobs<-trunc(runif(N,3,8))
	id<-factor(rep(seq(N),nobs))
	M<-sum(nobs)
	y<-sample(c(T,F),M,replace=T)
	x<-rnorm(M)
	k<-trunc(runif(1,1,4))

	lm=lf=100
	
	.pfda.debug.w <- NULL
	.pfda.debug.ww <- NULL
	.pfda.debug.aa <- NULL

	C.seed<-.Random.seed
	suppressWarnings(invisible(capture.output(
	c0<-pfda(y~x|id,k=k,penalties=c(100,100),response.class='binary',
		control=pfdaControl(Cdebug=pfda:::Cdebug(429),debug.capture.w=TRUE,trace=F))
	)))
	w<-.pfda.debug.w
	ww<-.pfda.debug.ww
	.aa.raw<-.pfda.debug.aa
	dim(.aa.raw)<-c(k,k,N)
	aa<-vector('list',N);for(i in seq(N))aa[[i]]<-matrix(.aa.raw[,,i],k,k)

	(tm<-c0@parameters@theta_mu)
	(tf<-c0@parameters@Theta_f)
	(Alpha<-c0@parameters@Alpha)
	(Da<-c0@parameters@Da)
	Sigma_aa<-c0@parameters@Sigma
	obase<-c0@Basis
	B<-evaluate(obase,x)
	K<-OuterProdSecondDerivative(obase)

	# debug(.bin_s_2)
	r1<-.single.b.4(tf,Alpha,aa)

	C.seed->.Random.seed
	# suppressWarnings(invisible(capture.output(
	c1<-pfda(y~x|id,k=k,penalties=c(100,100),response.class='binary',
		control=pfdaControl(Cdebug=pfda:::Cdebug(430),debug.capture.w=TRUE,trace=F))
	# )))
	tf2<-c1@parameters@Theta_f
	alpha<-c1@parameters@Alpha
	Da2<-c1@parameters@Da

	if(!isTRUE(all.equal(tf2,r1$tf)))return("tf does not match")
	if(!isTRUE(all.equal(Da2,r1$Da)))return("Da does not match")
	if(!isTRUE(all.equal(alpha,r1$alpha)))return("alpha does not match")
	TRUE




}
