# unittests.R

#  	library(pfda);i=1
#		set.seed(123)
#	  with(getNamespace('pfda'),with(new.env(),browser()))

{ # Simulate data Expressions
sim.tm<-function(B,K,pen=1){
	tmp<-matrix(rnorm(NROW(K)*1000),NROW(K),1000)
	p<-apply(tmp,2,function(x)crossprod(x,K%*%x))
	p<-p/sum(p)
	i=sample(1:1000,1,prob=p)
	
	i=which.min(p)
	t1<-tmp[,i]
	
	solve(crossprod(B)+pen*K,crossprod(B,B%*%t1))
}
sim.tf<-function(k,B,K,pen=1){
	tf<-replicate(k,sim.tm(B,K,pen))
	tf%*%solve(chol(crossprod(tf)))
}
sim.single.parameters<-expression({
	ka   = sample(1:3,1)
	N    = sample(7:25,1)*ka*ka
	nobs = sample(3:10,N,replace=T)
	M    = sum(nobs)
	times= runif(M,0,100)
	knots= expand.knots(seq(0,100,l=21))
	obase= OBasis(knots=knots)
	B    = evaluate(obase,times)
	K    = OuterProdSecondDerivative(obase)
	p    = NCOL(B)
	
	sigma= rchisq(1,1)
	
	tm = sim.tm(B,K)
	tf = sim.tf(ka,B,K)
	
	Da    <- sort(rchisq(ka,ka),T)
	
	alpha = matrix(rnorm(ka*N, sd = rep(sqrt(Da),each=N)),N,ka)
})
sim.single.c<-expression({
		eval(get('sim.single.parameters',envir=parent.env(environment())))
	subject = factor(rep(seq(N),nobs) )
	y = numeric(M)
	for(i in seq(N)){
		ix    = ( i == as.integer(subject) )
		y[ix] = B[ix,]%*%(tm +tf%*%alpha[i,])
	}
	y<-y+rnorm(M,0,sigma)
})
sim.dual.parameters<-expression({
	ka   = sample(1:3,1)
	kb   = sample(1:3,1)
	N    = sample(7:25,1)*ka*kb
	nobs = sample(3:10,N,replace=T)
	M    = sum(nobs)
	times= runif(M,0,100)
	knots= expand.knots(seq(0,100,l=21))
	obase= OBasis(knots=knots)
	B    = evaluate(obase,times)
	K    = OuterProdSecondDerivative(obase)
	p    = NCOL(B)
	
	s.eps= rchisq(1,1)
	s.xi = rchisq(1,1)
	
	tm = sim.tm(B,K)
	tn = sim.tm(B,K)
	tf = sim.tf(ka,B,K)
	tg = sim.tf(kb,B,K)
	
	while(T){ #break at end
		Da    <- sort(rchisq(ka,ka),T)
		Db    <- sort(rchisq(kb,kb),T)
		lambda<- matrix(rnorm(ka*kb),kb,ka)
		S     <-rbind(cbind(diag(Da,ka),t(lambda)),cbind(lambda,diag(Db,kb)))
	if(all(eigen(S)$values>0))break # check for positive definite
	}
	
	ab    = mvrnorm(N,mu=rep(0,ka+kb),Sigma=S)
	alpha = ab[,1:ka,drop=FALSE]
	beta  = ab[,-(1:ka),drop=FALSE]
})
sim.dual.cc<-expression({
	eval(get('sim.dual.parameters',envir=parent.env(environment())))
	subject = factor(rep(seq(N),nobs))
	w = numeric(M)
	z = numeric(M)
	for(i in seq(N)){
		ix    = (i == as.integer(subject))
		w[ix] = B[ix,]%*%(tm +tf%*%alpha[i,])
		z[ix] = B[ix,]%*%(tn +tg%*%beta[i,])
	}
	y<-w+rnorm(M,0,s.eps)
	z<-z+rnorm(M,0,s.xi)
})
sim.dual.bc<-expression({ 
	eval(get('sim.dual.cc',envir=parent.env(environment())))
	y<-w<mean(w)
	
	ww = vector('list',N)
	btb=array(0,c(p,p,N))
	for(i in seq(N)){
		ix      = i==as.integer(subject)
		ww[[i]] = tcrossprod(w[ix])
		btb[,,i]= crossprod(B[ix,])
	}
})
sim.dual.ca<-expression({
	kg   = sample(1:3,1)
	kd   = sample(1:3,1)
	n=N  = sample(7:25,1)*kg*kd
	nobs = sample(3:10,n,replace=T)
	M    = sum(nobs)
	times= runif(M,0,100)
	tk   = expand.knots(seq(0,100,l=21))
	tbase= OBasis(knots=tk)
	Bt   = evaluate(tbase,times)
	Kt   = OuterProdSecondDerivative(tbase)
	X    = runif(M,0,1)
	xk   = expand.knots(seq(0,1,l=11))
	xbase= OBasis(knots=xk)
	Bx   = evaluate(xbase,X)
	Kx   = OuterProdSecondDerivative(xbase)
	
	sigma= rchisq(1,1)
	
	tt = sim.tm(Bt,Kt)
	tx = sim.tm(Bx,Kx)
	tf = positive.first.row(sim.tf(kg,Bt,Kt))
	tg = positive.first.row(sim.tf(kd,Bx,Kx))
	
	while(T){ #break at end
		Dg    <- sort(rchisq(kg,kg),T)
		Dd    <- sort(rchisq(kd,kd),T)
		lambda<- matrix(rnorm(kg*kd),kd,kg)
		S     <-rbind(cbind(diag(Dg,kg),tcrossprod(diag(Dg,kg),lambda)),cbind(lambda%*%diag(Dg,kg),diag(Dd,kd)))
	if(all(eigen(S)$values>0))break # check for positive definite
	}
	
	gd    = mvrnorm(N,mu=rep(0,kg+kd),Sigma=S)
	gamma = gd[,1:kg,drop=FALSE]
	delta = gd[,-(1:kg),drop=FALSE]

	kz<-sample(1:3,1)
	Z<-matrix(rnorm(M*kz),M,kz)
	tz<-rnorm(kz)
	
	subject = factor(rep(seq(n),nobs))
	w = numeric(M)
	z = numeric(M)
	for(i in seq(N)){
		ix    = (i == as.integer(subject))
		w[ix] = Bt[ix,]%*%(tt +tf%*%gamma[i,]) + Bx[ix,]%*%(tx +tg%*%delta[i,])
	}
	y<-w+Z%*%tz +rnorm(M,0,sigma)
})
sim.aa<-function(N,ka) array(replicate(N,crossprod(matrix(rnorm(ka*ka),ka,ka))),c(ka,ka,N))
sim.ab<-function(N,ka,kb) array(rnorm(ka*kb*N),c(ka,kb,N))
sim.dual.aa_ab_bb<-expression({
aa<-sim.aa(N,ka)
ab<-sim.ab(N,ka,kb)
bb<-sim.aa(N,kb)
})
def.defaults<-expression({
	min.v<-1e-4
	tol<-1e-2
	max.I=1e5
	lt<-lx<-lf<-lg<-1
})
UT_generate<-function(X){
envir = environment()
function(n=100){
if(class(try(replicate(n,eval(X,envir=envir))))=="try-error"){
	cat("FAILED\n")
	return(invisible(FALSE))
}
invisible(TRUE)
}
}
test.equal<-function(name,envir=parent.frame()){ eval(substitute(
stopifnot(all.equal(as.vector(C$name),as.vector(R$name)))
),envir=envir)
}
test.equal.sa<-function(name, envir=parent.frame()){ eval(substitute({
u<-rep(upper.tri(R$name[,,1],TRUE),length.out=prod(dim(R$name)))
for(i in seq(N))stopifnot(all.equal(as.vector(C$name)[u],as.vector(R$name)[u]))
}),envir=envir)
}

}
{ # Utilities
UT_pfdaAlloc<-function(){
	cat('UNIT TEST - pfdaAlloc_i ... ')
	if(.C("test_pfdaAlloc_int",integer(1),0L)[[1]])cat('PASS\n') else cat('FAIL\n')
	cat('UNIT TEST - pfdaAlloc_d ... ')
	if(.C("test_pfdaAlloc_double",integer(1),0L)[[1]])cat('PASS\n') else cat('FAIL\n')
}
UT_pfda_transpose<-function(){
  cat('UNIT TEST - pfda_transpose\n')
  if(is.loaded('test_pfda_transpose')){
		cat("testing with random input...")
		n<-100
		a<-trunc(runif(n,max=7))+2
		b<-trunc(runif(n,max=7))+2
		for(i in seq(n)){	
			X <- matrix(runif(a[i]*b[i]),nrow=a[i],ncol=b[i])
			c_result<-{ .C("test_pfda_transpose",
				X=as.double(X),as.integer(a[i]),as.integer(b[i]),
				0L,
				# c_debug(310),
				double(a[i]*b[i]))}
			if(any(t(X)!=c_result$X)){
				cat('Failed on attempt ',i,'\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no transpose function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_transpose<-function(){
  cat('UNIT TEST - transpose\n')
  if(is.loaded('test_transpose')){
		cat("testing with random input...")
		n<-100
		a<-trunc(runif(n,max=7))+2
		b<-trunc(runif(n,max=7))+2
		for(i in seq(n)){	
			X <- matrix(runif(a[i]*b[i]),nrow=a[i],ncol=b[i])
			c_result<-.C("test_transpose",X=as.double(X),as.integer(a[i]),as.integer(b[i]))
			if(any(t(X)!=c_result$X)){
				cat('Failed on attempt ',i,'\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no transpose function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_checkdebug<-function(){
  cat('UNIT TEST - checkdebug\n')
  if(is.loaded('test_checkdebug')&&is.loaded("checkdebug")){
		{	cat("testing with sequence for individual to 1024...")
			n<-1024
			for(i in seq(n)){	
				debuglevel <- rep(0,i)
				debuglevel[i]<-TRUE
				c_result<-.C("test_checkdebug",as.integer(c(length(debuglevel),debuglevel)),as.integer(i),rtn=integer(1))
				if(!c_result$rtn){
					cat('Failed on ',i,'\n')
					return(invisible(FALSE))
				}
			}
			cat("PASS\n")
		}
		{ cat("testing with random input for debug indicators...")
			debuglevel <- rep(0,n)
			j<-sort(unique(trunc(runif(round(n/2),max=n))))
			debuglevel[j]<-TRUE
			for(i in seq(n)){
				c_result<-.C("test_checkdebug",as.integer(c(length(debuglevel),debuglevel)),as.integer(i),rtn=integer(1))
				if(c_result$rtn!=(i %in% j)){
					cat('Failed on checking',i,'\n')
					return(invisible(FALSE))
				}
			}
			cat("PASS\n")
		}
		{ cat("testing with deliberate small debug list...")
			debuglevel <- rep(0,round(n/2))
			j<-seq(n/10)*5
			debuglevel[j]<-TRUE
			for(i in seq(n)){
				c_result<-.C("test_checkdebug",as.integer(c(length(debuglevel),debuglevel)),as.integer(i),rtn=integer(1))
				if(c_result$rtn!=(i %in% j)){
					cat('Failed on checking',i,'\n')
					return(invisible(FALSE))
				}
			}
			cat("PASS\n")
		}
    { cat("testing for no debug...")
			for(i in seq(n)){
				c_result<-.C("test_checkdebug",0L,as.integer(i),rtn=1L)
				if(c_result$rtn){
					cat('Failed on checking',i,'\n')
					return(invisible(FALSE))
				}
			}
			cat("PASS\n")
		} 
		invisible(TRUE)
	} else {
    cat("Failed, no checkdebug function found.\n")
		invisible(FALSE)
  }
}
UT_pfda_cond_dd<-function(){
	cat('UNIT TEST - pfda_cond_dd\n')
	n<-100
	cat('testing for k<=10, N=',n,'...')
	for(k in seq(10)){
		delta<-matrix(rnorm(k*n),nrow=n,ncol=k)
		S_dd<-array(rnorm(n*k**2,sd=3),dim=c(k,k,n))
		for(i in seq(n)){
			S_dd[,,i]<-crossprod(S_dd[,,i])
			r_result<-tcrossprod(delta[i,])+S_dd[,,i]
			c_result<-.C("test_pfda_cond_dd",dd_hat=double(k**2),as.integer(i-1),as.integer(n),as.double(delta),as.double(S_dd), as.integer(k), 0L)
			eq_data<-all.equal(as.double(r_result),c_result$dd_hat)
			if(!isTRUE(eq_data)){
				cat("failed on k=",k," i=",i,"\n")
				return(invisible(FALSE))
			}
		}
	}
	cat("PASS\n")
	return(invisible(TRUE))
}
UT_pfda_sum_cond_dd<-function(){
	cat('UNIT TEST - pfda_sum_cond_dd\n')
	n<-100
	cat('testing for k<=20, N=',n,'...')
	for(k in seq(20)){
		delta<-matrix(rnorm(k*n),nrow=n,ncol=k)
		S_dd<-array(rnorm(n*k**2,sd=3),dim=c(k,k,n))
		Sum_d<-matrix(0,k,k)
		for(i in seq(n)){ Sum_d<-Sum_d+(S_dd[,,i]<-crossprod(S_dd[,,i]))}
		r_result<-Sum_d+crossprod(delta)
		c_result<-.C("pfda_sum_cond_dd",sum_d=double(k**2),as.integer(n),as.double(delta),as.double(S_dd), as.integer(k), 0L)
		eq_data<-all.equal(as.double(r_result),c_result$sum_d)
		if(!isTRUE(eq_data)){
			cat("failed on k=",k," i=",i,"\n")
			return(invisible(FALSE))
		}
	}
	cat("PASS\n")
	return(invisible(TRUE))
}
UT_pfda_cond_ab<-function(){
	cat('UNIT TEST - pfda_cond_ab\n')
	n<-100
	cat('testing for ka,kb<=10, N=',n,'...')
	for(ka in seq(10))for(kb in seq(10)){
		alpha<-matrix(rnorm(ka*n),nrow=n,ncol=ka)
		beta<-matrix(rnorm(kb*n),nrow=n,ncol=kb)
		S_ab<-array(rnorm(n*ka*kb,sd=3),dim=c(ka,kb,n))
		for(i in seq(n)){
			# S_ab[,,i]<-crossprod(S_ab[,,i])
			r_result<-tcrossprod(alpha[i,],beta[i,])+S_ab[,,i]
			c_result<-.C("test_pfda_cond_ab",
				ab_hat=double(ka*kb),
				as.integer(i-1),as.integer(n),
				as.double(alpha),as.integer(ka),
				as.double(beta),as.integer(kb),
				as.double(S_ab), 0L)
			eq_data<-all.equal(as.double(r_result),c_result$ab_hat)
			if(!isTRUE(eq_data)){
				cat("failed on ka=", ka,", kb=", kb,", i=",i,"\n")
				return(invisible(FALSE))
			}
		}
	}
	cat("PASS\n")
	return(invisible(TRUE))
}
UT_pfdaDual_m4<-function(){
	cat('UNIT TEST - pfdaDual_m4\n')
	cat('testing for random input...')
	N<-100
	for(ka in seq(10))for(kb in seq(10)){
		alpha<-matrix(rnorm(ka*N),nrow=N,ncol=ka)
		beta<-matrix(rnorm(kb*N),nrow=N,ncol=kb)
		S_aa<-array(rnorm(N*ka**2,sd=3),dim=c(ka,ka,N))
		S_ab<-array(rnorm(N*ka*kb,sd=3),dim=c(ka,kb,N))
		sum_aa<-matrix(0,ka,ka)
		sum_ba<-matrix(0,kb,ka)
		for(i in seq(N)){
			S_aa[,,i]<-crossprod(S_aa[,,i])
			sum_aa<-sum_aa+alpha[i,]%*%t(alpha[i,])+(S_aa[,,i])
			sum_ba<-sum_ba+beta[i,]%*%t(alpha[i,])+t(matrix(S_ab[,,i],ka,kb))
		}
		r_Lambda<-sum_ba%*%solve(sum_aa)
		c_result<-.C("pfdaDual_m4",
			as.integer(N),
			Lambda=double(kb*ka),
			as.double(alpha),as.integer(ka),
			as.double(beta), as.integer(kb),
			as.double(S_aa),as.double(S_ab),
			debuglevel=0L,
			# pfda:::c_debug(240),
			double(ka**2+ka*kb))
			eq_data<-all.equal(as.double(r_Lambda),c_result$Lambda)
		if(!isTRUE(eq_data)){
			cat("failed on ka=", ka,", kb=", kb,"\n")
			return(invisible(FALSE))
		}
	}
	cat("PASS\n")
	return(invisible(TRUE))
}
UT_pfda_matrix_outer_quadratic_form<-function(){
  cat('UNIT TEST - pfda_matrix_outer_quadratic_form\n')
  if(is.loaded('pfda_matrix_outer_quadratic_form')){
		cat("testing with random input...")
		n<-100
		a<-trunc(runif(n,max=7))+2
		b<-trunc(runif(n,max=7))+2
		for(i in seq(n)){	
			X <- matrix(runif(a[i]*b[i]),nrow=a[i],ncol=b[i])
			S <- crossprod(matrix(runif(b[i]**2),nrow=b[i],ncol=b[i]))		
			r_result<-X%*%S%*%t(X)
			c_result<-.C("pfda_matrix_outer_quadratic_form",
				Q=double(a[i]**2),
				X=as.double(X),as.integer(a[i]),as.integer(a[i]),
				as.double(S),as.integer(b[i]),
				0L,
				# pfda:::c_debug(306),
				dp=double(a[i]*b[i])
				)
			if(!isTRUE(all.equal(as.double(r_result),c_result$Q))){
				cat('Failed on attempt ',i,'\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no pfda_matrix_outer_quadratic_form function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_pfda_matrix_inner_quadratic_form<-function(){
  cat('UNIT TEST - pfda_matrix_inner_quadratic_form\n')
  if(is.loaded('pfda_matrix_inner_quadratic_form')){
		cat("testing with random input...")
		n<-100
		a<-trunc(runif(n,max=7))+2
		b<-trunc(runif(n,max=7))+2
		for(i in seq(n)){	
			X <- matrix(runif(a[i]*b[i]),nrow=a[i],ncol=b[i])
			S <- crossprod(matrix(runif(a[i]**2),nrow=a[i],ncol=a[i]))		
			r_result<-t(X) %*% S %*% X
			c_result<-{ .C("pfda_matrix_inner_quadratic_form",
				Q=double(b[i]**2),
				X=as.double(X),as.integer(b[i]),as.integer(a[i]),
				as.double(S),as.integer(a[i]),
				0L,
				# pfda:::c_debug(307),
				dp=double(a[i]*b[i])) }
			if(!isTRUE(all.equal(as.double(r_result),c_result$Q))){
				cat('Failed on attempt ',i,'\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no pfda_matrix_inner_quadratic_form function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_pfda_eigens<-function(){
  cat('UNIT TEST - pfda_eigens\n')
  if(is.loaded('pfda_eigens')){
		cat("testing with random input...")
		tryCatch({
		n<-100
		for(i in seq(n)){	
			(a<-trunc(runif(1,max=7))+4)
			(b<-trunc(runif(1,min=1,max=a-1)))
			S <- crossprod(matrix(runif(a**2),nrow=a,ncol=a))		
			r_eigen<-eigen(S)
			c_result<-.C("pfda_eigens",
				as.double(S), as.integer(a),
				vectors=double(a**2),
				values=double(a),
				k=as.integer(b),
				0L,
				# pfda:::c_debug(305),
				dp=double(8*a),ip=integer(6*a))
			if(!isTRUE(all.equal(as.double(r_eigen$values[1:b]),c_result$values[1:b]))){
				cat('Failed on attempt ',i,'\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(as.double(r_eigen$vectors[,1:b]*rep(sign(r_eigen$vectors[1,1:b]),each=a)),c_result$vectors[1:(b*a)]))){
				cat('Failed on attempt ',i,'\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n") } ,
		error=function(e){ cat('FAILED error encountered.\n');print(e) } )
    invisible(TRUE)
  } else {
    cat("no pfda_eigens function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_pfda_m5_1<-function(msgs=FALSE){
  cat('UNIT TEST - pfda_m5_1\n')
  if(is.loaded('pfda_m5_1')){
		cat("testing with random input...")
		Nsim<-100
		na<-trunc(runif(Nsim,min=21,max=101))
		pa<-trunc(runif(Nsim,min=5, max=25))
		ka<-trunc(runif(Nsim,min=1, max=pa-1))
		for(i in seq(Nsim)){	
			if(msgs)print(i)
			n<-na[i]
			p<-pa[i]
			k<-ka[i]
			newdelta<-delta<- matrix(runif(n*k),n,k)
			sigma_delta<-array(runif(k*k*n), dim=c(k,k,n)); for(j in seq(n))sigma_delta[,,j]=crossprod(sigma_delta[,,j])
			Theta <- matrix(runif(p*k),p,k)
			Sum_d <- crossprod(delta)+apply(matrix(sigma_delta,nrow=k**2),1,sum)
			Sum_d <- Sum_d/n
			r_eigen<-eigen(Theta%*%Sum_d%*%t(Theta))
			r_eigen$vectors[,r_eigen$vectors[1,]<0]<-r_eigen$vectors[,r_eigen$vectors[1,]<0]*-1
			trans<-t(r_eigen$vectors[,1:k])%*%Theta
			for(j in 1:n)newdelta[j,]<-trans%*%delta[j,]
			r_eigen$values[1:k]
			c_result<-.C("pfda_m5_1",
				as.integer(n),delta=as.double(delta),as.integer(k), as.double(sigma_delta),
				Theta=as.double(Theta), D=double(k), p=as.integer(p),trans=double(k*k),
				min_var=1e-4,
				0L,
				# pfda:::c_debug(251),
				double(k + 2*p^2 + p +  k*p + 8*p + 1000),
				integer(6*p+1000))
			gc()
			if(!isTRUE(all.equal(as.double(r_eigen$values[1:k]),c_result$D))){
				cat('Failed on attempt',i,' with D\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(as.double(r_eigen$vectors[,1:k]),c_result$Theta))){
				cat('Failed on attempt',i,' with Theta\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(as.double(newdelta),c_result$delta))){
				cat('Failed on attempt',i,' with delta\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(as.double(trans),c_result$trans))){
				cat('Failed on attempt',i,'with transformation\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no pfda_m5_1 function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_pfdaDual_m5_2<-function(){
  cat('UNIT TEST - pfdaDual_m5_2\n')
  if(is.loaded('pfdaDual_m5_2')){
	cat("testing with random input...")
	n<-100
	p<-20
	for(i in seq(100)){
		a<-trunc(runif(1,min=1,max=10))
		b<-trunc(runif(1,min=1,max=10))
		trans_f<-crossprod(matrix(rnorm(a**2),a))
		trans_g<-crossprod(matrix(rnorm(b**2),b))
		lambda<-matrix(rnorm(a*b),b,a)
		(newlambda<-trans_g%*%lambda%*%solve(trans_f))

		c_result<-.C("pfdaDual_m5_2",
			Lambda=as.double(lambda),
			trans_f=as.double(trans_f),as.integer(a),
			trans_g=as.double(trans_g),as.integer(b),
				0L,
				# pfda:::c_debug(252),
				double(a^2+a*b),
				integer(a))
		if(!all.equal(as.double(newlambda),c_result$Lambda)){
			cat('Failed on attempt',i,'\n')
			return(invisible(FALSE))
		}
	}
	cat("PASS\n")
	invisible(TRUE)
  } else {
    cat("no pfdaDual_m5_2 function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_pfda_dual_e1<-function(){
	cat('UNIT TEST - pfda_dual_e1\n')
	if(is.loaded('pfda_dual_e1')){
		cat("testing with random input...")
		N<-100
		for(i in 1:N){
			a<-trunc(runif(1, min=1, max = 10))
			b<-trunc(runif(1, min=1, max = 10))
			n<-trunc(runif(1, min=1, max = 20))
			Da<-(sort(rexp(a),T))
			Db<-(sort(rexp(b),T))
			L<-matrix(rnorm(b*a),b,a)
			epsilon<-rexp(1)
			xi<-rexp(1)
			phi<-matrix(rnorm(a*n),n,a)
			psi<-matrix(rnorm(b*n),n,b)
			
			Sni<-solve(diag(x=Db,b,b)-L%*%diag(x=Da,a,a)%*%t(L))
			
			zeta_aa<- solve(diag(x=Da,a,a))+t(L)%*%Sni%*%L+crossprod(phi)/epsilon
			zeta_ab<- -t(L)%*%Sni
			zeta_bb<- Sni + crossprod(psi)/xi
			
			c_result<-.C("pfda_dual_e1",
				as.integer(n),as.integer(n),as.double(L),
				as.double(Da),as.double(phi),as.double(epsilon),as.integer(a),
				as.double(Db),as.double(psi),as.double(xi)     ,as.integer(b),
				zeta_aa=double(a**2),zeta_ab=double(a*b),zeta_bb=double(b**2),
				0L,
				# pfda:::c_debug(270),
				double(b**2 + a**2 + max(a*b,2*(b**2))),
				integer(b))
			if(!isTRUE(all.equal(as.double(zeta_aa)[upper.tri(zeta_aa,T)],c_result$zeta_aa[upper.tri(zeta_aa,T)]))){
				cat('Failed on attempt',i,' for zeta_aa\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(as.double(zeta_ab),c_result$zeta_ab))){
				cat('Failed on attempt',i,' for zeta_ab\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(as.double(zeta_bb)[upper.tri(zeta_bb,T)],c_result$zeta_bb[upper.tri(zeta_bb,T)]))){
				cat('Failed on attempt',i,' for zeta_bb\n')
				return(invisible(FALSE))
			}
		}
	cat("PASS\n")
	} else {
	cat("no pfda_dual_e1 function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfda_gen_e2_1<-function(){
	cat('UNIT TEST - pfda_gen_e2_1\n')
	if(is.loaded('pfda_gen_e2_1')){
		cat("testing with random input...")
		N<-100
		for(i in 1:N){
			seed<-.Random.seed
			a<-trunc(runif(1, min=1, max = 10))
			c<-trunc(runif(1, min=1, max = 10))
			useouter<-FALSE
			A<-crossprod(matrix(rnorm(a**2),a,a))
			B<-rnorm(a*c)
			dim(B)<-c(c,a) 
			C<-crossprod(matrix(rnorm(c**2),c,c))

			D<-solve(A-t(B)%*%solve(C)%*%B)

			co<-.C("pfda_gen_e2_1",
				as.integer(useouter),
				as.integer(a),
				as.integer(c),
				A=as.double(A),
				B=as.double(B),
				C_inv=as.double(C),
				D=as.double(D),
				0L,
				# pfda:::c_debug(281),
				double(max( 2*a^2, a*c, 2*c^2)),
				integer(max(a,c)))
			if(!isTRUE(eqd<-all.equal(as.double(D),co$D))){
				write(seed,file='errseed.txt')
				cat('Failed on attempt',i,' for inner\n',eqd,'\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(eqd<-all.equal(as.double(solve(C)),co$C_inv))){
				write(seed,file='errseed.txt')
				cat('Failed on attempt',i,' for  C_inv (inner)\n',eqd,'\n')
				return(invisible(FALSE))
			}
			useouter<-TRUE
			dim(B) <- c(a,c) 
			D<-solve(A-B%*%solve(C)%*%t(B))
			co<-.C("pfda_gen_e2_1",
				as.integer(useouter),
				as.integer(a),
				as.integer(c),
				A=as.double(A),
				B=as.double(B),
				C_inv=as.double(C),
				D=as.double(D),
				0L,
				# pfda:::c_debug(281),
				double(max( 2*a^2, a*c, 2*c^2)),
				integer(max(a,c)))
			all.equal(as.double(D),co$D)
			if(!isTRUE(eqd<-all.equal(as.double(D),co$D,tolerance = .Machine$double.eps ^ 0.5 * 100))){
				write(seed,file='errseed.txt')
				cat('Failed on attempt',i,' for outer\n',eqd,'\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(eqd<-all.equal(as.double(solve(C)),co$C_inv,tolerance = .Machine$double.eps ^ 0.5 * 100))){
				write(seed,file='errseed.txt')
				cat('Failed on attempt',i,' for  C_inv (outer)\n',eqd,'\n')
				return(invisible(FALSE))
			}
		}
	cat("PASS\n")
	} else {
	cat("no pfda_gen_e2_1 function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfdaDual_e2<-function(){
	cat('UNIT TEST - pfdaDual_e2\n')
	if(is.loaded('pfdaDual_e2')){
		cat("testing with random input...")
		N<-100
		for(i in 1:N){
			a<-trunc(runif(1, min=1, max = 10))
			b<-trunc(runif(1, min=1, max = 10))
			zeta_aa<-crossprod(matrix(rnorm(a**2),a,a))
			zeta_bb<-crossprod(matrix(rnorm(b**2),b,b))
			zeta_ab<-matrix(rnorm(a*b),a,b)

			sigma_aa<- solve(zeta_aa-zeta_ab%*%solve(zeta_bb)%*%t(zeta_ab))
			sigma_ab<- -sigma_aa%*%zeta_ab%*%solve(zeta_bb)
			sigma_bb<- solve(zeta_bb-t(zeta_ab)%*%solve(zeta_aa)%*%zeta_ab)

			co<-.C("pfdaDual_e2",
				as.integer(a),
				as.integer(b),
				zeta_aa=as.double(zeta_aa),
				zeta_ab=as.double(zeta_ab),
				zeta_bb=as.double(zeta_bb),
				sigma_aa=as.double(sigma_aa),
				sigma_ab=as.double(sigma_ab),
				sigma_bb=as.double(sigma_bb),
				0L,
				# c(280L,pfda:::make_c_debug(280, 400)),
				double(b^2+2*max(a,b)^2),
				integer(max(a,b)))
			if(!isTRUE(eqd<-all.equal(as.double(sigma_aa),co$sigma_aa,tolerance = .Machine$double.eps ^ 0.5 * 100))){
				cat('Failed on attempt',i,' for  sigma_aa\n',eqd,'\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(eqd<-all.equal(as.double(sigma_ab),co$sigma_ab,tolerance = .Machine$double.eps ^ 0.5 * 100))){
				cat('Failed on attempt',i,' for sigma_ab\n',eqd,'\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(eqd<-all.equal(as.double(sigma_bb),co$sigma_bb,tolerance = .Machine$double.eps ^ 0.5 * 100))){
				cat('Failed on attempt',i,' for  sigma_bb\n',eqd,'\n')
				return(invisible(FALSE))
			}
		}
	cat("PASS\n")
	} else {
	cat("no pfdaDual_e2 function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfda_dual_e3_1<-function(){
	cat('UNIT TEST - pfda_dual_e3_1\n')
	if(is.loaded('pfda_dual_e3_1')){
		cat("testing with random input...")
		N<-100
		for(i in 1:N){
			a<-trunc(runif(1, min=1, max = 10))
			b<-trunc(runif(1, min=1, max = 10))
			n<-trunc(runif(1, min=4, max=20))
			M<-trunc(runif(1, min=n, max = 1000))

			u<-rnorm(n)
			v<-rnorm(n)

			epsilon=rexp(1)
			xi=rexp(1)

			AA<-crossprod(matrix(rnorm(a**2),a,a))
			AB<-matrix(rnorm(a*b),a,b)

			phi=matrix(rnorm(M*a),M,a)
			psi=matrix(rnorm(M*b),M,b)

			Xa=tcrossprod(AA,phi[1:n,,drop=FALSE])/epsilon
			Xb=tcrossprod(AB,psi[1:n,,drop=FALSE])/xi

			mu=as.vector(Xa%*%u+Xb%*%v)

			co<-.C("pfda_dual_e3_1",
				mu=double(a),
				1L,
				as.integer(M),
				as.integer(n),
				as.double(epsilon),
				as.double(AA),
				as.double(phi),
				as.double(u),
				as.integer(a),
				as.double(xi),
				as.double(AB),
				as.integer(FALSE),
				as.double(psi),
				as.double(v),
				as.integer(b),
				0L,
				# c(291L,pfda:::make_c_debug(291)),
				double(2*n*a**2)
				)	
			if(!isTRUE(eqd<-all.equal(as.double(mu),co$mu))){
				cat('Failed on attempt',i,' for  trans=FALSE\n',eqd,'\n')
				return(invisible(FALSE))
			}
			dim(AB)<-dim(AB)[2:1]
			Xb=tcrossprod(t(AB),psi[1:n,,drop=FALSE])/xi
			mu=Xa%*%u+Xb%*%v# 
			co<-.C("pfda_dual_e3_1",
				mu=double(a),
				1L,
				as.integer(M),
				as.integer(n),
				as.double(epsilon),
				as.double(AA),
				as.double(phi),
				as.double(u),
				as.integer(a),
				as.double(xi),
				as.double(AB),
				as.integer(TRUE),
				as.double(psi),
				as.double(v),
				as.integer(b),
				0L,
				# c(291L,pfda:::make_c_debug(291)),
				double(2*n*a**2)
				)	
			if(!isTRUE(eqd<-all.equal(as.double(mu),co$mu))){
				cat('Failed on attempt',i,' for  trans=TRUE\n',eqd,'\n')
				return(invisible(FALSE))
			}
		}
	cat("PASS\n")
	} else {
	cat("no pfda_dual_e3_1 function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfda_dual_e3<-function(){
	cat('UNIT TEST - pfda_dual_e3\n')
	if(is.loaded('pfda_dual_e3')){
		cat("testing with random input...")
		N<-100
		for(i in 1:N){
			a<-trunc(runif(1, min=1, max = 10))
			b<-trunc(runif(1, min=1, max = 10))
			n<-trunc(runif(1, min=4, max=20))
			M<-trunc(runif(1, min=n, max = 1000))

			u<-rnorm(n)
			v<-rnorm(n)

			epsilon=rexp(1)
			xi=rexp(1)

			sigma_aa<-crossprod(matrix(rnorm(a**2),a,a))
			sigma_bb<-crossprod(matrix(rnorm(b**2),b,b))
			sigma_ab<-matrix(rnorm(a*b),a,b)

			phi=matrix(rnorm(M*a),M,a)
			psi=matrix(rnorm(M*b),M,b)

			Xa=tcrossprod(sigma_aa,phi[1:n,,drop=FALSE])/epsilon
			Xb=tcrossprod(sigma_ab,psi[1:n,,drop=FALSE])/xi
			alpha=as.vector(Xa%*%u+Xb%*%v)

			Ya=tcrossprod(t(sigma_ab),phi[1:n,,drop=FALSE])/epsilon
			Yb=tcrossprod(sigma_bb,psi[1:n,,drop=FALSE])/xi
			beta=as.vector(Ya%*%u+Yb%*%v)

			co<-.C("pfda_dual_e3",
				alpha=double(a),
				beta=double(b),
				as.double(sigma_aa),
				as.double(sigma_ab),
				as.double(sigma_bb),
				as.double(phi),
				as.double(psi),
				as.double(epsilon),
				as.double(xi),
				as.double(u),
				as.double(v),
				as.integer(a),
				as.integer(b),
				as.integer(n),
				as.integer(1),
				as.integer(M),
				0L,
				#c(291L,pfda:::make_c_debug(290,291)),
				double(2*n*max(a,b))
				)
			if(!isTRUE(eqd<-all.equal(as.double(alpha),co$alpha))){
				cat('Failed on attempt',i,' on alpha\n',eqd,'\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(eqd<-all.equal(as.double(beta),co$beta))){
				cat('Failed on attempt',i,' on beta\n',eqd,'\n')
				return(invisible(FALSE))
			}
		}
	cat("PASS\n")
	} else {
	cat("no pfda_dual_e3 function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfdaDual_e<-function(){
	cat('UNIT TEST - pfdaDual_e\n')
	if(is.loaded('pfdaDual_e')){
		cat("testing with random input...")
		Nsim<-100
		for(i in 1:Nsim){
			N<-trunc(runif(1,20,100))
			nobs<-trunc(runif(N,4,10))
			M<-sum(nobs)

			y<-rnorm(M)
			z<-rnorm(M)
			ka<-trunc(runif(1,1,6))
			kb<-trunc(runif(1,1,6))
			p <-trunc(runif(1,10,20))
			B <-matrix(rnorm(p*M),M,p)
			tm<-rnorm(p)
			tn<-rnorm(p)
			tf<-matrix(rnorm(p*ka),p,ka);tf<-tf*rep(ifelse(tf[1,]<0,-1,1),each=p)
			tg<-matrix(rnorm(p*kb),p,kb);tg<-tg*rep(ifelse(tg[1,]<0,-1,1),each=p)
			Da<-(sort(rexp(ka),T))
			Db<-(sort(rexp(kb),T))
			Lambda<-matrix(rnorm(kb*ka),kb,ka)
			epsilon<-rexp(1)
			xi<-rexp(1)
			Alpha<-matrix(rnorm(ka*N),N,ka)
			Beta <-matrix(rnorm(kb*N),N,kb)
			Sigma_aa<-array(0,dim=c(ka,ka,N))
			Sigma_ab<-array(0,dim=c(ka,kb,N))
			Sigma_bb<-array(0,dim=c(kb,kb,N))

			btf<-B%*%tf
			btg<-B%*%tg

			Ry<- as.double(y-B%*%tm)
			Rz<- as.double(z-B%*%tn)

			for(i in 1:N){ #R computations
				a<-ka
				b<-kb
				n<-nobs[i]
				L<-Lambda
				phi<-btf[c(0,cumsum(nobs))[i]+1:n,,drop=FALSE]
				psi<-btg[c(0,cumsum(nobs))[i]+1:n,,drop=FALSE]
				u<-Ry[c(0,cumsum(nobs))[i]+1:n,drop=FALSE]
				v<-Rz[c(0,cumsum(nobs))[i]+1:n,drop=FALSE]

				Sni<-solve(diag(x=Db,b,b)-L%*%diag(x=Da,a,a)%*%t(L))

				zeta_aa<- solve(diag(x=Da,a,a))+t(L)%*%Sni%*%L+crossprod(phi)/epsilon
				zeta_ab<- -t(L)%*%Sni
				zeta_bb<- Sni + crossprod(psi)/xi

				sigma_aa<- solve(zeta_aa-zeta_ab%*%solve(zeta_bb)%*%t(zeta_ab))
				sigma_ab<- -sigma_aa%*%zeta_ab%*%solve(zeta_bb)
				sigma_bb<- solve(zeta_bb-t(zeta_ab)%*%solve(zeta_aa)%*%zeta_ab)

				Xa=tcrossprod(sigma_aa,phi[1:n,,drop=FALSE])/epsilon
				Xb=tcrossprod(sigma_ab,psi[1:n,,drop=FALSE])/xi
				alpha=as.vector(Xa%*%u+Xb%*%v)

				Ya=tcrossprod(t(sigma_ab),phi[1:n,,drop=FALSE])/epsilon
				Yb=tcrossprod(sigma_bb,psi[1:n,,drop=FALSE])/xi
				beta=as.vector(Ya%*%u+Yb%*%v)

				Alpha[i,]<-alpha
				Beta[i,] <-beta
				Sigma_aa[,,i]<-sigma_aa
				Sigma_ab[,,i]<-sigma_ab
				Sigma_bb[,,i]<-sigma_bb
			}
			co<-{ .C("pfdaDual_e", #C Computations
				as.double(y),
				as.double(z),
				as.integer(nobs),
				as.integer(M),
				as.integer(N),
				as.integer(ka),
				as.integer(kb),
				as.double(B),
				as.integer(p),
				as.double(tm),
				as.double(tn),
				as.double(tf),
				as.double(tg),
				as.double(Da),
				as.double(Db),
				as.double(Lambda),
				as.double(epsilon),
				as.double(xi),
				Alpha=double(N*ka),
				Beta =double(N*kb),
				Sigma_aa=double(ka*ka*N),
				Sigma_ab=double(ka*kb*N),
				Sigma_bb=double(kb*kb*N),
				0L,
				# pfda:::c_debug(260,261,270,280,281,290,291),
				double(2*M+2*M*max(ka,kb)+3*max(ka,kb)^2 +  max( 2*max(nobs)*max(ka,kb),4*max(ka,kb)^2) ),
				integer(max(ka,kb))
				)
			}
			if(!isTRUE(eqd<-all.equal(as.double(Alpha),co$Alpha))){
				cat('Failed on attempt',i,' on alpha\n',eqd,'\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(eqd<-all.equal(as.double(Beta),co$Beta))){
				cat('Failed on attempt',i,' on beta\n',eqd,'\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(eqd<-all.equal(as.double(Sigma_aa),co$Sigma_aa))){
				cat('Failed on attempt',i,' on Sigma_aa\n',eqd,'\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(eqd<-all.equal(as.double(Sigma_ab),co$Sigma_ab))){
				cat('Failed on attempt',i,' on Sigma_ab\n',eqd,'\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(eqd<-all.equal(as.double(Sigma_bb),co$Sigma_bb))){
				cat('Failed on attempt',i,' on Sigma_bb\n',eqd,'\n')
				return(invisible(FALSE))
			}
		}
	cat("PASS\n")
	} else {
	cat("no pfdaDual_e function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfda_sym_inverse<-function(){
	cat('UNIT TEST - pfda_sym_inverse\n testing with random input ... ')
	for(i in 1:100){
		n<-trunc(runif(1,1,21))
		A<-crossprod(matrix(runif(n**2),n,n))

		cr<-.C('pfda_sym_inverse',Ai=as.double(A),as.integer(n),sr=0L, 0L, double(2*n**2), integer(n))

		if(cr$sr||!isTRUE(all.equal(as.double(solve(A)),cr$Ai))){
				cat('Failed on attempt',i,'\n')
				return(invisible(FALSE))
		}
	}
	cat("PASS\n")
}
UT_pfda_computebtb<-function(){
  cat('UNIT TEST - pfda_computebtb\n')
  if(is.loaded('pfda_computebtb')){
		cat("testing with random input...")
		Nsim<-100
		for(i in seq(Nsim)){	
			p     = trunc(runif(1,min=6,max=20))
			N     = trunc(runif(1,min=20,max=100))
			nobs  = trunc(runif(N,min=3,max=10))
			M     = sum(nobs)
			Bl    = sapply(nobs,function(n) matrix(rnorm(p*n),n,p))
			B     = NULL ; for(i in 1:N)B<-rbind(B, Bl[[i]])

			unlist(lapply(Bl,function(x){
			crossprod(x)->y
			y[lower.tri(y)]<-0
			y }))->btb
			dim(btb)<-c(p,p,N)
				
			cr<-.C('pfda_computebtb',
				btb=double(N*p*p),
				as.integer(N),
				as.double(B),
				as.integer(M),
				as.integer(p),
				as.integer(nobs),
				0L)$btb
			dim(cr)<-c(p,p,N)

			if(!isTRUE(all.equal(btb,cr))){
				cat('Failed on attempt',i,'\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no pfda_computebtb function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_pfda_uppersymadd<-function(){
  cat('UNIT TEST - pfda_uppersymadd\n')
  if(is.loaded('pfda_uppersymadd')){
		cat("testing with random input...")
		Nsim<-100
		for(i in seq(Nsim)){	
			n<-5#trunc(runif(1,min=6,max=20))
			x<-crossprod(matrix(rnorm(n**2),n,n))
			x[lower.tri(x)]<-0

			y<-crossprod(matrix(rnorm(n**2),n,n))
			y[lower.tri(y)]<-0

			alpha<-rnorm(1)

			r<-x+alpha*y

			cr<-.C('pfda_uppersymadd',
				as.integer(n),
				x=as.double(x),as.integer(n),
				as.double(alpha),
				y=as.double(y),as.integer(n),
				0L)
			dim(cr$x)<-c(n,n)

			if(!isTRUE(all.equal(r,cr$x))){
				cat('Failed on attempt',i,'\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no pfda_uppersymadd function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
{ #pfda_gen_e_eta
X1_pfda_gen_e_eta<-expression({
	eval(sim.dual.parameters)

	C<-.C('pfda_gen_e_eta',	double(kb**2),Db,lambda,Da,ka,kb,0L,double(ka**2+ka*kb))

	Seta <- diag(Db,length(Db))-lambda %*% tcrossprod(diag(Da,ka),lambda)

	if(!isTRUE(all.equal(matrix(C[[1]],kb,kb),Seta)))stop("Failed")
})
X_pfda_gen_e_eta<-expression({
	replicate(100,eval(X1_pfda_gen_e_eta))
})
UT_pfda_gen_e_eta<-function(){
	if(class(try(eval(X_pfda_gen_e_eta)))=="try-error")return('Failed')
}
}

X_gen_orthog<-expression({
eval(sim.single.c)
avgaa<-crossprod(matrix(rnorm(ka**2),ka,ka))

R<-.gen.orthog(tf,alpha,avgaa*N)


dp = double(9*p + 3*p^2 + N*ka)
ip = integer(6*p)
dl=0L
C<-.C('gen_orthog', tf=tf, alpha=alpha, Da=Da, trans=double(ka^2), avgaa=avgaa, N=N, ka=ka, p=p, dl=dl, dp=dp, ip=ip)

if(!isTRUE(all.equal(matrix(C$trans,ka,ka),R$transformation)))stop('failed on transformation')
if(!isTRUE(all.equal(as.vector(C$alpha),as.vector(R$d))))stop('failed on alpha')
if(!isTRUE(all.equal(as.vector(C$tf),as.vector(R$tf))))stop('failed on tf')
if(!isTRUE(all.equal(as.vector(C$Da),as.vector(R$D))))stop('failed on Da')
})
UT_gen_orthog<-UT_generate(X_gen_orthog)

}
{ # single steps
UT_pfda_m1<-function(){
  cat('UNIT TEST - pfda_m1\n')
  if(is.loaded('pfda_m1')){
		cat("testing with random input...")
		Nsim<-100
		for(i in seq(Nsim)){	
			p     = trunc(runif(1,min=6,max=20))
			N     = trunc(runif(1,min=20,max=100))
			nobs  = trunc(runif(N,min=3,max=10))
			M     = sum(nobs)
			y     = rnorm(M)
			ka    = trunc(runif(1,1,p))
			Bl    = sapply(nobs,function(n) matrix(rnorm(p*n),n,p))
			B     = NULL ; for(i in 1:N)B<-rbind(B, Bl[[i]])
			lm    = 1
			ISD   = crossprod(matrix(rnorm(p**2),p,p))
			tf    = matrix(rnorm(p*ka),p,ka)
			seps  = rexp(1)
			alpha = matrix(rnorm(N*ka), N, ka)

			delta = 1e-4
			tm = rnorm(p) 
			Sigma_aa = array(rnorm(ka*ka*N),dim=c(ka,ka,N)) ; for(i in 1:N) Sigma_aa[,,i]<- crossprod(Sigma_aa[,,i])

			tr<-function(x)if(is.matrix(x)) sum(diag(x)) else sum(x)
			sum_tr=0
			Ry=numeric(length(y))
			phi=B%*%tf
			for(i in 1:N){
				idx=seq(c(0,cumsum(nobs))[i]+1,cumsum(nobs)[i])
				Ry[idx] = y[idx] - Bl[[i]] %*% tf %*% alpha[i,]
				Spp = Sigma_aa[,,i] %*%t(phi[idx,]) %*% phi[idx,] 
				sum_tr = sum_tr + tr(Spp)
			}
			Ry<- Ry-B%*%tm
			seps=(crossprod(Ry)+sum_tr)/M

			cr<-.C('pfda_m1',
				seps     = double(1), 
				y        = as.double(y),
				nobs     = as.integer(nobs),
				M        = as.integer(M),
				N        = as.integer(N),
				ka       = as.integer(ka),
				B        = as.double(B),
				p        = as.integer(p),
				delta    = as.double(delta), 
				tm       = as.double(tm),
				tf       = as.double(tf),
				Alpha    = as.double(alpha),
				Sigma_aa = as.double(Sigma_aa), 
				0L,
				# pfda:::c_debug(110,111,112,306,307),
				dp=double(M + M*ka + 2*ka^2))$seps

			if(!isTRUE(all.equal(as.double(seps),cr)) || cr==1e-4){
				cat('Failed on attempt',i,'\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no pfda_m1 function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_pfda_m2<-function(){
  cat('UNIT TEST - pfda_m2\n')
  if(is.loaded('pfda_m2')){
		cat("testing with random input...")
		Nsim<-100
		for(i in seq(Nsim)){	
			p     = trunc(runif(1,min=6,max=20))
			N     = trunc(runif(1,min=20,max=100))
			nobs  = trunc(runif(N,min=3,max=10))
			M     = sum(nobs)
			y     = rnorm(M)
			ka    = trunc(runif(1,1,p))
			Bl    = sapply(nobs,function(n) matrix(rnorm(p*n),n,p))
			B     = NULL ; for(i in 1:N)B<-rbind(B, Bl[[i]])
			lm    = 1
			ISD   = crossprod(matrix(rnorm(p**2),p,p))
			tf    = matrix(rnorm(p*ka),p,ka)
			seps  = rexp(1)
			alpha = matrix(rnorm(N*ka), N, ka)

			Ry=numeric(length(y))
			for(i in 1:N){
				idx=seq(c(0,cumsum(nobs))[i]+1,cumsum(nobs)[i])
				Ry[idx] = y[idx] - Bl[[i]] %*% tf %*% alpha[i,]
			}

			right = crossprod(B,Ry)
			BtB=crossprod(B)+seps*lm*ISD
			tm=solve(BtB,right)

			cr<-.C('pfda_m2',
				tm    = double(p),
				y     = as.double(y),
				nobs  = as.integer(nobs),
				M     = as.integer(M),
				N     = as.integer(N),
				ka    = as.integer(ka),
				B     = as.double(B),
				p     = as.integer(p),
				lm    = as.double(lm),
				ISD   = as.double(ISD),
				tf    = as.double(tf),
				seps  = as.double(seps),
				Alpha = as.double(alpha),
				0L,
				# pfda:::c_debug(120,121),
				dp=double(p^2+ M+ M*ka))$tm
			
			if(!isTRUE(all.equal(as.double(tm),cr))){
				cat('Failed on attempt',i,'\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no pfda_m2 function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_pfda_m3<-function(msgs = FALSE){
  cat('UNIT TEST - pfda_m3\n')
  if(is.loaded('pfda_m3')){
		cat("testing with random input...")
		Nsim<-100
		for(i in seq(Nsim)){
			if(msgs)cat('attempt',i,"\n")
			p     = trunc(runif(1,min=6,max=20))
			N     = trunc(runif(1,min=20,max=100))
			nobs  = trunc(runif(N,min=3,max=10))
			M     = sum(nobs)
			y     = rnorm(M)
			ka    = trunc(runif(1,1,p))
			Bl    = lapply(nobs,function(n) matrix(rnorm(p*n),n,p))
			B     = NULL ; for(i in 1:N)B<-rbind(B, Bl[[i]])
			lf    = rchisq(1,1)
			ISD   = crossprod(matrix(rnorm(p**2),p,p))
			tf    = matrix(rnorm(p*ka),p,ka)
			seps  = rexp(1)
			alpha = matrix(rnorm(N*ka), N, ka)

			unlist(lapply(Bl,function(x){
				crossprod(x)->y
				y[lower.tri(y)]<-0
				y }))->btb
			dim(btb)<-c(p,p,N)

			tm = rnorm(p) 
			Sigma_aa = array(rnorm(ka*ka*N),dim=c(ka,ka,N)) ; for(i in 1:N) Sigma_aa[,,i]<- crossprod(Sigma_aa[,,i])

			tfold<-tf
			Ry<- y-B%*%tm
			for(ec in 1:ka){
				left = matrix(0,p,p)
				right = numeric(p)
				for(sn in 1:N){				
					idx=seq(c(0,cumsum(nobs))[sn]+1,cumsum(nobs)[sn])
					rho = tcrossprod(alpha[sn,])+Sigma_aa[,,sn]
					left = left + rho[ec,ec]*crossprod(Bl[[sn]])
					
					sum_cols = numeric(p)
					for(cn in 1:ka)if(cn!=ec){
						sum_cols = sum_cols + rho[ec,cn]*tf[,cn]
					}
					right = right + t(Bl[[sn]])%*%Ry[idx]*alpha[sn,ec]-crossprod(Bl[[sn]])%*%sum_cols
				}
				left = left + seps * lf * ISD
				tf[,ec]<-solve(left,right)
			}


			if(msgs)cat("entering compiled code\n")
			cr<-.C('pfda_m3',
				tf       = as.double(tfold),
				y        = as.double(y),
				nobs     = as.integer(nobs),
				M        = as.integer(M),
				N        = as.integer(N),
				ka       = as.integer(ka),
				B        = as.double(B),
				p        = as.integer(p),
				lf       = as.double(lf),
				ISD      = as.double(ISD),
				tm       = as.double(tm),
				seps     = as.double(seps), 
				Alpha    = as.double(alpha),
				Sigma_aa = as.double(Sigma_aa), 
				btb      = as.double(btb), 
				0L,
				# pfda:::c_debug(301,134:137),
				dp=double(M + ka^2 * N + p^2 + p ),
				ip=double(p))$tf
			if(msgs)cat("returned from compiled code\n")
			dim(cr)<-c(p,ka)
			
			if(!isTRUE(all.equal(tf,cr))){
				cat('Failed on attempt',i,'\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no pfda_m3 function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_pfdaSingle_e <-function(){
  cat('UNIT TEST - pfdaSingle_e\n')
  if(is.loaded('pfdaSingle_e')){
		cat("testing with random input...")
		Nsim<-100
		for(i in seq(Nsim)){	
			p     = trunc(runif(1,min=6,max=20))
			N     = trunc(runif(1,min=20,max=100))
			nobs  = trunc(runif(N,min=3,max=10))
			M     = sum(nobs)
			y     = rnorm(M)
			ka    = trunc(runif(1,1,p))
			Bl    = lapply(nobs,function(n) matrix(rnorm(p*n),n,p))
			B     = NULL ; for(i in 1:N)B<-rbind(B, Bl[[i]])
			lf    = rchisq(1,1)
			ISD   = crossprod(matrix(rnorm(p**2),p,p))
			tm    = rnorm(p) 
			tf    = eigen(tcrossprod(matrix(rnorm(p*ka),p,ka)))$vectors[,seq(ka),drop=FALSE]
			tf    = tf*rep(sign(tf[1,]),each=p)
			Da    = sort(rexp(ka),T)
			seps  = rexp(1)
			alpha = matrix(0, N, ka)
			Sigma_aa = array(0,dim=c(ka,ka,N))

			Ry    <- y-B%*%tm
			phi   <- B%*%tf
			sn=2
			for(sn in 1:N){				
				idx=seq(c(0,cumsum(nobs))[sn]+1,cumsum(nobs)[sn])
				
				Sigma_aa[,,sn] <- solve(solve(diag(x=Da,ka,ka)) + crossprod(phi[idx,])/seps)
				alpha   [sn,]  <- Sigma_aa[,,sn]%*%crossprod(phi[idx,],Ry[idx,])/seps
			}

			dpl = M + M*ka + N*ka + 2* ka^2 
			cr<-.C('pfdaSingle_e',
				alpha   =double(N*ka),
				Sigma_aa=double(N*ka*ka),
				as.double  (y),
				as.integer (nobs),
				as.integer (M),
				as.integer (N),
				as.integer (ka),
				as.double  (B),
				as.integer (p),
				as.double  (tm),
				as.double  (tf),
				as.double  (Da),
				as.double  (seps),
				0L,
				# pfda:::c_debug(160,161),
				double(dpl) , integer(ka)
				)

			if(!isTRUE(all.equal(as.double(alpha),cr$alpha))){
				cat('Failed on attempt',i,' with Alpha\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(as.double(Sigma_aa),cr$Sigma_aa))){
				cat('Failed on attempt',i,' with Sigma_aa\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
    invisible(TRUE)
  } else {
    cat("no pfda_m2 function found.\n")
    cat("Failed\n")
		invisible(FALSE)
  }
}
UT_pfda_s_i<-function(){
	fname='pfda_s_i'
	cat('UNIT TEST -',fname,'\n')
	if(is.loaded(fname)){
		cat("testing with random input...")
		Nsim<-100
		# p=print
		p=invisible
		
		for(simnum in 1:Nsim){
			N<-trunc(runif(1,10,50))
			nobs<-trunc(runif(N,3,10))
			M<-sum(nobs)
			y<-rbinom(M,1,.5)
			t<-runif(M,0,100)
			knots<-expand.knots(0:10)*10
			id<-factor(rep(seq_len(N),nobs))
			control=pfdaControl()
			k=1
			penalties=c(10,100)
			kr=10
			delta=1e-3
			response	<- y
			domain	<- t
			subject	<- id 
			obase<-OBasis(knots)
			B	<-Bmatrix<-evaluate(obase,domain)
			K	<-Kmatrix <- OuterProdSecondDerivative(obase)
			p	<- as.integer(dim(obase)[2])
			w<-response[[1]]
			R <- .single.c.i(y=y, B=B, subject=id, k=k, min.v=delta)
			btb<-tapply(1:M,subject,function(indx)crossprod(B[indx,]))
			C <- .C('pfda_s_i',
				tm=double(p),
				tf=double(p*k),
				alpha=double(N*k),
				Da=double(k),
				aa=double(k*k*N),
				sigma=double(0),
				as.double(y), 	as.integer(nobs),	as.integer(M), 	as.integer(N), 	as.integer(k), 
				as.double(B), 	as.double(unlist(btb)),	as.integer(p), 	as.double(delta),
				pfda:::Cdebug(), double(2*p^2+M+p*N+8*p),	integer(6*p))

			if(!isTRUE(all.equal(C$tm,as.double(R$tm)))){
				cat('Failed on attempt',simnum,'for values for tm\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(C$tf,as.double(R$tf)))){
				cat('Failed on attempt',simnum,'for values for tf\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(C$alpha,as.double(R$alpha)))){
				cat('Failed on attempt',simnum,'for values for alpha\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(C$Da,R$Da))){
				cat('Failed on attempt',simnum,'for values for Da\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(unlist(R$aa),C$aa))){
				cat('Failed on attempt',simnum,'for values for aa\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
	} else {
	cat("no",fname,"function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
}
{ # Dual Functions
{ # dual_gen_sigmas
X_dual_gen_sigmas<-expression({
{ # setup
# library(pfda);i=1
# set.seed(123)
# source("unittests.R")
# source("R_steps.R")
eval(sim.dual.cc)
s.eps  = rexp(1,1)
s.xi   = rexp(1,1)
}

ix = subject==1
R<-.gen.dual.sigmas(B[ix,],B[ix,],tf,tg,lambda,Da,Db,s.eps,s.xi)

dpl =  7 * max(ka,kb)^2
dp  = double(dpl)
C<-{ .C('dual_gen_sigmas',
	Saa	  =   double(ka*ka),
	Sab   =   double(ka*kb),
	Sbb   =   double(kb*kb),
	phi   =as.double(B%*%tf),
	psi   =as.double(B%*%tg),
	lambda=as.double(lambda),
	Da    =as.double(Da),
	Db    =as.double(Db),
	sep   =as.double(s.eps),
	sxi   =as.double(s.xi),
	M     =as.integer(M),
	ni    =as.integer(nobs),
	ka    =as.integer(ka),
	kb    =as.integer(kb),
	pfda:::Cdebug(0),	as.double(dp), integer(max(ka,kb)),
	DUP=FALSE)}
gc()

if(!isTRUE(all.equal(C$Saa,as.vector(R$Saa))))stop("failed")
if(!isTRUE(all.equal(C$Sab,as.vector(R$Sab))))stop("failed")
if(!isTRUE(all.equal(C$Sbb,as.vector(R$Sbb))))stop("failed")

})
UT_dual_gen_sigmas<-function(n=100){ if(class(try(replicate(n,eval(X_dual_gen_sigmas))))=='try-error'){ cat("Failed\n"); return(invisible(FALSE))} ; invisible(TRUE) }
}

}
{ # binary Conditionals
UT_pfda_bin_cond_aa<-function(){
	fname='pfda_bin_cond_aa'
	cat('UNIT TEST -',fname,'\n')
	if(is.loaded(fname)){
		cat("testing with random input...")
		Nsim<-100
		# p=print
		p=invisible
		
		for(simnum in 1:Nsim){

			N=trunc(runif(1,min=6,max=50))
			nobs=trunc(runif(N,min=2,max=8))
			M=sum(nobs)
			ka=trunc(runif(1,min=1,max=5))

			w = rnorm(M)
			ww = sapply(nobs,function(n)crossprod(matrix(rnorm(n**2),n,n)))
			# ww = array(rnorm(4**2*N),dim=c(4,4,N))
			# for(obsnum in 1:N)ww[,,obsnum]<-crossprod(ww[,,obsnum])
			n = nobs[1]
			S1 = matrix(rnorm(ka*N),ka,N)
			phi = matrix(rnorm(ka*M),M,ka)
			Sigma_aa = array(rnorm(ka**2*N),dim=c(ka,ka,N))
			for(obsnum in 1:N)Sigma_aa[,,obsnum]<-crossprod(Sigma_aa[,,obsnum])

			for(obsnum in 1:N){
				(Sigma_aa_i = matrix(Sigma_aa[,,obsnum],ka,ka))
				(phii = phi[(sum(nobs[0:(obsnum-1)])+1):(sum(nobs[0:obsnum])),])
				(wi = w[(sum(nobs[0:(obsnum-1)])+1):(sum(nobs[0:obsnum]))])
				(wwi = ww[[obsnum]])
				(S1_i = S1[,obsnum])

				aa_hat = 	{
						p(p(p(p(Sigma_aa_i %*% p(t(phii) %*% wwi %*% phii) %*% t(Sigma_aa_i) )
						+ Sigma_aa_i%*%t(phii)%*%wi%*%t(S1_i) + S1_i%*%t(Sigma_aa_i%*%t(phii)%*%wi) ) 
						+ tcrossprod(S1_i) )
						+ Sigma_aa_i )
						}

				n=max(nobs)
				dpl = 	ka^2 +  max(ka * n, ka^2)+2*ka
				cr<-.C('test_pfda_bin_cond_aa',
					aa_hat=double(ka**2),
					w_i= as.double(wi), 
					ww_i= as.double(wwi), 
					n = as.integer(nobs[obsnum]),
					M = as.integer(M),
					S1_i = as.double(S1_i),
					phi = as.double(phi),
					Sigma_aa_i = as.double(Sigma_aa_i),
					ka = as.integer(ka),	
					offset = as.integer(sum(nobs[0:(obsnum-1)])),
					0L,
					# pfda:::c_debug(355),
					double(dpl)
					)

				cbind(as.double(aa_hat),cr$aa_hat)
					
				if(!isTRUE(eqd<-all.equal(as.double(aa_hat),cr$aa_hat))){
					cat('Failed on attempt',simnum,"with observation",obsnum,'\n',eqd,'\n')
					return(invisible(FALSE))
				}
			}
		}
	cat("PASS\n")
	} else {
	cat("no",fname,"function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfda_bin_cond_bb<-function(){
	fname='pfda_bin_cond_bb'
	cat('UNIT TEST -',fname,'\n')
	if(is.loaded(fname)){
		cat("testing with random input...")
		Nsim<-100
		# p=print
		p=invisible
		
		for(simnum in 1:Nsim){
			N=trunc(runif(1,min=6,max=50))
			nobs=trunc(runif(N,min=2,max=8))
			M=sum(nobs)
			ka=trunc(runif(1,min=1,max=5))
			kb=trunc(runif(1,min=1,max=5))

			w = rnorm(M)
			ww = sapply(nobs,function(n)crossprod(matrix(rnorm(n**2),n,n)))
			# ww = array(rnorm(4**2*N),dim=c(4,4,N))
			# for(obsnum in 1:N)ww[,,obsnum]<-crossprod(ww[,,obsnum])
			n = nobs[1]
			S2 = matrix(rnorm(kb*N),kb,N)
			phi = matrix(rnorm(ka*M),M,ka)
			Sigma_ab = array(rnorm(ka*kb*N),dim=c(ka,kb,N))
			Sigma_bb = array(rnorm(kb**2*N),dim=c(kb,kb,N))
			for(obsnum in 1:N)Sigma_bb[,,obsnum]<-crossprod(Sigma_bb[,,obsnum])

			for(obsnum in 1:N){
				(Sigma_ab_i = matrix(Sigma_ab[,,obsnum],ka,kb))
				(Sigma_bb_i = matrix(Sigma_bb[,,obsnum],kb,kb))
				(phii = phi[(sum(nobs[0:(obsnum-1)])+1):(sum(nobs[0:obsnum])),])
				(wi = w[(sum(nobs[0:(obsnum-1)])+1):(sum(nobs[0:obsnum]))])
				(wwi = ww[[obsnum]])
				(S2_i = S2[,obsnum])

				bb_hat = 	{
						p(p(p(p(t(Sigma_ab_i) %*% p(t(phii) %*% wwi %*% phii) %*% Sigma_ab_i )
						+ t(Sigma_ab_i)%*%t(phii)%*%wi%*%t(S2_i) + S2_i%*%t(t(Sigma_ab_i)%*%t(phii)%*%wi) )
						+ tcrossprod(S2_i) )
						+ Sigma_bb_i )
						}

				n=max(nobs)
				dpl = ka^2 +  max(ka * n, ka * kb) + ka + kb
				cr<-.C('test_pfda_bin_cond_bb',
					bb_hat=double(kb**2),
					w_i= as.double(wi),
					ww_i= as.double(wwi),
					n = as.integer(nobs[obsnum]),
					M = as.integer(M),
					S2_i = as.double(S2_i),
					phi = as.double(phi),
					Sigma_ab_i = as.double(Sigma_ab_i),
					Sigma_bb_i = as.double(Sigma_bb_i),
					ka = as.integer(ka),
					kb = as.integer(kb),
					offset = as.integer(sum(nobs[0:(obsnum-1)])),
					0L,
					# pfda:::c_debug(357,358),
					double(dpl)
					)

				if(!isTRUE(eqd<-all.equal(as.double(bb_hat),cr$bb_hat))){
					cat('Failed on attempt',simnum,"with observation",obsnum,'\n',eqd,'\n')
					return(invisible(FALSE))
				}
			}
		}
	cat("PASS\n")
	} else {
	cat("no",fname,"function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfda_bin_cond_ab<-function(){
	fname='pfda_bin_cond_ab'
	cat('UNIT TEST -',fname,'\n')
	if(is.loaded(fname)){
		cat("testing with random input...")
		Nsim<-100
		# p=print
		p=invisible
		
		for(simnum in 1:Nsim){
			N=trunc(runif(1,min=6,max=50))
			nobs=trunc(runif(N,min=2,max=8))
			M=sum(nobs)
			ka=trunc(runif(1,min=1,max=5))
			kb=trunc(runif(1,min=1,max=5))

			w = rnorm(M)
			ww = sapply(nobs,function(n)crossprod(matrix(rnorm(n**2),n,n)))
			# ww = array(rnorm(4**2*N),dim=c(4,4,N))
			# for(obsnum in 1:N)ww[,,obsnum]<-crossprod(ww[,,obsnum])
			n = nobs[1]
			S1 = matrix(rnorm(ka*N),ka,N)
			S2 = matrix(rnorm(kb*N),kb,N)
			phi = matrix(rnorm(ka*M),M,ka)
			Sigma_ab = array(rnorm(ka*kb*N),dim=c(ka,kb,N))
			Sigma_aa = array(rnorm(ka**2*N),dim=c(ka,ka,N))
			for(obsnum in 1:N)Sigma_aa[,,obsnum]<-crossprod(Sigma_aa[,,obsnum])

			for(obsnum in 1:N){
				(Sigma_ab_i = matrix(Sigma_ab[,,obsnum],ka,kb))
				(Sigma_aa_i = matrix(Sigma_aa[,,obsnum],ka,ka))
				(phii = phi[(sum(nobs[0:(obsnum-1)])+1):(sum(nobs[0:obsnum])),])
				(wi = w[(sum(nobs[0:(obsnum-1)])+1):(sum(nobs[0:obsnum]))])
				(wwi = ww[[obsnum]])
				(S1_i = S1[,obsnum])
				(S2_i = S2[,obsnum])

				ab_hat = 	{
						p(p(p(p(p(Sigma_aa_i %*% p(t(phii) %*% wwi %*% phii) %*% Sigma_ab_i )
						+ S1_i%*%t(t(Sigma_ab_i)%*%t(phii)%*%wi) )
						+ Sigma_aa_i%*%t(phii)%*%wi%*%t(S2_i) )
						+ tcrossprod(S1_i,S2_i) )
						+ Sigma_ab_i )
						}

				n=max(nobs)
				dpl = ka^2 +  max(ka * n, ka^2 + 2 * ka + kb)
				cr<-.C('test_pfda_bin_cond_ab',
					ab_hat=double(ka*kb),
					w_i= as.double(wi),
					ww_i= as.double(wwi),
					n = as.integer(nobs[obsnum]),
					M = as.integer(M),
					S1_i = as.double(S1_i),
					S2_i = as.double(S2_i),
					phi = as.double(phi),
					Sigma_aa_i = as.double(Sigma_aa_i),
					Sigma_ab_i = as.double(Sigma_ab_i),
					ka = as.integer(ka),
					kb = as.integer(kb),
					offset = as.integer(sum(nobs[0:(obsnum-1)])),
					0L,
					# pfda:::c_debug(359),
					double(dpl)
					)

				if(!isTRUE(eqd<-all.equal(as.double(ab_hat),cr$ab_hat))){
					cat('Failed on attempt',simnum,"with observation",obsnum,'\n',eqd,'\n')
					return(invisible(FALSE))
				}
			}
		}
	cat("PASS\n")
	} else {
	cat("no",fname,"function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
}
{ # Binary Single
UT_pfda_bin_single_e<-function(){
	fname='pfda_bin_single_e'
	cat('UNIT TEST -',fname,'\n')
	if(is.loaded(fname)){
		cat("testing with random input...")
		Nsim<-100
		# p=print
		p=invisible
		
		for(simnum in 1:Nsim){
			p     = trunc(runif(1,min=6,max=20))
			N     = trunc(runif(1,min=20,max=100))
			nobs  = trunc(runif(N,min=3,max=10))
			M     = sum(nobs)
			w     = rnorm(M)
			ww    = vector('list',N) ; wwc=NULL ; for(i in 1:N){ ww[[i]] = crossprod(matrix(rnorm(nobs[i]^2),nobs[i]));wwc = c(wwc,as.double(ww[[i]]))}
			ka    = trunc(runif(1,1,p))
			Bl    = lapply(nobs,function(n) matrix(rnorm(p*n),n,p))
			B     = NULL ; for(i in 1:N)B<-rbind(B, Bl[[i]])
			lf    = rchisq(1,1)
			ISD   = crossprod(matrix(rnorm(p**2),p,p))
			tm    = rnorm(p)
			tf    = matrix(rnorm(p*ka),p,ka)
			seps  = rexp(1)
			alpha = matrix(rnorm(N*ka), N, ka)
			Sigma_aa = array(rnorm(ka*ka*N),dim=c(ka,ka,N)) ; for(i in 1:N) Sigma_aa[,,i]<- crossprod(Sigma_aa[,,i])
			Da    = sort(rexp(ka),T)
			Dam = diag(Da,ka)

			id=rep(1:N,nobs)

			unlist(lapply(Bl,function(x){
				crossprod(x)->y
				y[lower.tri(y)]<-0
				y }))->btb
			dim(btb)<-c(p,p,N)

			Saanew<-array(0,dim=c(ka,ka,N))
			alphanew<-matrix(0,N,ka)
			aahat<-array(0,dim=c(ka,ka,N))
			subject=1
			for(subject in 1:N) {
				wi = w[id==subject]
				phii<-Bl[[subject]]%*%tf
				Saanew[,,subject]<-S<- solve(crossprod(phii)+solve(Dam))
				alphanew[subject,]<-S%*%crossprod(phii,w[id==subject]-Bl[[subject]]%*%tm)
				
				s1 = -Saanew[,,subject]%*%crossprod(phii,Bl[[subject]]%*%tm)
				aahat[,,subject] = { S + tcrossprod(s1) +
								 + tcrossprod(S%*%crossprod(phii,wi),s1)+
								 + tcrossprod(s1,S%*%crossprod(phii,wi))+
								 + tcrossprod(tcrossprod(S,phii)%*%ww[[subject]],tcrossprod(S,phii)) }
			}	
			dp = M + M*ka + N*ka + 2* ka^2 + 2 * ka + ka*max(nobs)
			ip = ka

			cr<- .C( 'pfda_bin_single_e',
				aa_hats = double(ka^2*N),
				Alpha = as.double(alpha),
				Sigma_aa = as.double(Sigma_aa),
				w = as.double(w),
				ww = as.double(wwc),
				nobs = as.integer(nobs),
				M = as.integer(M),
				N = as.integer(N),
				k = as.integer(ka),
				B = as.double(B),
				p = as.integer(p),
				tm = as.double(tm),
				tf = as.double(tf),
				Da = as.double(Da),
				dl = 0L,
				dp = double (dp),
				ip = integer(ip)
				)

			if(!isTRUE(all.equal(alphanew,matrix(cr$Alpha,N,ka)))){
				cat('Failed on attempt',simnum,'for values for alpha\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(as.double(Saanew),cr$Sigma_aa))){
				cat('Failed on attempt',simnum,'for values for Sigma_aa\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(as.double(aahat),cr$aa_hats))){
				cat('Failed on attempt',simnum,'for values for aa_hats\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
	} else {
	cat("no",fname,"function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfda_bin_single_generate_w_parms1<-function(){
	fname='pfda_bin_single_generate_w_parms1'
	cat('UNIT TEST -',fname,'\n')
	if(is.loaded(fname)){
		cat("testing with random input...")
		Nsim<-100
		# p=print
		p=invisible
		
		for(simnum in 1:Nsim){
			ni  = trunc(runif(1,4,15))
			p   = trunc(runif(1,7,15))
			k   = trunc(runif(1,1,4))
			wij = rnorm(ni)
			Bij = matrix(rnorm(p*ni),ni,p)
			tf  = matrix(rnorm(p*k),p,k)
			tm  = rnorm(p)
			Si  = crossprod(matrix(rnorm(k^2),k,k))
			Da  = sort(rexp(k),T)

			Ss <- solve(crossprod(Bij%*%tf)+solve(diag(Da,k)))
			mu <- Ss%*%t(Bij%*%tf)%*%(wij-Bij%*%tm)

			C<-.C('pfda_bin_single_generate_w_parms1',
				mu    = double(k),
				sigma = double(k^2),
				rw    = as.double(wij-Bij%*%tm),
				Da    = as.double(Da),
				phii  = as.double(Bij%*%tf),
				k     = as.integer(k),
				ni    = as.integer(ni),
				dl    = 0L,
				dp=double(2*k^2+k), ip=integer(k))

			if(!isTRUE(all.equal(C$mu,as.double(mu)))){
				cat('Failed on attempt',simnum,'for values for mu\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(C$sigma,as.double(Ss)))){
				cat('Failed on attempt',simnum,'for values for sigma\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
	} else {
	cat("no",fname,"function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfda_bin_single_generate_w_parms2<-function(){
	fname='pfda_bin_single_generate_w_parms2'
	cat('UNIT TEST -',fname,'\n')
	if(is.loaded(fname)){
		cat("testing with random input...")
		Nsim<-100
		# p=print
		p=invisible
		
		for(simnum in 1:Nsim){
			ni  = trunc(runif(1,4,15))
			p   = trunc(runif(1,7,15))
			k   = trunc(runif(1,1,4))
			wi = rnorm(ni)
			Bi  = matrix(rnorm(p*ni),ni,p)
			tf  = matrix(rnorm(p*k),p,k)
			tm  = rnorm(p)
			Si  = crossprod(matrix(rnorm(k^2),k,k))
			Da  = sort(rexp(k),T)

			j = sample(seq_len(ni),1)

			Bij <- Bi[-j,,drop=F]
			wij <- wi[-j]

			Ss <- solve(crossprod(Bij%*%tf)+solve(diag(Da,k)))
			mu <- Ss%*%t(Bij%*%tf)%*%(wij-Bij%*%tm)

			a <- t(Bi[j,])%*%tm + t(Bi[j,])%*%tf%*%mu
			s <- 1+t(Bi[j,])%*%tf%*%Ss%*%t(tf)%*%Bi[j,]
				
			
			C<-.C('pfda_bin_single_generate_w_parms2',
				a        = double(1),
				s        = double(1),
				mu       = as.double(mu),
				sigma    = as.double(Ss),
				B        = as.double(Bi[j,]),
				M        = as.integer(1),
				p        = as.integer(NCOL(Bi)),
				tm       = as.double(tm),
				tf       = as.double(tf),
				k        = as.integer(k),
				dl       = 0L, #pfda:::Cdebug(1:500),#0L,
				dp       = double(p+p**2+k*p))
			
			
			if(!isTRUE(all.equal(as.double(a),C$a))){
				cat('Failed on attempt',simnum,'for values for a\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(as.double(s),C$s))){
				cat('Failed on attempt',simnum,'for values for s\n')
				return(invisible(FALSE))
			}
		}
		cat("PASS\n")
	} else {
	cat("no",fname,"function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
X_pfda_bin_s_gen_w<-expression({
	ni  = trunc(runif(1,4,15))
	p   = trunc(runif(1,7,15))
	k   = trunc(runif(1,1,4))
	wi  = rnorm(ni)
	Bi  = matrix(rnorm(p*ni),ni,p)
	tf  = matrix(rnorm(p*k),p,k)
	tm  = rnorm(p)
	Da  = sort(rexp(k),T)
	kr=10
	# compute y
	yi = wi>0
	j = sample(seq_len(ni),1)
	# drop the Bij and wij
	Bij <- Bi[-j,,drop=F]
	wij <- wi[-j]
	# R computations
	saved.seed<-.Random.seed
	R<-pfda:::.single.b.w.genw(yi,wi,Bi,tm,tf,Da,kr,j)
	# Compute C arguments
	saved.seed->>.Random.seed
	dpl = ni* k + ni+ k^2 + k + p + p^2 + k*p + 3*k
	C<-.C('pfda_bin_s_gen_w',
		w_sim =   double(kr),
		Nsim  = as.integer(kr),
		Yi    = as.integer(yi), 
		RWi   = as.double(wi-Bi%*%tm),
		ni    = as.integer(ni),
		Bi    = as.double(Bi),
		M     = as.integer(NROW(Bi)), 
		p     = as.integer(p), 
		tm    = as.double(tm),
		tf    = as.double(tf),
		k     = as.integer(k),
		Da    = as.double(Da),
		j     = as.integer(j-1),
		dl    = 0L,
		dp=double(dpl), ip=integer(k))
	gc()
	# Test equality
	stopifnot(all.equal(C$w_sim,R))
})
UT_pfda_bin_s_gen_w<-UT_generate(X_pfda_bin_s_gen_w)

UT_pfda_bin_single_approximate_moments_forobs<-function(){
	fname='pfda_bin_single_approximate_moments_forobs'
	cat('UNIT TEST -',fname,'\n')
	if(is.loaded(fname)){
		cat("testing with random input...")
		Nsim<-100
		# p=print
		p=invisible

		if(!exists('.Random.seed'))runif(1)
		save.seed<-.Random.seed
		for(simnum in 1:Nsim){
			ni  = trunc(runif(1,4,15))
			p   = trunc(runif(1,7,15))
			k   = trunc(runif(1,1,4))
			wi  = rnorm(ni)
			Bi  = matrix(rnorm(p*ni),ni,p)
			tf  = matrix(rnorm(p*k),p,k)
			tm  = rnorm(p)
			Da  = sort(rexp(k),T)

			yi = wi>0

			kr = sample(10:100,1)

			save.seed<-.Random.seed
			seed<-trunc(runif(1,1e6,1e7))
			set.seed(seed)
			R<-pfda:::.single.b.w.1(yi,Bi,wi,tm,tf,Da,kr)

			set.seed(seed)
			pool = double(kr*ni + kr + (ni* k + ni+ k^2 + k + p + p^2 + k*p + 3*k)+100)
			C<-.C('pfda_bin_single_approximate_moments_forobs',
				wi  = as.double(wi),
				wwi = double(ni^2),
				#as.double(wi - Bi%*%tm),
				as.integer(yi),
				as.integer(ni),
				as.integer(NROW(Bi)),
				as.integer(p),
				as.integer(k),
				as.double(tm),
				as.double(Bi),
				as.double(tf),
				as.double(Da),
				as.integer(kr),
				as.double(1),
				debug = 0L, 
				dp=as.double(pool),ip=integer(k+10))
			gc()
			
			(cww<-matrix(C[[2]],ni,ni))
			(rww<-R[[2]])
			idx<-upper.tri(cww)
			
			if(!isTRUE(all.equal(C[[1]],R[[1]]))){
				cat('Failed on attempt',simnum,'for w\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(cww[idx],rww[idx]))){
				cat('Failed on attempt',simnum,'for ww\n')
				return(invisible(FALSE))
			}
		}
		.Random.seed<-save.seed
		cat("PASS\n")
	} else {
	cat("no",fname,"function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
UT_pfda_bin_single_approximate_moments<-function(){
	fname='pfda_bin_single_approximate_moments'
	cat('UNIT TEST -',fname,'\n')
	if(is.loaded(fname)){
		cat("testing with random input...")
		Nsim<-100
		# p=print
		p=invisible

		if(!exists('.Random.seed'))runif(1)
		save.seed<-.Random.seed
		for(simnum in 1:Nsim){
			N   = trunc(runif(1,4,15))
			no  = trunc(runif(N,4,15))
			M   = sum(no)
			p   = trunc(runif(1,7,15))
			k   = trunc(runif(1,1,4))
			w  = rnorm(M)
			ww = lapply(no,function(n)crossprod(matrix(rnorm(n*n),n,n)))
			B   = matrix(rnorm(p*M),M,p)
			tf  = matrix(rnorm(p*k),p,k)
			tm  = rnorm(p)
			Da  = sort(rexp(k),T)
			obs = factor(rep(letters[1:N],no))

			y = w>0

			kr = 10
			weight = 1 #runif(1)

			save.seed<-.Random.seed
			seed<-trunc(runif(1,1e6,1e7))

			set.seed(seed)
			R<-pfda:::.single.b.w(y,B,obs,w,ww,tm,tf,Da,weight,kr)

			ni=max(no)
			pool = double(kr*ni + kr + (ni* k + ni+ k^2 + k + p + p^2 + k*p + 3*k))
			set.seed(seed)
			C<-.C('pfda_bin_single_approximate_moments',
				w=as.double(w),
				ww=as.double(unlist(ww)),
				as.integer(y),
				as.integer(no),
				as.integer(M),
				as.integer(N),
				as.double(B),
				as.integer(p),
				as.integer(k),
				as.double(tm),
				as.double(tf),
				as.double(Da),
				as.double(weight),
				as.integer(kr),
				0L, #pfda:::Cdebug(412:414),
				dp=as.double(pool), 
				integer(k),
				DUP=TRUE)
			gc()
			
			if(!isTRUE(all.equal(C$w,R$w))){
				cat('Failed on attempt',simnum,'for w\n')
				return(invisible(FALSE))
			}
			if(!isTRUE(all.equal(C$ww,unlist(R$ww)))){
				cat('Failed on attempt',simnum,'for ww\n')
				return(invisible(FALSE))
			}
		}
		.Random.seed<-save.seed
		cat("PASS\n")
	} else {
	cat("no",fname,"function found.\n")
	cat("Failed\n")
		invisible(FALSE)
	}
}
}
{ # Dual binary/continuous
X_dual_bc_i<-expression({
	eval(get('sim.dual.bc',envir=parent.env(environment())))
	ldp = 2*p^2 + 2*M + M*max(ka,kb)+ p*N + 8*p
	lip = 6*p
	C<- { .C( 'dual_bc_i',
			tm=double(p),
			tn=double(p),
			tf=double(p*ka),
			tg=double(p*kb),
			alpha=double(ka*N),
			beta=double(kb*N),
			lambda=double(ka*kb),
			Da=double(ka),
			Db=double(kb),
			s.xi=double(1),
			aa=double(ka*ka*N),
			bb=double(kb*kb*N),
		as.integer(y),
		as.double(z),
		as.double(B),
		as.integer(nobs),
		as.double(btb),
		as.integer(N),
		as.integer(M),
		as.integer(ka),
		as.integer(kb),
		as.integer(p),
		as.double(1e-4),
		Cdebug(),
		dp=double(ldp), ip=integer(lip))
	}
	gc()	
	R<-.dual.bc.i(y,z,B,subject,ka,kb,1e-4)
	{
		if(!isTRUE(all.equal(as.vector(R$tm   )     ,C$tm   )))stop("failed on tm"   )
		if(!isTRUE(all.equal(as.vector(R$tf   )     ,C$tf   )))stop("failed on tf"   )
		if(!isTRUE(all.equal(as.vector(R$alpha)     ,C$alpha)))stop("failed on alpha")
		if(!isTRUE(all.equal(as.vector(R$Da   )     ,C$Da   )))stop("failed on Da"   )
		if(!isTRUE(all.equal(as.vector(R$tn   )     ,C$tn   )))stop("failed on tn"   )
		if(!isTRUE(all.equal(as.vector(R$tg   )     ,C$tg   )))stop("failed on tg"   )
		if(!isTRUE(all.equal(as.vector(R$beta )     ,C$beta )))stop("failed on beta" )
		if(!isTRUE(all.equal(as.vector(R$Db   )     ,C$Db   )))stop("failed on Db"   )
		if(!isTRUE(all.equal(as.vector(R$s.xi )     ,C$s.xi )))stop("failed on s.xi" )
	}
})
UT_dual_bc_i<-UT_generate(X_dual_bc_i)

X_dual_bc_1a<-expression({
	eval(get('sim.dual.bc',envir=parent.env(environment())))
	Saa = crossprod(matrix(rnorm(ka**2),ka,ka))
	Sab = matrix(rnorm(ka*kb),ka,kb)
	
	ldp = ka+kb
	dp = double(ldp)
	C<-.C('dual_bc_1a',
		mu  =    double(ka),
		Rw  = as.double(w-B%*%tm),
		Rz  = as.double(z-B%*%tn),
		phi = as.double(B%*%tf),
		psi = as.double(B%*%tg),
		s.xi= as.double(s.xi),
		Sa  = as.double(Saa),
		Sab = as.double(Sab),
		M   = as.integer(M),
		ni  = as.integer(nobs),
		ka  = as.integer(ka),
		kb  = as.integer(kb),
		0L,
		dp,
		DUP=TRUE)
	gc()
	
	ix=1:nobs[1]
	zi =z[ix]
	Bi = B[ix,]
	wi = w[ix]
	R<-.dual.bc.1a(zi,Bi,wi,tm,tn,tf,tg,s.xi,Saa,Sab)
	if(!isTRUE(all.equal(C$mu,as.vector(R))))stop("failed.")
})
UT_dual_bc_1a<-UT_generate(X_dual_bc_1a)

X_dual_bc_1b<-expression({
	eval(get('sim.dual.bc',envir=parent.env(environment())))
	Sbb = crossprod(matrix(rnorm(kb**2),kb,kb))
	Sab = matrix(rnorm(ka*kb),ka,kb)
	
	ldp = ka+kb
	dp = double(ldp)
	C<-.C('dual_bc_1b',
		mu  =    double(kb),
		Rw  = as.double(w-B%*%tm),
		Rz  = as.double(z-B%*%tn),
		phi = as.double(B%*%tf),
		psi = as.double(B%*%tg),
		s.xi= as.double(s.xi),
		Sab = as.double(Sab),
		Sbb = as.double(Sbb),
		M   = as.integer(M),
		ni  = as.integer(nobs),
		ka  = as.integer(ka),
		kb  = as.integer(kb),
		0L,
		dp,
		DUP=FALSE)
	gc()
	
	ix=1:nobs[1]
	zi =z[ix]
	Bi = B[ix,]
	wi = w[ix]
	R<-.dual.bc.1b(zi,Bi,wi,tm,tn,tf,tg,s.xi,Sab,Sbb)
	if(!isTRUE(all.equal(C$mu,as.vector(R))))stop("failed.")
})
UT_dual_bc_1b<-UT_generate(X_dual_bc_1b)

X_dual_bc_1cde<-expression({
	eval(get('sim.dual.bc',envir=parent.env(environment())))
Saa = crossprod(matrix(rnorm(ka**2),ka,ka))
Sbb = crossprod(matrix(rnorm(kb**2),kb,kb))
Sab = matrix(rnorm(ka*kb),ka,kb)

i=1
ni=nobs[1]
ldp = 3*ka +3*kb + ka*ni + kb*ni + ni*max(ka,kb)
dp = double(ldp)
C<-.C("dual_bc_1cde",
	aa =    double(ka**2),
	ab =    double(ka*kb),
	bb =    double(kb**2),
	Rz = as.double(z-B%*%tn),
	wi = as.double(w),
	rho= as.double(-B%*%tm),
	phi= as.double(B%*%tf),
	psi= as.double(B%*%tg),
	wwi= as.double(ww[[i]]),
	tf = as.double(tf),
	tg = as.double(tg),
	sxi= as.double(s.xi),
	Saa= as.double(Saa),
	Sab= as.double(Sab),
	Sbb= as.double(Sbb),
	M  = as.integer(M),
	ni = as.integer(nobs),
	ka = as.integer(ka),
	kb = as.integer(kb),
	0L,  as.double(dp),
	DUP=F)

i=1
ix = subject==i
zi<-z[ix]
Bi<-B[ix,]
wi<-w[ix]
wwi<-ww[[i]]
R<-.dual.bc.1cde(z[ix],B[ix,],w[ix],ww[[i]],tm,tn,tf,tg,s.xi,Saa,Sab,Sbb)


uta<-upper.tri(R$aa,T)
if(!isTRUE(all.equal(as.vector(R$aa)[uta],C$aa[uta])))return('failed')
utb<-upper.tri(R$bb,T)
if(!isTRUE(all.equal(as.vector(R$bb)[utb],C$bb[utb])))return('failed')
if(!isTRUE(all.equal(as.vector(R$ab),C$ab)))return('failed')
})
UT_dual_bc_1cde<-UT_generate(X_dual_bc_1cde)

X_dual_bc_1<-expression({
	eval(get('sim.dual.bc',envir=parent.env(environment())))
	R<-.dual.bc.1(z,B,subject,w,ww,tm,tn,tf,tg,lambda,Da,Db,s.xi)
	C<-{ # C related computations for dual_bc_1
	ldp = M*(ka + kb + 3) + max(7 * max(ka,kb)^2 , 3*(ka+kb) + (ka + kb+ max(ka,kb))*max(nobs))
	lip = max(ka,kb)
	dp = double(ldp)
	ip = integer(lip)
	.C("dual_bc_1",
		alpha =matrix(0,N,ka),
		beta  =matrix(0,N,kb),
		aa    =double(N*ka*ka),
		ab    =double(N*ka*kb),
		bb    =double(N*kb*kb),
		Saa   =double(N*ka*ka),
		Sab   =double(N*ka*kb),
		Sbb   =double(N*kb*kb),
		z     =as.double(z),
		B     =as.double(B),
		w     =as.double(w),
		ww    =as.double(unlist(ww)),
		tm    =as.double(tm),
		tn    =as.double(tn),
		tf    =as.double(tf),
		tg    =as.double(tg),
		lambda=as.double(lambda),
		Da    =as.double(Da),
		Db    =as.double(Db),
		sxi   =as.double(s.xi),
		nobs  =as.integer(nobs),
		N     =as.integer(N),
		M     =as.integer(M),
		ka    =as.integer(ka),
		kb    =as.integer(kb),
		p     =as.integer(p), 
		Cdebug(200), 
		dp, ip,
		DUP=FALSE)
	}
	stopifnot( all.equal(R$alpha,C$alpha))
	stopifnot( all.equal(R$beta,C$beta))
	stopifnot( all.equal(as.vector(unlist(R$Saa)),C$Saa))
	stopifnot( all.equal(as.vector(unlist(R$Sab)),C$Sab))
	stopifnot( all.equal(as.vector(unlist(R$Sbb)),C$Sbb))
	stopifnot( all.equal(as.vector(R$ab),C$ab))
})
UT_dual_bc_1<-UT_generate(X_dual_bc_1)

X_dual_bc_2<-expression({
	eval(get('sim.dual.bc',envir=parent.env(environment())))
Sbb<-vector('list',N)
for(i in seq(N)){ Sbb[[i]]<-crossprod(matrix(rnorm(kb*kb),kb,kb))}

R<-.dual.bc.2(z,B,subject,tn,tg,beta,Sbb)

dp =double(M + M*kb + 2*kb^2 + 10000)
C<-.C('dual_bc_2',s.xi=double(1),z,B,tn,tg,beta,unlist(Sbb),
	as.integer(nobs),as.integer(N),as.integer(M),as.integer(kb),as.integer(p),
	pfda:::Cdebug(123),dp)
stopifnot(all.equal(as.vector(R$s.xi),C$s.xi))
})
UT_dual_bc_2<-UT_generate(X_dual_bc_2)

X_dual_bc_3<-expression({
	eval(get('sim.dual.bc',envir=parent.env(environment())))
lm=runif(1,1,100)
ln=runif(1,1,100)

R<-.dual.bc.3(z,B,subject,w,tf,tg,alpha,beta,s.xi,lm,ln,K)

dp =double(p^2+ M+ M*max(ka,kb))
C<-{ .C('dual_bc_3', tm=double(p), tn=double(p), z, B, w, tf, tg, alpha, beta, s.xi, nobs, N, M, ka, kb, p, lm, ln, K, 0L, dp)}
gc()

all.equal(as.vector(R$tm),C$tm)
all.equal(as.vector(R$tn),C$tn)
})
UT_dual_bc_3<-UT_generate(X_dual_bc_3)

X_dual_bc_4<-expression({
runif(1);.saved.seed<-.Random.seed
	eval(get('sim.dual.bc',envir=parent.env(environment())))
lf = runif(1,1,100)
lg = runif(1,1,100)

aa<-vector('list',N)
bb<-vector('list',N)
for(i in seq(N)){ 
	aa[[i]]<-crossprod(matrix(rnorm(ka**2),ka,ka))
	bb[[i]]<-crossprod(matrix(rnorm(kb**2),kb,kb))
}
aa<-array(unlist(aa),c(ka,ka,N))
bb<-array(unlist(bb),c(kb,kb,N))

R<-.dual.bc.4(z,B,subject,w,tm,tn,tf,tg,alpha,beta,s.xi,aa,bb,lf,lg,K)

dp = double(M + max(ka,kb)^2 * N + p^2 + p )
ip = integer(p)
C<-.C('dual_bc_4', tf=tf, tg=tg, z, B, w, tm, tn, alpha, beta, s.xi, aa, bb, as.integer(nobs), as.integer(N), as.integer(M), as.integer(ka), as.integer(kb), as.integer(p), lf, lg, K, unlist(btb), 0L,dp, ip)

if(!isTRUE(all.equal(R$tf,C$tf*rep(sign(C$tf[1,]),each=p))))stop('tf failed')
if(!isTRUE(all.equal(R$tg,C$tg*rep(sign(C$tg[1,]),each=p))))stop('tg failed')
TRUE
})
UT_dual_bc_4<-UT_generate(X_dual_bc_4)

X_dual_bc_5<-expression({
	eval(get('sim.dual.parameters',envir=parent.env(environment())))


aa<-vector('list',N)
ab<-vector('list',N)
for(i in seq(N)){ 
	aa[[i]]<-tcrossprod(matrix(rnorm(ka*ka),ka,ka))
	ab[[i]]<-matrix(rnorm(ka*kb),ka,kb)
}
aa<-array(unlist(aa),c(ka,ka,N))
ab<-array(unlist(ab),c(ka,kb,N))

R<-.dual.cc.5(aa,ab)

dp = double(ka^2 + ka*kb + 10*max(ka,kb))
ip = integer(ka)
C<-.C('dual_bc_5', lambda=lambda, aa, ab, N, ka, kb, 0L, dp, ip)

if(!isTRUE(all.equal(R$lambda,C$lambda)))stop('failed')
})
UT_dual_bc_5<-UT_generate(X_dual_bc_5)

X_dual_bc_6<-expression({
	eval(get('sim.dual.bc',envir=parent.env(environment())))
	eval(get('sim.dual.aa_ab_bb',envir=parent.env(environment())))

R<-.dual.cc.6(tf,tg,alpha,beta,lambda,aa,bb)

dpl=4*p^2 + 2*p +  max(  ka, kb, 8 ) * p 
ipl=6*p
C<-.C('dual_bc_6', tf=tf, tg=tg, alpha=alpha, beta=beta, lambda=lambda, Da=Da, Db=Db, aa=aa, bb=bb, N=N, ka=ka, kb=kb, p=p, dl=0L, dp=double(dpl), ip=integer(ipl))

if(!isTRUE(all.equal(R$tf,C$tf)))stop("failed on tf")
if(!isTRUE(all.equal(as.vector(R$tg),as.vector(C$tg))))stop("failed on tg")
if(!isTRUE(all.equal(as.vector(R$alpha),as.vector(C$alpha))))stop("failed on alpha")
if(!isTRUE(all.equal(as.vector(R$beta),as.vector(C$beta))))stop("failed on beta")
if(!isTRUE(all.equal(as.vector(R$lambda),as.vector(C$lambda))))stop("failed on lambda")
if(!isTRUE(all.equal(as.vector(R$Da),as.vector(C$Da))))stop("failed on lambda")
if(!isTRUE(all.equal(as.vector(R$Db),as.vector(C$Db))))stop("failed on lambda")
})
UT_dual_bc_6<-UT_generate(X_dual_bc_6)

X_dual_bc_genw<-expression({
	eval(get('sim.dual.bc',envir=parent.env(environment())))

ni=nobs[1]
j=sample(seq(ni),1)
Rw = w-B%*%tm
Rz = z-B%*%tn
Nsim=10L
phi=B%*%tf
psi=B%*%tg
i=1
ix = subject ==i
	Bi  <- B[ix,]
	Bij<-Bi;Bij[j,]=0
phii = Bij%*%tf

saved.seed<-.Random.seed
dpl = 2 + 5*ka + kb + 10 * max(ka,kb)^2
ipl = max(ka,kb)
C<-.C('test_dual_bc_genw', wsim = double(Nsim), y=as.integer(y), Rz=Rz, Rw=Rw, rho = (B%*%tm), phi=phi, psi=psi, lambda=lambda, Da=Da, Db=Db,ni=ni, M=M, ka=ka, kb=kb, Nsim=Nsim,  p=p, j=as.integer(j-1), dl=0L, dp=double(dpl), ip=integer(ipl))

saved.seed->>.Random.seed
R <- .dual.bc.genw(j,Nsim,y[ix],w[ix],z[ix],B[ix,],tm,tn,tf,tg,Da,Db,lambda)

cbind(C=C$wsim, R)
stopifnot(all.equal(C$wsim,R))
})
UT_dual_bc_genw<-function()eval(X_dual_bc_genw)#UT_generate(X_dual_bc_genw)

X_dual_bc_w_1<-expression({
	eval(get('sim.dual.bc',envir=parent.env(environment())))
kr = 100L
i=1
ix = i ==subject
ni=nobs[i]
weight = 1.0

saved.seed<-.Random.seed
R<-.dual.bc.w.1(kr,y[ix],w[ix],z[ix],B[ix,],tm,tn,tf,tg,Da,Db,lambda)

saved.seed->.Random.seed
dpl = max(nobs)*kr + kr +   (2 + 5*ka + kb + 10 * max(ka,kb)^2)
ipl = max(ka,kb)
C <- .C('dual_bc_w_1', 
	w=double(ni), ww=double(ni^2),
	y=as.integer(y), Rw=(w-B%*%tm), Rz=(z-B%*%tn), 
	pi=(B%*%tm), phi=(B%*%tf), psi=(B%*%tg), 
	lambda, Da, Db, 
	ni, M, ka, kb, 
	weight, kr, 
	p, 0L, dp=double(dpl), ip=integer(ipl))

stopifnot(t.test(C$w,R$wi,paired=TRUE)$p.value>0.05/100)
stopifnot(t.test(C$ww,as.vector(R$wwi),paired=TRUE)$p.value>0.05/100)
})

}
{ # Calcium Model

X_dual_ca_resid<-expression({
	eval(get('sim.dual.ca',envir=parent.env(environment())))
pt=NCOL(Bt)
px=NCOL(Bx)

Ry = y-Z%*%tz -Bt%*%tt-Bx%*%tx
for(i in seq(nlevels(subject))){
	ix = i==subject
	Ry[ix]<-Ry[ix] - Bt[ix,]%*%tf%*%gamma[i,] 
	Ry[ix]<-Ry[ix] - Bx[ix,]%*%tg%*%delta[i,]
}

C<-.C('dual_ca_resid', Ry=double(M), y, Z, Bt, Bx, tz, tt, tx, tf, tg, gamma, delta, nobs, N, M, as.integer(kz), kd, kg, pt, px, 0L, double(max(kg,kd)*M))
stopifnot(all.equal(as.vector(Ry),as.vector(C$Ry)))
})
UT_dual_ca_resid<-UT_generate(X_dual_ca_resid)

X_u_orthogonalize<-expression({
a = sample(5:10,1)
b = sample(1:(a-1),1)
U = matrix(rnorm(a*b),a,b)
U = U %*% solve(chol(crossprod(U)))

zapsmall(crossprod(U))

v = rnorm(a)

R<-.u.orthogonalize(U,v)
C<-.C('u_orthogonalize', v,  U, a, b)[[1]]

stopifnot(all.equal(R,C))
})
UT_u_orthogonalize<-UT_generate(X_u_orthogonalize)

X_gen_symblock_solve<-expression({
a<-as.integer(sample(1:20,1))
c<-as.integer(sample(1:20,1))

A<-crossprod(matrix(rnorm(a**2),a,a))
B<-(matrix(rnorm(a*c),a,c))
C<-crossprod(matrix(rnorm(c**2),c,c))

ABC<-rbind(cbind(A,B),cbind(t(B),C))

X<-solve(ABC)
R<-.gen.symblock.solve(A,B,C)
Y<-.C('gen_symblock_solve', A, B,  C, a, c, dp=double(a^2 + c^2 + max(a,c)*10 + 2*a*c) , ip=integer(max(a,c)))
gc()

all.equal(R,Y[1:3],tolerance = max(a,c)*.Machine$double.eps ^ 0.5, check.attributes = FALSE)

# stopifnot(isTRUE(all.equal(R[[3]],X[-seq(a),-seq(a),drop=FALSE])))
# stopifnot(isTRUE(all.equal(R[[3]],Y[[3]])))

# stopifnot(isTRUE(all.equal(R[[2]],X[seq(a),-seq(a),drop=FALSE])))
# stopifnot(isTRUE(all.equal(R[[2]],Y[[2]])))

# stopifnot(isTRUE(all.equal(R[[1]],X[seq(a),seq(a),drop=FALSE])))
# stopifnot(isTRUE(all.equal(R[[1]],Y[[1]],)))

})
UT_gen_symblock_solve<-UT_generate(X_gen_symblock_solve)

X_dual_ca_i<-expression({
	eval(get('sim.dual.ca',envir=parent.env(environment())))
	eval(get('def.defaults',envir=parent.env(environment())))

R<-.dual.ca.i(y,Z,Bt,Bx,subject,kg,kd,min.v)
names(R)


pt=NCOL(Bt)
px=NCOL(Bx)
p=max(pt,px)
k=max(kd,kg)
dpl = pt*pt*N + px*px*N + 2*p^2 + M + M*k+ p*N + 8*p
ipl = kz + 6*p
dp = double(dpl)
ip = integer(ipl)
C<-.C('dual_ca_i',tz=tz, tt=tt, tx=tx, tf=tf, tg=tg, gamma=gamma, delta=delta, lambda=lambda, Dg=Dg, Dd=Dd, sigma=sigma, gg=double(kg*kg*N), gd=double(kg*kd*N), dd=double(kd*kd*N), y=y, z=Z, Bt=Bt, Bx=Bx, nobs=nobs, N=N, M=M, kz=kz, kg=kg, kd=kd, pt=pt, px=px, minV=min.v, dl=0L, dp=dp, ip=ip)
gc()

test.equal(tz)
test.equal(tt)
test.equal(tx)
test.equal(tf)
test.equal(tg)
test.equal(gamma)
test.equal(delta)
test.equal(Dg)
test.equal(Dd)
test.equal(lambda)
test.equal(sigma)
test.equal.sa(gg)
test.equal.sa(gd)
test.equal.sa(dd)
})
UT_dual_ca_i<-UT_generate(X_dual_ca_i)

X_dual_ca_E1<-expression({
	eval(get('sim.dual.ca',envir=parent.env(environment())))
i=1
ix = i==subject

O<-.gen.symblock.solve(
	diag(Dg,length(Dg)),
	tcrossprod(diag(Dg,length(Dg)),lambda),
	diag(Dd,length(Dd)))

phi=Bt%*%tf
psi=Bx%*%tg

R<-.dual.ca.E.1(y[ix],tt,tx,sigma,phi[ix,,drop=F],psi[ix,,drop=F],O)
{
gg <-matrix(0,kg,kg)
gd <-matrix(0,kg,kd)
dd <-matrix(0,kd,kd)
Sgg<-matrix(0,kg,kg)
Sgd<-matrix(0,kg,kd)
Sdd<-matrix(0,kd,kd)
}
C<-.C('dual_ca_E1', gamma=gamma, delta=delta, gg=gg, gd=gd, dd=dd, Sgg=Sgg, Sgd=Sgd, Sdd=Sdd, Ri=y, phi=phi,  psi=psi, sigma=sigma,  O1=O[[1]], O2=O[[2]], O3=O[[3]],ni=nobs[1], N=N, M=M, kg=kg, kd=kd, dp=double(kg^2 + kd^2 + max(kg,kd)*10 + 2*kg*kd), ip=integer(max(kg,kd)))

stopifnot(isTRUE(all.equal(R$Sigma,list(C$Sgg,C$Sgd,C$Sdd))))
stopifnot(isTRUE(all.equal(R$mu.gamma,C$gamma[1,])))
stopifnot(isTRUE(all.equal(R$mu.delta,C$delta[1,])))
})
UT_dual_ca_E1<-UT_generate(X_dual_ca_E1)

X_dual_ca_E<-expression({
	eval(get('sim.dual.ca',envir=parent.env(environment())))
R<-.dual.ca.E(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,lambda,Dg,Dd,sigma)

C<-.C('dual_ca_E',gamma=gamma, delta=delta, 
	gg=array(0,c(kg,kg,N)), gd=array(0,c(kg,kd,N)), dd=array(0,c(kd,kd,N)), Sgg=array(0,c(kg,kg,N)), Sgd=array(0,c(kg,kd,N)), Sdd=array(0,c(kd,kd,N)), 
	y=y, Z=Z, Bt=Bt,  Bx=Bx,  tz=tz,  tt=tt,  tx=tx,  phi=tf,  psi=tg, lambda=lambda,  Dg=Dg,  Dd=Dd,  sigma=sigma,  
	nobs=nobs, N=N, M=M, kz=kz, kg=kg, kd=kd, pt=ncol(Bt), px=ncol(Bx),
	dl=0L, dp=double(2^kg^2 + 2*kd^2 + kg*kd+ max(kg,kd)*10 + 2*kg*kd+(1+kg+kd)*M), ip=integer(max(kg,kd)))

O<-.gen.symblock.solve(
	diag(Dg,length(Dg)),
	tcrossprod(diag(Dg,length(Dg)),lambda),
	diag(Dd,length(Dd)))
		
stopifnot(all.equal(O[[1]] ,matrix(C$dp                            [seq(kg^2)],kg,kg)))
stopifnot(all.equal(O[[3]] ,matrix(C$dp[-seq(kg^2)]                [seq(kd^2)],kd,kd)))
stopifnot(all.equal(O[[2]] ,matrix(C$dp[-seq(kg^2+kd^2)]          [seq(kg*kd)],kg,kd)))
stopifnot(all.equal(Bt%*%tf,matrix(C$dp[-seq(kg^2+kd^2+kg*kd)]     [seq(M*kg)],M ,kg)))
stopifnot(all.equal(Bx%*%tg,matrix(C$dp[-seq(kg^2+kd^2+kg*kd+M*kg)][seq(M*kd)],M ,kd)))
	
test.equal(gamma)	
test.equal(delta)	
test.equal(gg)
test.equal(gd)	
test.equal(dd)	
test.equal(Sgg)
test.equal(Sgd)
test.equal(Sdd)
})
UT_dual_ca_E<-UT_generate(X_dual_ca_E)

X_dual_ca_unpenalized<-expression({
	eval(get('sim.dual.ca',envir=parent.env(environment())))
R<-.dual.ca.unpenalized(y,Z,Bt,Bx,subject,tt,tx,tf,tg,gamma,delta,sigma)
C<-.C('dual_ca_unpenalized', tz=tz, y=y, Z=Z, Bt=Bt, Bx=Bx, tt=tt, tx=tx, tf=tf, tg=tg, gamma=gamma, delta=delta, sigma=sigma, nobs=nobs, N=N, M=M, kz=kz, kg=kg, kd=kd, pt=ncol(Bt), px=ncol(Bx), dl=0L, 
dp=double(M + kz^2 + 10*kz + M*max(kg, kd)), ip=integer(kz))
gc()

test.equal(tz)
})
UT_dual_ca_unpenalized<-UT_generate(X_dual_ca_unpenalized)

X_dual_ca_penalized<-expression({
eval(get('sim.dual.ca',envir=parent.env(environment())))
pt=ncol(Bt)
px=ncol(Bx)
lt = rexp(1,pt)
lx = rexp(1,px)

R<-.dual.ca.penalized(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,gamma,delta,sigma,lt,lx,Kt,Kx)
C<-.C('dual_ca_penalized', tt=tt, tx=tx, y=y, Z=Z, Bt=Bt, Bx=Bx, tz=tz, tf=tf, tg=tg, gamma=gamma, delta=delta, sigma=sigma, nobs=nobs, N=N, M=M, kz=kz, kg=kg, kd=kd, pt=pt, px=px, lt=lt, lx=lx, Kt=Kt, Kx=Kx, dl=0L, dp=double(10*max(pt,px) + max(pt,px)^2 +M + M*max(kg,kd)), ip=integer(max(pt,px)))

test.equal(tt)

R<-.dual.ca.penalized(y,Z,Bt,Bx,subject,tz,C$tt,tx,tf,tg,gamma,delta,sigma,lt,lx,Kt,Kx)
test.equal(tx)
})
UT_dual_ca_penalized<-UT_generate(X_dual_ca_penalized)

X_dual_ca_princcomp<-expression({
{ eval(get('sim.dual.ca',envir=parent.env(environment())))
pt=ncol(Bt)
px=ncol(Bx)
lf = rexp(1,pt)
lg = rexp(1,px)
gg<-sim.aa(N,kg)
gd<-sim.ab(N,kg,kd)
dd<-sim.aa(N,kd)
}

R<-.dual.ca.princcomp(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,gamma,delta,sigma,gg,gd,dd, lf,lg,Kt,Kx)
C<-.C('dual_ca_princcomp', tf=tf, tg=tg, y, Z, Bt, Bx, tz, tt, tx, gamma, delta, sigma, gg, gd, dd, nobs, N, M, kz, kg, kd, pt, px, lf, lg, Kt, Kx, 0L, dp=double(M + max(pt,px)^2 + max(pt,px) + max( 10*max(pt,px) ,  M + N*(kg+kd) + M*max(kg,kd) )), ip=integer(max(pt,px)))
gc()

test.equal(tf)
test.equal(tg)
})
UT_dual_ca_princcomp<-UT_generate(X_dual_ca_princcomp)

X_dual_ca_variances<-expression({
{	eval(get('sim.dual.ca',envir=parent.env(environment())))
	eval(get('sim.dual.ca',envir=parent.env(environment())))
	pt=ncol(Bt)
	px=ncol(Bx)
	p= max(pt,px)
	k= max(kg,kd)
	gg<-sim.aa(N,kg)
	gd<-sim.ab(N,kg,kd)
	dd<-sim.aa(N,kd)
}
R<-.dual.ca.variances(y,Z,Bt,Bx,subject,tz,tt,tx,tf,tg,gamma,delta,sigma,gg,gd,dd)
C<-.C('dual_ca_variances', tf=tf, tg=tg, gamma=gamma, delta=delta, lambda=lambda, Dg=Dg, Dd=Dd, sigma=sigma, y=y, Z=Z, Bt=Bt, Bx=Bx, tz=tz, tt=tt, tx=tx, gg=gg, gd=gd, dd=dd, nobs=nobs, N=N, M=M, kz=kz, kg=kg, kd=kd, pt=pt, px=px, 0L, dp=double(2*kg^2 + 2*kd^2 + M + 10*k + 9*p + 3*p^2 + N*k + M*k), ip=integer(6*p))
gc()

test.equal(sigma)
test.equal(Dg)
test.equal(Dd)
test.equal(lambda)
test.equal(tf)
test.equal(tg)
test.equal(gamma)
test.equal(delta)
})
UT_dual_ca_variances<-UT_generate(X_dual_ca_variances)

}
{ # Random Number Generation
UT_pfda_gen_truncnorm<-function(n=100){
m=100
all(
	sapply(seq(-3,3,length=n),function(lower){
		seed<-.Random.seed
		(R<-pfda:::.rtruncnormlower(m,0,1,lower))
		seed->>.Random.seed
		(C<-.C("pfda_gen_truncnorm",x=double(m),as.integer(m),lower,0L)$x)
		all.equal(R,C)
	})
)
}
}
