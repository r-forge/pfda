# tests/single.binary.R
# testing and debugging code for the single binary 
library(pfda)
library(debug)
library(ggplot2)
qmatplot<-function(x,y,...){
	X<-rep(x, NCOL(y))
	Y<-as.vector(y)
	COL<-factor(rep(seq_len(NCOL(y)),each=NROW(y)))
	qplot(X,Y,group=COL,color=COL,...)
}

{ # Testing the base function
	# mtrace(single.b)
	set.seed(20100611)
	N<-10
	id<-factor(subject<-rep(1:N,each=11))
	x<-rep(-5:5,N)
	w<-pnorm(rnorm(n=length(x),mean=x+(subject-4.5)/2,sd=4))
	y<- w>.5
	# qplot(x,w,color=subject)+facet_wrap(~subject)
	m<-pfda(y~x%|%id,driver='single.binary', k=1, df=c(3,3))
}
{ # AIC
	{ # Can AIC even compute?
	# mtrace(AIC.pfda.single.b)
	# mtrace(.single.b.n2L)
	# mtrace(.single.b.n2L.1)
	# mtrace(.single.b.estimatePCS)
	# mtrace(.single.b.updatePCS)
	AIC(m)
	}
	{ # Is AIC smooth?
		dfm<-seq(2.1,5,length=50)
		results<-lapply(dfm,function(dfm){(pfda(y~x%|%id,driver='single.binary', k=1, df=c(dfm,3)))})
		aics<-sapply(results,function(object){AIC(updatelist(m,object['Da']))})
		qplot(dfm, aics, geom='line', xlab=expression('df'[mu]),ylab="AIC")
		matplot(dfm,t(sapply(results,'[[','tm')),type='l', ylab='tm')
		matplot(dfm,t(sapply(results,'[[','tf')),type='l',ylab='tf')
		matplot(dfm,t(sapply(results,'[[','alpha')),type='l',ylab='alpha')
		plot(dfm,t(sapply(results,'[[','Da')),type='l',ylab=expression('D'[alpha]))
		Das <-  sapply(results,'[[','Da')
		plot(Das, aics)
		
		ealpha<-sapply(results, function(m){ with(m,pfda:::.single.b.estimatePCS(y, Bt, subject, tm, tf, Da))})
		matplot(dfm,t(ealpha), type='l')
		plot(t(sapply(results,'[[','alpha')), t(ealpha), col=rep(1:10,each=50),pch='.')
		
		n2L<-sapply(results, function(m){ with(m,pfda:::.single.b.n2L(y, subject, Bt, tm, tf, Da))})
		plot(dfm,n2L,type='l')
		plot(n2L,aics)
		
		n2L.bySubject<-sapply(results,function(x)with(x,pfda:::.single.b.n2L.bySubject(y, subject, Bt, tm, tf, Da, alpha)))
		qmatplot(dfm, t(n2L.bySubject),geom='line')+ylim(c(0,60))
		
		n2L.parts<-function(si){ sapply(results,function(x)with(x,{
			#si=10
			ix=subject==si
			pfda:::.single.b.n2L.parts(y[ix], Bt[ix,], tm, tf, alpha[si,], Da)
		}))
		}
		#qplot(matrix(dfm,ncol=NROW(n2L.1.parts)), t(n2L.1.parts),geom='line',group=factor(rep(1:NCOL(n2L.1.parts),each=)
		sapply(seq_len(nlevels(subject)),function(i){
			parts<-t(n2L.parts(i))
			qmatplot(dfm, cbind(parts,apply(parts,1,sum)),geom='line')
		})
		
		#this finally shows smoothness in the AIC.
		# the non smoothness above shows the variation in the estimation.
		infl<-seq(.1, 2.5, length=100)
		aic.Da<-sapply(infl, function(infl){AIC(updatelist(m,list(Da=m$Da*infl)))})
		qplot(infl,aic.Da,geom=c('point','line'))
		
		aic.tm<-sapply(infl, function(infl){AIC(updatelist(m,list(tm=m$tm*infl)))})
		qplot(infl,aic.tm,geom=c('point','line'))
		
		aic.tf<-sapply(infl, function(infl){AIC(updatelist(m,list(tf=m$tf*infl)))})
		qplot(infl,aic.tf,geom=c('point','line'))
	}
}



mtrace.off()
