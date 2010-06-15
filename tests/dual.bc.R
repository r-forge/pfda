# dual.bc.R
# tests dual.bd code
rm(list=ls())
library(pfda)
library(debug)
library(ggplot2)

 # test data
	set.seed(20100615)
	N<-10
	id<-factor(subject<-rep(1:N,each=11))
	x<-rep(-5:5,N)
	w<-pnorm(rnorm(n=length(x),mean=x+(subject-4.5)/2,sd=4))
	y<- w>.5
	sxi<-12^2
	beta=rnorm(N, sd=50)
	z<- exp(x)+rnorm(length(x),sd=sqrt(sxi))+as.vector(sapply(beta,'*',(-5:5)))
	qplot(x,z,color=subject, geom=c('point','line'))+facet_wrap(~subject)
	
	m<-pfda(y%&%z~x%|%id,driver='dual.mixed', k=c(1,1), df=c(3,3,3,3), control=pfdaControl(C.debug=pfda:::Cdebug(200)))




