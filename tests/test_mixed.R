library(pfda)
d<-data.frame(
y=as.logical(rbinom(110,1,.5)),
z=rnorm(110,0,1),
x=rep(0:10,10),
id=factor(rep(1:10,each=11)))

# Basic
pfda(y%&%z~x%|%id, data=d,k=c(1,1),df=c(2.1,2.1,2.1,2.1),driver="dual.mixed",
		#	control=pfdaControl(C.debug=pfda:::Cdebug(200:299)))->testbase
    control=pfdaControl()) -> testbase
		
# AIC
with(testbase, pfda:::.dual.bc.n2L(y, z, B, subject, tm, tn, tf, tg, lambda, Da, Db, sxi))
AIC(testbase)

# check smoothness
x<-seq(.5,3,length=250)
library(ggplot2)
#for tm
aic.tm<-sapply(x,function(x){ AIC( 	within.list(testbase,tm<-tm*x) 	) })
qplot(x,aic.tm,xlab="Inflation factor",ylab='AIC',main=expression(theta[mu]),geom='line')
#for tn
aic.tn<-sapply(x,function(x){ AIC( 	within.list(testbase,tn<-tn*x) 	) })
qplot(x,aic.tn,xlab="Inflation factor",ylab='AIC',main=expression(theta[nu]),geom='line')

# Can optimize penalties?
pfda(y%&%z~x%|%id, data=d,k=c(1,1),df=c(NA,2.1,2.1,2.1),driver="dual.mixed") -> testpen1
system.time({ pfda(y%&%z~x%|%id, data=d,k=c(1,1),df=NULL,driver="dual.mixed") -> testpen })

#can optimize k?
system.time({ pfda(y%&%z~x%|%id, data=d,k=NULL,df=c(2.1,2.1,2.1,2.1),driver="dual.mixed") -> testk })
system.time({ pfda(y%&%z~x%|%id, data=d,k=NULL,df=c(2.1,2.1,2.1,2.1),driver="dual.mixed") -> testk })

# can optimize both?
system.time({ pfda(y%&%z~x%|%id, data=d,k=NULL,df=NULL,driver="dual.mixed") -> testboth })


# Debugging
options(error=recover)
traceback()
library(debug)
mtrace(dual.bc)
mtrace(AIC.pfda.dual.bc)

