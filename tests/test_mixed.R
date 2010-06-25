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
sapply(x,function(x){ AIC(within(testbase,tm=tm*x)) })


# Can optimize?


# Debugging
options(error=recover)
traceback()
library(debug)
mtrace(dual.bc)
mtrace(AIC.pfda.dual.bc)

