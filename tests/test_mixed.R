rm(list=ls())
library(pfda)
set.seed(12345)
with(new.env(parent=getNamespace('pfda')),{
	eval(sim.dual.bc)
	d<<-data.frame(y=y,z=z,x=times,id=subject)
})

# Basic
try({pfda(y~x%|%id, data=d,k=c(1,1),df=c(2.1,2.1,2.1,2.1),driver="single.binary")->testy})
try({pfda(z~x%|%id, data=d,k=c(1,1),df=c(2.1,2.1,2.1,2.1),driver="single.continuous")->testz})
try({pfda(y%&%z~x%|%id, data=d,k=c(1,1),df=c(2.1,2.1,2.1,2.1),driver="dual.mixed",
    control=pfdaControl()) -> testbase})
try({testbase})
# plot(testbase)	
	
# AIC
try({with(testbase, pfda:::.dual.bc.n2L(y, z, B, subject, tm, tn, tf, tg, lambda, Da, Db, sxi))})
try(AIC(testbase))

# check smoothness
# x<-seq(.5,3,length=250)
# library(ggplot2)
#for tm
# aic.tm<-sapply(x,function(x){ AIC( 	within.list(testbase,tm<-tm*x) 	) })
# qplot(x,aic.tm,xlab="Inflation factor",ylab='AIC',main=expression(theta[mu]),geom='line')
#for tn
# aic.tn<-sapply(x,function(x){ AIC( 	within.list(testbase,tn<-tn*x) 	) })
# qplot(x,aic.tn,xlab="Inflation factor",ylab='AIC',main=expression(theta[nu]),geom='line')

# Can optimize penalties?
# pfda(y%&%z~x%|%id, data=d,k=c(1,1),df=c(NA,2.1,2.1,2.1),driver="dual.mixed") -> testpen1
system.time(try({ pfda(y%&%z~x%|%id, data=d,k=c(1,1),df=NULL,driver="dual.mixed") -> testpen }))

#can optimize k?
system.time(try({ pfda(y~x%|%id,     data=d, k=NULL, penalties=testpen$penalties, driver="single.binary") -> testb }))
try(testb$k)
system.time(try({ pfda(z~x%|%id,     data=d, k=NULL, df=c(2.1,2.1), driver="single.continuous") -> testc }))
try(testc$k)
system.time(try({ pfda(y%&%z~x%|%id, data=d, k=c(3,1), df=c(2.1,2.1,2.1,2.1),driver="dual.mixed") -> testk }))

system.time(try({ pfda(y%&%z~x%|%id, data=d, k=NULL, df=c(2.1,2.1,2.1,2.1),driver="dual.mixed") -> testk }))

# can optimize both?
system.time(try({ pfda(y%&%z~x%|%id, data=d,k=NULL,df=NULL,driver="dual.mixed") -> testboth }))

# Debugging
if(F){
	options(error=recover)
	traceback()
	library(debug)
	mtrace(dual.bc)
	mtrace(.F.single.optimize.npc)
	mtrace(AIC.pfda.dual.bc)
}
