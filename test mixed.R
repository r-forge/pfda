library(pfda)
d<-data.frame(
y=as.logical(rbinom(110,1,.5)),
z=rnorm(110,0,1),
x=rep(0:10,10),
id=factor(rep(1:10,each=11)))

pfda(y%&%z~x%|%id, data=d,k=1,df=c(2.1,2.1,2.1,2.1),driver="dual.mixed",
	control=pfdaControl(C.debug=pfda:::Cdebug(200:299)))->testbase


# Debugging
options(error=recover)
traceback()
library(debug)
mtrace(dual.bc)



