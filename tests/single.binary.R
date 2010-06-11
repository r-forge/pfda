# tests/single.binary.R
# testing and debugging code for the single binary 
library(pfda)
library(debug)
library(ggplot2)

# Testing the base function
# mtrace(single.b)
set.seed(20100611)
N<-10
id<-factor(subject<-rep(1:N,each=11))
x<-rep(-5:5,N)
w<-pnorm(rnorm(n=length(x),mean=x+(subject-4.5)/2,sd=4))
y<- w>.5
# qplot(x,w,color=subject)+facet_wrap(~subject)
m<-pfda(y~x%|%id,driver='single.binary', k=1, df=c(3,3))
bp(6)

#AIC
mtrace(AIC.pfda.single.b)
# mtrace(.single.b.n2L)
# mtrace(.single.b.n2L.1)
# mtrace(.single.b.estimatePCS)
# mtrace(.single.b.updatePCS)
AIC(m)

mtrace.off()
