.plotMeanWithPC1<-function(tp,tt,tf,Dg,hcol=0,name.t,name.y,i=1){
	sup <- if(i==1) "st" else if(i==2) "nd" else if(i==3) "rd" else "th"
	{ qplot(tp,tt,geom='line',main=substitute(paste("Mean", phantom(.)%+-% i^sup," component for ",name.t),list(i=i,sup=sup,name.t=name.t)))+
	labs(y=name.y,x=name.t)+
	geom_ribbon(aes(ymin=tt-2*tf*sqrt(Dg),ymax=tt+2*tf*sqrt(Dg)),fill=hcl(hcol,l=90))+
	geom_ribbon(aes(ymin=tt-1*tf*sqrt(Dg),ymax=tt+1*tf*sqrt(Dg)),fill=hcl(hcol,l=70))+
	geom_line(aes(tp,tt))}
}
plotMeanWithPC.pfda.additive<-function(m,ask=TRUE,show=TRUE){
	with(m,{
		if (ask) {
			if (.Device != "null device") {
				oldask <- grDevices::devAskNewPage(ask = TRUE)
			if (!oldask) 
				on.exit(grDevices::devAskNewPage(oldask), add = TRUE)
			}
			op <- options(device.ask.default = TRUE)
			on.exit(options(op), add = TRUE)
		}
		p<-if(show) print else invisible
		total.npc<-length(Dg)+length(Dd)
		col.inc<-trunc(360/total.npc)
		graphs<-vector('list',total.npc)
		tp<-seq(min(tbase@knots),max(tbase@knots),length=100)
		Bt<-evaluate(tbase,tp)
		for(i in seq_along(Dg)){
			p(graphs[[i]]<-.plotMeanWithPC1(tp,Bt%*%tt,Bt%*%tf[,i],Dg[i],hcol=(i-1)*col.inc,attr(m,'name.t'),attr(m,'name.y'),i))
		}
		xp<-seq(min(xbase@knots),max(xbase@knots),length=100)
		Bx<-evaluate(xbase,xp)[,-1]
		for(i in seq_along(Dd)){
			p(graphs[[i+length(Dg)]]<-.plotMeanWithPC1(xp,Bx%*%tx,Bx%*%tg[,i],Dd[i],hcol=(i-1+length(Dg))*col.inc,attr(m,'name.x'),attr(m,'name.y'),i))
		}
		return(graphs)
	})
}
plotMeanWithPC<-function(m,...)UseMethod("plotMeanWithPC")

mtrace(plotMeanWithPC.pfda.additive)
