setGeneric("getPrior", function(x,y,...) standardGeneric("getPrior"))
setMethod("getPrior",sig=c("fcObject","character"),definition=function(x,y){
							
			prior<-x@prior
			priorNames<-names(prior)
#			if(is.null(y))
#			{
#
#				if(length(priorNames)==1)
#					prior[[1]]
#				else
#					stop("Need to specify which prior:",paste(priorNames,collapse="or"))
#			}else
#			{
				
				ind<-match(y,priorNames)
				if(is.na(ind))
					stop("prior not found for:",y)
				else
					prior[[ind]]		
#			}
				
		})
setGeneric("getFcFilter", function(x,...) standardGeneric("getFcFilter"))
setMethod("getFcFilter",sig=c("fcObject"),definition=function(x){
			x@fcFilters
			
})


#
setMethod("plot",sig=c("fcObject","ANY"),definition=function(x,y,samples=NULL,posteriors=FALSE,...){
			
			priorNames<-names(x@prior)
			if(is.null(priorNames))
				stop("no valid flowClust prior for this population!")
			if(is.null(y))
			{
				if(length(priorNames)==1)
					y<-priorNames[1]
				else
					stop("Need to specify which channel:",paste(priorNames,collapse="or"))
			}			
#		browser()
		prior<-getPrior(x=x,y=y)
		thisFcFilters<-getFcFilter(x)
		#try to set x,y lim from 
		minX<-min(unlist(lapply(thisFcFilters,function(curFilter)posteriors(x=curFilter,y=y)$min)))
		maxX<-max(unlist(lapply(thisFcFilters,function(curFilter)posteriors(x=curFilter,y=y)$max)))
		x_dens <- seq(minX, maxX, length = 1000)
		
		K<-nrow(prior$Mu0)
		
		
		
		prior_density <-lapply(seq_len(K),function(k)dnorm(x_dens, mean = prior$Mu0[k], sd = sqrt(prior$Omega0[k])))
		minY<-min(unlist(lapply(prior_density,min)))
		maxY<-max(unlist(lapply(prior_density,max)))
		
		plot(x=NULL,type="n"
				,xlim=c(minX,maxX)
				,ylim=c(minY,maxY)
				,xlab=y
				,ylab=""
				,...)
		for(k in seq_len(K)) {
			lines(x_dens, prior_density[[k]], col = rainbow(K)[k], lty = 2, lwd = 1)
#		browser()
			#plot post
			if(posteriors)
			{
				if(is.null(samples))
					samples<-names(thisFcFilters)
				for(samp in samples){
					curFilter<-thisFcFilters[[samp]]
					curPost<-posteriors(x=curFilter,y=y)
					
#					curPrior<-prior[[samp]]
#						browser()
					posterior_density <- dmvt(x_dens, mu = curPost$mu[k,],
							sigma = curPost$sigma[k,,], nu = curPost$nu,
							lambda = curPost$lamdda)$value
					lines(x_dens, posterior_density, col = rainbow(K)[k], lwd = 1)
					#plot gate

					abline(v=curFilter@min[y],col="grey")
				}
				
			}
		}
		
	})	
setMethod("plot",sig=c("fcObject2d","ANY"),definition=function(x,y,samples=NULL,posteriors=FALSE,...){
	message("To be implemented!")			
		})