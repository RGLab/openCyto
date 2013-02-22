#
#
###############################################################################

setMethod("plot",sig=c("fcObject","missing"),definition=function(x,y="missing",samples=NULL,posteriors=FALSE,...){
		
		prior<-x@prior$yChannel
		
		post<-x@posteriors
		
		minX<-min(unlist(lapply(post,function(curPost)curPost$min)))
		maxX<-max(unlist(lapply(post,function(curPost)curPost$max)))
		x_dens <- seq(minX, maxX, length = 1000)
		
		K<-nrow(prior$Mu0)
		
		
		
		prior_density <-lapply(seq_len(K),function(k)dnorm(x_dens, mean = prior$Mu0[k], sd = sqrt(prior$Omega0[k])))
		minY<-min(unlist(lapply(prior_density,min)))
		maxY<-max(unlist(lapply(prior_density,max)))
		
		plot(x=NULL,type="n"
				,xlim=c(minX,maxX)
				,ylim=c(minY,maxY)
				,xlab=prior$channelName
				,ylab="",...)
		for(k in seq_len(K)) {
			lines(x_dens, prior_density[[k]], col = rainbow(K)[k], lty = 2, lwd = 1)
#		browser()
			#plot post
			if(posteriors)
			{
				if(is.null(samples))
					samples<-names(post)
				for(samp in samples){
					curPost<-post[[samp]]
					posterior_density <- dmvt(x_dens, mu = curPost$mu[k,],
							sigma = curPost$sigma[k,,], nu = curPost$nu,
							lambda = curPost$lamdda)$value
					lines(x_dens, posterior_density, col = rainbow(K)[k], lwd = 1)
				}
				
			}
		}
		
	})	
			
