
#
setMethod("plot",sig=c("fcFilterList","ANY"),definition=function(x,y,samples=NULL,posteriors=FALSE,...){
#			browser()
			prior1<-priors(x[[1]],y)
            
#			if(is.null(priorNames))
#				stop("no valid flowClust prior for this population!")
			if(is.null(y))
			{
                priorNames<-names(prior1)
				if(length(prior1)==1)
					y<-priorNames
				else
					stop("Need to specify which channel:",paste(priorNames,collapse="or"))
			}			
        prior1<-priors(x[[1]],y)#refetch prior according to modified y	
		#try to set x,y lim from 
		minX<-min(unlist(lapply(x,function(curFilter)posteriors(x=curFilter,y=y)$min)))
		maxX<-max(unlist(lapply(x,function(curFilter)posteriors(x=curFilter,y=y)$max)))
		x_dens <- seq(minX, maxX, length = 1000)
#		browser()
        #TODO:right now the priors from the first sample is used to determine y scale
		K<-nrow(prior1$Mu0)
		prior_density <-lapply(seq_len(K),function(k)dnorm(x_dens
                                                            , mean = prior1$Mu0[k]
                                                            , sd = sqrt(prior1$Omega0[k])))
		minY<-min(unlist(lapply(prior_density,min)))
		maxY<-max(unlist(lapply(prior_density,max)))
		
		plot(x=NULL,type="n"
				,xlim=c(minX,maxX)
				,ylim=c(minY,maxY)
				,xlab=y
				,ylab=""
				,...)
		
#		browser()
			#plot post
			if(posteriors)
			{
				if(is.null(samples))
					samples<-names(x)
				for(samp in samples){
					curFilter<-x[[samp]]
                    
					curPost<-posteriors(x=curFilter,y=y)
					
					curPrior<-priors(x=curFilter,y=y)
                    K<-nrow(curPrior$Mu0)
                    prior_density <-lapply(seq_len(K),function(k)dnorm(x_dens
                                                                         , mean = curPrior$Mu0[k]
                                                                         , sd = sqrt(curPrior$Omega0[k])
                                                                        )
                                           )

                    for(k in seq_len(K)) {
                      lines(x_dens, prior_density[[k]], col = rainbow(K)[k], lty = 2, lwd = 1)
                      
					posterior_density <- dmvt(x_dens, mu = curPost$mu[k,],
							sigma = curPost$sigma[k,,], nu = curPost$nu,
							lambda = curPost$lamdda)$value
					lines(x_dens, posterior_density, col = rainbow(K)[k], lwd = 1)
					#plot gate

		
					g<-c(curFilter@min[y],curFilter@max[y])
					ind<-which(!is.infinite(g))
					abline(v=g[ind],col="grey")
				  }
				
			  }
		  }
		
	})	
setMethod("plot",sig=c("filterList","ANY"),definition=function(x,y,samples=NULL,posteriors=FALSE,...){
      message("Not valid flowClust filter results!")      
    })
