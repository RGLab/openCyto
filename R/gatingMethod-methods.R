# TODO: Add comment
# 
# Author: wjiang2
###############################################################################

setMethod("show",signature=c("gatingMethod"),definition=function(object)
		{
			
			cat("Gating Method: ")
			cat(names(object))
			cat("(")
#			browser()
			chnls<-dims(object)
			cat(paste(paste(names(chnls),chnls,sep="="),collapse=","))
			cat(",")
			cat(parameters(object))
			cat(") \n");
		})


setMethod("names",signature=c("gatingMethod"),definition=function(x)
		{
			
			x@name
		})

setMethod("dims",signature=c("gatingMethod"),definition=function(object)
		{
			
			dims<-strsplit(object@dims,",")[[1]]
			if(length(dims)==1)
				dims<-c(NA,dims)
			else if(length(dims)!=2)
				stop("invalid dimensions!")
			names(dims)<-c("xChannel","yChannel")			
			dims
		})


setMethod("parameters",signature=c("gatingMethod"),definition=function(object)
		{
			
			object@args
		})

setMethod("gating", signature = c("gatingMethod", "GatingSet"), definition = function(x, gs,names, parent,...) {

})