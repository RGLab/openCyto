setMethod("show",sig=c("fcFilter"),definition=function(object){
			cat("fcFilter:")
			callNextMethod(object)
			
		})

setGeneric("posteriors", function(x,y,...) standardGeneric("posteriors"))
setMethod("posteriors",sig=c("fcFilter","ANY"),definition=function(x,y="missing"){
			x@posteriors
			
		})
setMethod("posteriors",sig=c("fcFilter","character"),definition=function(x,y){
			post<-posteriors(x)
			
					
			postNames<-names(post)
#			if(is.null(y))
#			{
#				if(length(postNames)==1)
#					post[[1]]
#				else
#					stop("Need to specify which prior:",paste(postNames,collapse="or"))
#			}else
#			{
				ind<-match(y,postNames)
				if(is.na(ind))
					stop("FcFilter not found for:",y)
				else
					post[[ind]]		
#			}
			
		})