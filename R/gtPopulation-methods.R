


setMethod("names",signature=c("gtPopulation"),definition=function(x)
		{
			
			x@name
		})

setMethod("alias",signature=c("gtPopulation"),definition=function(object)
		{
			
			object@alias
		})


