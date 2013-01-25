setMethod("gating", signature = c("gatingTemplate", "GatingSet"), definition = function(x, gs, ...) {
			
			
			#start from root
			parent <-"root"
			
			while(length(parent)>0){
				#get entries by parent
			browser()
			children<-x[x[,"parent"]==parent,]
			#there may be multiple entries for each gate				
			gates<-unique(children[,c("dims","method","args")])
			apply(gates,1,function(gate){
						#deal with each gate applied to the current pop
						browser()
						gm<-unname(gate["method"])
						pop<-unname(gate["pop"])
						dims<-unname(gate["dims"])
						args<-as.symbol(unname(gate["args"]))
						#pass the arguments to the gating function
						thisCall<-substitute(
										get(gm)(x=x
												,wf=gs
												,pViewName=parent
												,name=pop
												,xChannel = dims[1]
												, yChannel = dims[2]
												,args
												)
										,list(args=args))
						eval(thisCall)
					})
			}
			
			
			
			message("finished.")
		})

#setMethod("getChildren",signature=c("gatingTemplate","character"),definition=function(obj,y)
#{
#	
#	as.character(getGate(obj,y)[,"pop"])		
#})
#setMethod("getGate",signature=c("gatingTemplate","character"),definition=function(obj,y)
#{
#			subset(obj,parent==y)		
#		})
setMethod("show",signature=c("gatingTemplate"),definition=function(object)
		{

			cat("--- Gating Template: ")
			cat(object@name)
#			cat("\ntransformations:")
#			print(getTransformations(object))
#			cat("\ncompensations:")
#			print(getCompensationMatrices(object))
#			cat(" ---")
			cat("\n")
#			cat("\nquadrantGate",fill=T,labels="TH,NL:")
#			cat(tsubMarkers(object),fill=T,labels="Markers:",sep=",")
#			cat("\n")
		})
#setMethod("getTransformations",signature=c("gatingTemplate"),definition=function(x)
#		{
##			browser()
#			x@transformation
#		})
#setMethod("getCompensationMatrices",signature=c("gatingTemplate"),definition=function(x)
#		{
#			
#			x@compensation
#		})


