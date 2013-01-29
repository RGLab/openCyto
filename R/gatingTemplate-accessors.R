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
												,parent=parent
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
			cat("\n")
#			browser()
			cat("\twith ",length(object@nodes)," populations\n");
		})
setMethod("plot",signature=c("gatingTemplate"),definition=function(x,y=missing)
		{
			x<-as(x,"graphNEL")
#			browser()
			
			nAttrs<-list(label=unlist(nodeData(x,attr="label")))
			#get gating method name attr from edges
			gm<-unlist(lapply(edgeData(x,attr="gatingMethod"),names))
			eCols<-gm
			#fix the name attr 
			e.colnames<-gsub("\\|","~",names(eCols))
			#get the categorical method names
			gm.txt<-unique(eCols)
			nMethods<-length(gm.txt)
			#encode the method name with color
			gm.col<-RColorBrewer::brewer.pal(nMethods,name="Dark2")
			names(gm.col)<-gm.txt
			eCols<-gm.col[eCols]
			names(eCols)<-e.colnames#restore names 
			eAttrs<-list(color=eCols)
#			browser()
			plot(x
					,nodeAttrs=nAttrs
					,edgeAttrs=eAttrs
					,attrs=list(graph=list(rankdir="LR",page=c(8.5,11))
											,node=list(fixedsize=FALSE
														,fillcolor="gray"
#														,color="white"
														,shape="ellipse"
														)
							
#								,edge=list(labelfontsize="27")
							
								)
			)
			
#			browser()
			plot.space = par()[["usr"]]
			x1 = plot.space[1]
			x2 = plot.space[2]
			y1 = plot.space[3]
			y2 = plot.space[4]
			
			
#			plot(lay)
			legend(x1+100,y2-100
					,legend=gm.txt
					,title="Gating Methods"
					,col=gm.col
#					,pch="l"
					,lty=1
					,cex=0.8
					)
		})
setMethod("names",signature=c("gatingMethod"),definition=function(x)
		{
			
			x@name
		})


