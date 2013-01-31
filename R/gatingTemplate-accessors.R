setGeneric("gating", function(x, gs, ...) standardGeneric("gating"))
setMethod("gating", signature = c("gatingTemplate", "GatingSet"), definition = function(x, gs, ...) {
			
			
			#gate each node by the topological order
			#maintain the mapping between template node ID and gating set node ID
			#in order to refer back to the template and find the parent 
			gtNodeIDs<-tsort(x)
			for(gtParentID in gtNodeIDs)
			{
				browser()
				curGtNode<-getNodes(gt,gtParentID)
				parentName<-curGtNode["name"]
				parentAlias<-curGtNode["label"]
				curChildren.gt.IDs<-getChildren(gt,gtParentID)
				#check if gates are the same and merge them
				gates<-lapply(curChildren.gt.IDs,function(curChildren.gt.ID){
							getGate(gt,gtParentID,curChildren.gt.ID)
						})
				if(isSame(gates))
					gates<-gates[[1]]
				lapply(gates,function(gate){
							
						})
				gt.child.node<-getNodes(gt,curChildren.gt.ID)
				
				#use curNode as a parent population
				#and try to detect all the branches sourced from curNode
				#and merge them to one gating method if they share the same
				#dims and method names
				#name the resulting gated population by the dest of each branch
				
				#if the current parent doesn't exsit in gs
				#then need to go back to gate the parent first  
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

setMethod("getNodes",signature=c("gatingTemplate"),definition=function(x,y)
{
	
	nodeData(x,y)[[1]]
})

setMethod("getChildren",signature=c("gatingTemplate","character"),definition=function(obj,y)
{
	
	edges(gt,y)[[1]]		
})
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
			
			
			#get gating method name attr from edges
			gm.names<-unlist(lapply(edgeData(x,attr="gatingMethod"),names))
			gm.types<-unique(gm.names)
			#fix the name attr 
			e.colnames<-gsub("\\|","~",names(gm.names))
	
			#encode the method name with color
			nMethods<-length(gm.types)
			gm.col<-RColorBrewer::brewer.pal(nMethods,name="Dark2")
			names(gm.col)<-gm.types
			eCols<-gm.col[gm.names]
			names(eCols)<-e.colnames#restore names 
			
			#encode edge style
			gm.isPoly<-unlist(lapply(gm.names,function(y)ifelse(y=="polyfunctions","poly","regular")))
			edge.styles<-c("solid","dashed")
			names(edge.styles)<-unique(gm.isPoly)
			eStyles<-edge.styles[gm.isPoly]
			names(eStyles)<-e.colnames
			eAttrs<-list(color=eCols
						,lty=eStyles
						)
#			browser()
			
			#encode the node shape and color with isSubsets
			nodeTypes<-unlist(lapply(nodeData(x,attr="isSubsets"),function(y){ifelse(y,"subset","pop")}))
			n.colnames<-names(nodeTypes)
			nodeType<-unique(nodeTypes)
			
			#color
			subset.col<-RColorBrewer::brewer.pal(9,name="Set3")[c(9,7)]
			names(subset.col)<-nodeType
			nFillCol<-subset.col[nodeTypes]
			names(nFillCol)<-n.colnames
			
			#line type(somehow it doesn't work here)
			LineType<-c("solid","dotted")
			names(LineType)<-nodeType
			nLtys<-LineType[nodeTypes]
			names(nLtys)<-n.colnames
			
			nLabels<-unlist(nodeData(x,attr="label"))
			
			nAttrs<-list(label=nLabels
						,fillcolor=nFillCol
						,lty=nLtys)
#			browser()
			plot(x
					,nodeAttrs=nAttrs
					,edgeAttrs=eAttrs
					,attrs=list(graph=list(rankdir="LR",page=c(8.5,11))
											,node=list(fixedsize=FALSE
														,shape="ellipse"
														)
								)
				)
			
#			browser()
			plot.space = par()[["usr"]]
			x1 = plot.space[1]
			x2 = plot.space[2]
			y1 = plot.space[3]
			y2 = plot.space[4]
			
			legend.lty<-ifelse(gm.types=="polyfunctions","poly","regular")
			legend(x1+100,y2-100
					,legend=gm.types
					,title="Gating Methods"
					,col=gm.col
#					,pch="l"
					,lty=edge.styles[legend.lty]
					,cex=0.8
					)
		})
setMethod("names",signature=c("gatingMethod"),definition=function(x)
		{
			
			x@name
		})


