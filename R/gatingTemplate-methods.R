
setMethod("getNodes",signature=c("gatingTemplate"),definition=function(x,y)
{
#	browser()
	if(missing(y))
		nodeData(x,attr="pop")
	else
		nodeData(x,y,"pop")[[1]]
})

setMethod("getChildren",signature=c("gatingTemplate","character"),definition=function(obj,y)
{
#	browser()
	edges(obj,y)[[1]]	
})
setMethod("getGate",signature=c("gatingTemplate","character"),definition=function(obj,y,z)
{
		edgeData(obj,from=y,to=z,attr="gtMethod")[[1]]
	
})




setMethod("show",signature=c("gatingTemplate"),definition=function(object)
		{

			cat("--- Gating Template: ")
			cat(object@name)
			cat("\n")
#			browser()
			cat("\twith ",length(object@nodes)," populations defined\n");
		})
setMethod("plot",signature=c("gatingTemplate"),definition=function(x,y=missing)
		{
			
			#get gating method name attr from edges
			gm.names<-unlist(lapply(edgeData(x,attr="gtMethod"),names))
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
			gm.isPoly<-unlist(lapply(gm.names,function(y)ifelse(y=="polyFunctions","poly","regular")))
			edge.styles<-c("solid","dashed")
			names(edge.styles)<-unique(gm.isPoly)
			eStyles<-edge.styles[gm.isPoly]
			names(eStyles)<-e.colnames
			eAttrs<-list(color=eCols
						,lty=eStyles
						)
#			browser()
			
			#encode the node shape and color with isSubsets
			nodeTypes<-unlist(lapply(nodeData(x,attr="pop"),function(y){ifelse(class(y)=="gtSubsets","subset","pop")}))
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
			
#			browser()
			nLabels<-unlist(lapply(nodeData(x,attr="pop"),alias))
			
			nAttrs<-list(label=nLabels
						,fillcolor=nFillCol
						,lty=nLtys)
			
			plot(as(x,"graphNEL")
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
			
			legend.lty<-ifelse(gm.types=="polyFunctions","poly","regular")
			legend(x1+100,y2-100
					,legend=gm.types
					,title="Gating Methods"
					,col=gm.col
#					,pch="l"
					,lty=edge.styles[legend.lty]
					,cex=0.8
					)
		})

#extend flowWorkspace to store gatingSet 
#setMethod("gatingTemplate",signature(x="GatingSet"),function(x){
#			
#			gatingTemplate(x[[1]])			
#			
#		})
#setMethod("gatingTemplate",signature(x="GatingHierarchy"),function(x){
#			
#				
#			r<-nodeDataDefaults(x@tree,"data")[["data"]];
#			r$gt
#
#		})
#setReplaceMethod("gatingTemplate",signature(x="GatingSet"),function(x,value){
#	
#			
#	r<-nodeDataDefaults(x[[1]]@tree,"data")[["data"]];
#	r$gt<-value
#	x
#	
#})

