setGeneric("gating", function(x, gs, ...) standardGeneric("gating"))
setMethod("gating", signature = c("gatingTemplate", "GatingSet"), definition = function(x, gs, ...) {
			
			
			#gate each node by the topological order
			#maintain the mapping between template node ID and gating set node ID
			#in order to refer gating set node ID back to the template ID and find the parent gs node ID 
			gt_node_ids<-tsort(x)
			node_ids<-cbind(gt=gt_node_ids,gs=NA)
			node_ids[1,"gs"]<-1#fill out default gsid for root node
			for(i in 1:nrow(node_ids))
			{
				browser()
				gt_parent_id<-node_ids[i,"gt"]
				gs_parent_id<-node_ids[i,"gs"]
				gt_parent_pop<-getNodes(gt,gt_parent_id)
				parent_name<-names(gt_parent_pop)
				parent_alias<-alias(gt_parent_pop)
				#detect all the branches/gates sourced from gt_parent_id
				gt_children_ids<-getChildren(gt,gt_parent_id)
				#locate the gt node id in the map and check if gs node id already exsits				
#				ind<-match(gt_children_ids,node_ids[,"gt"])
#				gs_node_ids<-node_ids[ind,"gs"]
#				if(all(is.na(gs_node_ids)))
#				{	
					#chec
					gates<-lapply(gt_children_ids,function(i){
								getGate(gt,gt_parent_id,i)
							})
					isSame<-duplicated(gates)
					if(all(!isSame))
					{
						#do the gating for each unique gate
						for(gate in gates)
						{
									
							
							pops<-lapply(gt_children_ids
											,function(gt_children_id){
												getNodes(gt,gt_children_id)
											})

							browser()
							gs_node_ids<-gating(x=gate
												,gs=gs
												,parent=as.integer(gs_parent_id)
												,gtPops=pops
											)
							#upodate gs node ids
							ind<-match(gt_children_ids,node_ids[,"gt"])
							node_ids[ind,"gs"]<-gs_node_ids
						}	
					}else
					{
						##quadgate
						##TODO:currently each parent should either have one quadgate or multiple separate gates
						#we may need to figure out how to extend it to a more generic scenario
						if(length(which(isSame))<=3)#quadgate should only have 3 duplicates
						{
							
						}else stop("don't know how to handle quadgate with more than 4 sub-populations!")
					}
					
					
					
				
#				}else if(all(!is.na(gs_node_ids)))
#				{
#					
#					message("'",alias(getNodes(gt,gt_children_ids)),"' already exists,skip gating!")
#				}else
#					stop("Don't know how to handle partially gated children node yet!")
			}
			
			
			
			
			
			
			message("finished.")
		})

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
	
	edges(gt,y)[[1]]	
})
setMethod("getGate",signature=c("gatingTemplate","character"),definition=function(obj,y,z)
{
		edgeData(obj,from=y,to=z,attr="gatingMethod")[[1]]
	
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


