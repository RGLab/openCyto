setMethod("gating", signature = c("gatingTemplate","GatingSet"), definition = function(x,y, ...) {
			
#			browser()
			#gate each node by the topological order
			#maintain the mapping between template node ID and gating set node ID
			#in order to refer gating set node ID back to the template ID and find the parent gs node ID
			gt<-x
			gt_node_ids<-tsort(gt)
			node_ids<-cbind(gt=gt_node_ids,gs=NA)
			node_ids[1,"gs"]<-1#fill out default gsid for root node
			for(i in 1:nrow(node_ids))
			{
				
				#get parent node to gate
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
				
				gates<-lapply(gt_children_ids,function(i){
							getGate(gt,gt_parent_id,i)
						})
#				browser()
#					isSame<-unique(gates)
#					if(all(!isSame))
#					{
				#do the gating for each unique gate
				for(gate in unique(gates))
				{
					

					#select the children that are associated with this gate
					children_ind<-which(unlist(lapply(gates,"identical",gate)))
					cur_gt_children_ids<-gt_children_ids[children_ind]
					#get population info
					pops<-lapply(cur_gt_children_ids
							,function(gt_children_id){
								getNodes(gt,gt_children_id)
							})
#							browser()					
					#pass the pops and gate to gating routine
					gs_node_ids<-gating(x=gate
										,y
										,parent=as.integer(gs_parent_id)
										,gtPops=pops
										)		
					#upodate gs node ids
					ind<-match(gt_children_ids,node_ids[,"gt"])
#					browser()
					node_ids[ind,"gs"]<-gs_node_ids
				}	
#					}else
#					{
#						
#						##quadgate
#						##TODO:currently each parent should either have one quadgate or multiple separate gates
#						#we may need to figure out how to extend it to a more generic scenario
#						if(length(which(isSame))<=3)#quadgate should only have 3 duplicates
#						{
#							gate<-gates[[1]]
#							pops<-lapply(gt_children_ids
#									,function(gt_children_id){
#										getNodes(gt,gt_children_id)
#									})
#							pops
#						}else 
#							stop("don't know how to handle quadgate with more than 4 sub-populations!")
#					}
				
				
				
				
#				}else if(all(!is.na(gs_node_ids)))
#				{
#					
#					message("'",alias(getNodes(gt,gt_children_ids)),"' already exists,skip gating!")
#				}else
#					stop("Don't know how to handle partially gated children node yet!")
			}
			
			
			
			
			
			
			message("finished.")
		})

		
setMethod("gating", signature = c("gtMethod", "GatingSet")
		, definition = function(x, y,gtPops, parent
				,num_nodes = 1, parallel_type = c("multicore", "sock")
				,plot = FALSE, xbin = 128, ...) 
		{
			
			require('parallel')
			
			
			args<-parameters(x)
			
#			browser()
			gm<-paste(".",names(x),sep="")
			
			dims<-dims(x)
			xChannel<-unname(dims["xChannel"])
			yChannel<-unname(dims["yChannel"])
			if(!is.na(xChannel))
				xChannel_s<-paste("'",xChannel,"'",sep="")
			else
				xChannel_s<-xChannel
			if(!is.na(yChannel))
				yChannel_s<-paste("'",yChannel,"'",sep="")
			
			popAlias<-unlist(lapply(gtPops,alias))
			popNames<-unlist(lapply(gtPops,names))
			popIds<-unlist(lapply(gtPops,"slot","id"))
			gs_nodes<-getChildren(y[[1]],getNodes(y[[1]],parent))
			

#			browser()
			
			if (!any(grepl(popAlias, gs_nodes))) 
			{
				message("Population '",paste(popAlias,collapse=","),"'")
				
				parent_data <- getData(y, parent)
				parallel_type <- match.arg(parallel_type)
				
				# Splits the flow set into a list, where each element in the list is a
				# flowSet containg one flow frame.
				# Below, we access manually the individual flow frame with current_frame[[1]].
				fslist <- split(parent_data, sampleNames(parent_data))
				
				#construct method call
				thisCall<-substitute(f1())
				thisCall[["X"]]<-quote(fslist)#set data
				thisCall[["FUN"]]<-as.symbol(gm)#set gating method
				thisCall[["xChannel"]]<-xChannel#set x,y channel
				thisCall[["yChannel"]]<-yChannel
				
#						browser()
				#append positive argument based on population name
				if(length(popNames)==1) #1d gate
				{
					#parse the +.- from pop name
					if(grepl("-$",popNames))
						positive=FALSE
					else if(grepl("+$",popNames))
						positive=TRUE
					else
						stop("invalid population name!Name should end with '+' or '-' symbol.")
					
					thisCall[["positive"]]<-positive
					
				}else if(length(popNames)<=4)#QuadGate
				{
					message("quadGate...")
#				grepl(popNames				
				}else
					stop("don't know how to handle quadgate with more than 4 sub-populations!")
#			browser()
				#parse the string into named vector
				paired_args<-strsplit(args,split="\\,")[[1]]
				named_args<-unlist(lapply(paired_args
										,function(paired_arg){
											curPair<-strsplit(paired_arg,split="=")[[1]]
											curPair_value<-curPair[2]
											#strip while spaces and quote symbols
											curPair_value<-sub("[ \t\n\r]*$", "", sub("^[ \t\n\r]*", "", curPair_value))
											curPair_value<-gsub("\"","",curPair_value)
#											browser()
											names(curPair_value)<-sub("[ \t\n\r]*$", "", sub("^[ \t\n\r]*", "", curPair[1]))
											curPair_value
											}
										)
									)
				named_args<-as.list(named_args)
				#try to convert to numeric if applicable
				named_args<-lapply(named_args,function(cur_arg){
							cur_arg_new<-as.numeric(cur_arg)
							if(!is.na(cur_arg_new))
								cur_arg<-cur_arg_new
							cur_arg
						})
#				browser()
				#parse args for flowClust 
				#prior estimation is done separately from flowClust routine 
				#because piror_flowClust1d requires the entire parent flowSet yet flowClust only takes one flowFrame
				#####################################################################
				if(grepl("^\\.flowClust\\.[12]d$",gm))
				{
#				browser()
					if(gm==".flowClust.1d")
					{
						#get the value of neg and pos
						
						if(all(is.element(c("pos","neg"),names(named_args))))
						{
							neg_cluster<-as.integer(named_args["neg"])			
							K<-as.integer(named_args["pos"])+neg_cluster
						}else
						{
							message("either 'neg' or 'pos' argument is missing!Using default setting:neg=1,pos=1")
							neg_cluster<-as.integer(1)			
							K<-2
						}
						
						# Elicitation of priors for flowClust
						prior <- list()
						if(!is.na(xChannel))
							prior$xChannel <- prior_flowClust1d(flow_set = parent_data,channel = xChannel, K = K)
						
						
						prior$yChannel <- prior_flowClust1d(flow_set = parent_data,channel = yChannel, K = K)
						
#						browser()
						#replace neg and pos and convert the named vector back to string
						named_args[["K"]]<-K
						named_args[["neg_cluster"]]<-neg_cluster
						named_args[["prior"]]<-prior
						neg_ind<-match("neg",names(named_args))
						if(!is.na(neg_ind))
							named_args<-named_args[-neg_ind]	
						pos_ind<-match("pos",names(named_args))
						if(!is.na(pos_ind))
							named_args<-named_args[-pos_ind]
					}else
					{
						#get the value of neg and pos
						K<-as.integer(named_args["K"])			
						named_args[["K"]]<-K
						
						
						
						# Elicitation of priors for flowClust
						prior_list <- prior_flowClust2d(fr = parent_data[[1]]
								,xChannel = xChannel
								,yChannel = yChannel
								, K = K)
						
						named_args[["prior_list"]]=prior_list
						thisCall[["positive"]]<-NULL #remove positive arg since 2D gate doesn't understand it
					}
					
				}
				#add args to thisCall
				for(arg in names(named_args))
					thisCall[[arg]]<-named_args[[arg]]
				
				
#			browser()
				##choose serial or parallel mode
				if (num_nodes > 1) 
				{
					message("Running in parallel mode with ", num_nodes, " nodes.")
					if (parallel_type == "multicore") 
					{
						thisCall[[1]]<-quote(mclapply)
						thisCall[["mc.cores"]]<-num_nodes
						flist<-eval(thisCall)
						
					}else 
					{
						cl <- makeCluster(num_nodes, type = "SOCK")
						thisCall[[1]]<-quote(parLapply)
						thisCall[["cl"]]<-cl
						#replace FUN with fun for parLapply
						thisCall["fun"]<-thisCall["FUN"]
						thisCall["FUN"]<-NULL
						flist<-eval(thisCall)
						stopCluster(cl)
					}
				}else
				{
					thisCall[[1]]<-quote(lapply)#select loop mode
					flist<-eval(thisCall)
				}
				
				#we expect a filter/filters from gm
				
				if(class(flist[[1]])=="filters")#reconstruct filterlist out of filterslist
				{
					#right now we consider quadgate as the only use case for filters
					gs_node_id<-NULL
#					browser()
					#clock-wise starting from top-left quadrant											
					quadPatterns<-c(".+-.+\\+$"   #top left  -+
									,".+\\+.+\\+$" #top right ++
									,".+\\+.+-$"  #bottom right +-
									,".+-.+-$")   #bottom left	--									
					for(i in 1:length(popNames))
					{

						curAlias<-popAlias[i]
						curPop<-popNames[i]
						curPopId<-popIds[i]

						quadInd<-which(unlist(lapply(quadPatterns,grepl,curPop)))
						#fetch appropriate filter based on the quadrant ind
						curFlist <- lapply(flist, function(curFilters) {
#								curFilters[[quadInd]]@filterId <- curAlias
									curFilters[[quadInd]]
								})
						curFlist <- filterList(curFlist)
						cur_gs_node_id <- add(gs, curFlist, parent = parent,name=curAlias)	
						recompute(y, cur_gs_node_id)
						gs_node_id<-c(gs_node_id,cur_gs_node_id)
					}
					
					
					
				}else
				{
					gs_node_id <- add(y, flist, parent = parent,name=popAlias)
					recompute(gs, gs_node_id)
					
				}
				
				message("done.")
				
			}else
			{
				message("Skip gating!Population '",paste(popAlias,collapse=","),"' already exists.")
#			browser()
				gs_node_id<-getChildren(y[[1]],parent)
				#select the corresponding gs node id by matching the node names
				gs_node_name<-getNodes(y[[1]])[gs_node_id]
				gs_node_id<-gs_node_id[match(popAlias,gs_node_name)]
			}
			
			if (plot) {
				print(plotGate(y, gs_node_id, xbin = xbin, pos = c(0.5, 0.8)))
			}
			names(gs_node_id)<-popIds
			gs_node_id
		})
		
setMethod("gating", signature = c("polyFunctions", "GatingSet")
				, definition = function(x, y,gtPops, parent
						,num_nodes = 1, parallel_type = c("multicore", "sock")
						,plot = FALSE, xbin = 128, ...) 
{
	
	require('parallel')
	
	args<-parameters(x)
	gm<-paste(".",names(x),sep="")
	popAlias<-unlist(lapply(gtPops,alias))
	popNames<-unlist(lapply(gtPops,names))
	gs_nodes<-getChildren(y[[1]],getNodes(y[[1]],parent))
	
	
	message("Population '",paste(popAlias,collapse=","),"'")
	
	
	refNodes<-strsplit(args,split=":")[[1]]
	
	nMarkers <- length(refNodes)
	## all the comibnations of A & B & C
	# The 'permutations' function is from the 'gregmisc' package on CRAN.
	require('gregmisc')
	opList <- permutations(n = 1, r = nMarkers - 1, c("&"), repeats = TRUE)
	isNotList <- permutations(n = 2, r = nMarkers, c("!", ""), repeats = TRUE)
	polyExprsList <- apply(opList, 1, function(curOps) {
				apply(isNotList, 1, function(curIsNot) {
							polyExprs <- curIsNot
							polyExprs[-1] <- paste0(curOps, curIsNot[-1])
							
							paste(paste0(polyExprs, refNodes), collapse = "")
						})
			})
	polyExprsList <- as.vector(polyExprsList)
	
	#actual gating
	lapply(polyExprsList, function(polyExpr) {
				
				tNodes <-polyExpr
				polyExpr <- as.symbol(polyExpr)
#					browser()
				if (!any(tNodes %in% gs_nodes)) {
					message(tNodes, " gating...")
					bf <- eval(substitute(booleanFilter(x), list(x = polyExpr)))
					bf@filterId <- tNodes
					invisible(gs_node_id <- add(gs, bf, parent = parent))
					invisible(recompute(gs, gs_node_id))
				}else
					message("Skip gating!Population '",tNodes,"' already exists.")
			})					
	
	
	
	message("done.")
		
	
	
	return (-1)
})
		
## wrappers for the different gating routines
.singletGate <- function(fs, xChannel = "FSC-A",yChannel = "FSC-H", prediction_level = 0.99,...) 
{
	require('flowStats')
	# Creates a list of polygon gates based on the prediction bands at the minimum and maximum
	# x_channel observation using a robust linear model trained by flowStats.
	
	
	singletGate(fs[[1]], area = xChannel
			, height = yChannel
			,prediction_level = prediction_level
	)$gate
	
}
.flowClust.1d<-function(fs, xChannel = NA,yChannel,tol=1e-3,prior=NULL,filterId=""
						,usePrior="yes",split=TRUE,...)
{
	
	require("flowClust")
	
	
	fr<-fs[[1]]
#		browser()			
	if(is.na(xChannel))#1d gate			
		flowClust.1d(fr = fr
				,params = yChannel
				,tol = tol
				,filterId = filterId
				,prior = prior$yChannel
				,usePrior=usePrior
				
				,...)
	else#quadgate 
	{
		xParam <- getChannelMarker(fr, xChannel)
		yParam <- getChannelMarker(fr, yChannel)
		xChannel<-xParam$name
		xStain<-xParam$desc
		yChannel<-yParam$name
		yStain<-yParam$desc
#			browser()
		filter1 <- flowClust.1d(fr = fr
				, params = xChannel, tol = tol
				,filterId = as.character(xStain)
				,prior = prior$xChannel
				,usePrior=usePrior
				, ...)
		################################      
		#flowClust on xParam
		################################
		
		cut.x <- filter1@min
		if (split) {
			################################################################
			## split the data for the further gating
			################################################################
			newFrList <- flowCore::split(fr, filter1)
			
			negFr <- newFrList[[paste0(xStain, "-")]]
			posFr <- newFrList[[paste0(xStain, "+")]]
			
			filter2 <- flowClust.1d(negFr, params = yChannel,
					usePrior = usePrior
					,filterId = as.character(yStain)
					,prior = prior$yChannel
					, ...)
			cut.y.l <- filter2@min
			
			filter3 <- flowClust.1d(posFr, params = yChannel,
					usePrior = usePrior
					,filterId = as.character(yStain)
					,prior = prior$yChannel
					, ...)
			
			cut.y.r <- filter3@min
		} else {
			filter2 <- flowClust.1d(fr, params = yChannel,
					usePrior = usePrior
					,filterId = as.character(yStain)
					,prior = prior$yChannel
					, ...)
			cut.y.l <- cut.y.r <- filter2@min
		}
		
		###############################################################     
		#construct rectangleGates based on the cuts and popNames,clock-wise
		###############################################################
		gateList <- new("filters")
		
		chnls<-c(xChannel, yChannel)
		markers<-c(xStain,yStain)
		
		coord <- list(c(-Inf, cut.x), c(cut.y.l, Inf))
		names(coord) <- as.character(chnls)
		gateList[[paste(paste0(markers, c("-", "+")), collapse="")]] <- rectangleGate(coord)
		
		coord <- list(c(cut.x, Inf), c(cut.y.r, Inf))
		names(coord) <- as.character(chnls)
		gateList[[paste(paste0(markers, c("+", "+")), collapse="")]] <- rectangleGate(coord)
		
		coord <- list(c(cut.x, Inf), c(-Inf, cut.y.r))
		names(coord) <- as.character(chnls)
		gateList[[paste(paste0(markers, c("+", "-")), collapse="")]] <- rectangleGate(coord)
		
		coord <- list(c(-Inf, cut.x), c(-Inf, cut.y.l))
		names(coord) <- as.character(chnls)
		gateList[[paste(paste0(markers, c("-", "-")), collapse="")]] <- rectangleGate(coord)
		
		gateList
	}
	
}
.flowClust.2d<-function(fs, xChannel,yChannel,usePrior="yes",...)
{
	
	require("flowClust")
	
	
	fr<-fs[[1]]
#	browser()
	flowClust.2d(fr = fr, xChannel = xChannel, yChannel = yChannel,usePrior=usePrior, ...)
	
}
.rangeGate<-function(fs, xChannel = NA,yChannel,absolute = FALSE,filterId="",...)
{
	require('flowStats')
	rangeGate(fs[[1]], stain = yChannel, inBetween = TRUE, absolute = absolute,
			filterId = filterId, ...)	
}


.quantileGate<-function(fs, xChannel = NA,yChannel,probs = 0.999,filterId="",...)
{
	quantileGate(fr = fs[[1]], probs = probs, stain = yChannel, filterId = filterId, ...)	
}

.quadrantGate<-function(fs, xChannel = NA,yChannel,...)
{
	require('flowStats')
	quadrantGate(fs[[1]], stain = c(xChannel,yChannel), absolute = FALSE, inBetween = TRUE, ...)			
}

