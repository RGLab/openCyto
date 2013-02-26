
setGeneric("gating", function(x,y,...) standardGeneric("gating"))
setMethod("gating", signature = c("gatingTemplate","GatingSetInternal"), definition = function(x,y,env_fct=NULL, ...) {
			gt<-x
#			browser()
			#gate each node by the topological order
			#maintain the mapping between template node ID and gating set node ID
			#in order to refer gating set node ID back to the template ID and find the parent gs node ID
			
			if(!is.null(env_fct))
			{
				#use the fcTree if already exists 
				if(exists("fct",env_fct))
				{
					fct<-get("fct",env_fct)
				}else
				{
					#create one from gt if not
					fct<-fcTree(gt)
					assign("fct",fct,env_fct)		
				}
				
			}
			
			
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
#				browser()
				gt_children_ids<-getChildren(gt,gt_parent_id)
				
				
				gates<-lapply(gt_children_ids,function(i){
							getGate(gt,gt_parent_id,i)
						})
				#do the gating for each unique gate
				for(gate in unique(gates))
				{
					
#					browser()
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
					res<-gating(x=gate
										,y
										,parent=as.integer(gs_parent_id)
										,gtPops=pops
										,...
										)
					gs_node_ids<-res[["gs_node_id"]]										
					fcObjList<-res[["fcObjList"]]													
					for(curGtId in names(gs_node_ids))
					{
#						browser()
						#upodate gs node ids
						curGsId<-gs_node_ids[[curGtId]]
						ind<-match(curGtId,node_ids[,"gt"])
						node_ids[ind,"gs"]<-curGsId
						#update fct
						if(!is.null(env_fct))
						{
#						browser()
							if(!is.null(fcObjList)&&length(fcObjList)>0)
							{
								nodeData(env_fct$fct,curGtId,"fcObj")<-fcObjList[[curGtId]]
							}
							
							
						}	
					}
					
					
				}	
			}
			
			message("finished.")
			
		})

		
setMethod("gating", signature = c("gtMethod", "GatingSet")
		, definition = function(x, y,gtPops, parent
				,num_nodes = 1, parallel_type = c("multicore", "sock")
				,plot = FALSE, xbin = 128,...) 
		{
			
			require('parallel')
			
			
			args<-parameters(x)
			
#			browser()
			gm<-paste(".",names(x),sep="")
			
			dims<-dims(x)
			xChannel<-unname(dims["xChannel"])
			yChannel<-unname(dims["yChannel"])
			
			popAlias<-unlist(lapply(gtPops,alias))
			popNames<-unlist(lapply(gtPops,names))
			popIds<-unlist(lapply(gtPops,"slot","id"))
			gs_nodes<-getChildren(y[[1]],getNodes(y[[1]],parent))
			

#			browser()
			fcObjList<-list()
			
			
			if (!any(sapply(popAlias,function(a)any(grepl(a, gs_nodes)))))
			{
				message("Population '",paste(popAlias,collapse=","),"'")
				
				parent_data <- getData(y, parent)
				parallel_type <- match.arg(parallel_type)
				##get the accurate channel name by matching to the fr
				if(!is.na(xChannel))
				{
					xParam <- getChannelMarker(parent_data[[1]], xChannel)
					xChannel<-as.character(xParam$name)
				}
				yParam <- getChannelMarker(parent_data[[1]], yChannel)
				yChannel<-as.character(yParam$name)
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

				paired_args<-as.list(as.list(parse(text=paste("c(",args,")")))[[1]])[-1]
				names(paired_args)<-tolower(names(paired_args))

				prior_list<- list()
				#prior estimation is done separately from flowClust routine 
				#because piror_flowClust1d requires the entire parent flowSet yet flowClust only takes one flowFrame
				#####################################################################
				if(grepl("^\\.flowClust\\.[12]d$",gm))
				{
#				browser()
					
					if(gm==".flowClust.1d")
					{
						#get the value of neg and pos
						
						if(all(is.element(c("pos","neg"),names(paired_args))))
						{
							neg_cluster<-as.integer(paired_args["neg"])			
							K<-as.integer(paired_args["pos"])+neg_cluster
						}else
						{
							message("Using default setting:neg=1,pos=1")
							neg_cluster<-as.integer(1)			
							K<-2
						}
#						browser()
						# Elicitation of priors for flowClust
						if(!is.na(xChannel))
							prior_list[[xChannel]] <- prior_flowClust1d(flow_set = parent_data,channel = xChannel, K = K)
						
						prior_list[[yChannel]] <- prior_flowClust1d(flow_set = parent_data,channel = yChannel, K = K)
						
#						browser()
						#replace neg and pos and convert the named vector back to string
						paired_args[["K"]]<-K
						paired_args[["neg_cluster"]]<-neg_cluster
						paired_args[["prior"]]<-prior_list
						neg_ind<-match("neg",names(paired_args))
						if(!is.na(neg_ind))
							paired_args<-paired_args[-neg_ind]
						
						pos_ind<-match("pos",names(paired_args))
						if(!is.na(pos_ind))
							paired_args<-paired_args[-pos_ind]
						
						
						
					}else
					{
						#get the value of neg and pos
						if(is.element(c("k"),names(paired_args)))
						{
							K<-as.integer(paired_args["k"])			
						}else
						{
							message("'K' argument is missing!Using default setting:K=2")
							K<-2
						}	
						paired_args[["k"]]<-K						
						names(paired_args)[match("k",names(paired_args))]<-"K"#restore K to capital letter
						
						# Elicitation of priors for flowClust
						prior_list <- prior_flowClust2d(fr = parent_data[[1]]
								,xChannel = xChannel
								,yChannel = yChannel
								, K = K)
						
						paired_args[["prior_list"]]<-prior_list
						thisCall[["positive"]]<-NULL #remove positive arg since 2D gate doesn't understand it
					}
									
				}
				#update arg_names
				
				for(arg in names(paired_args))
					thisCall[[arg]]<-paired_args[[arg]]
				
				
				
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
#					browser()
					flist<-eval(thisCall)
				}
				
				#we expect a filter/filters from gm
#				browser()
				if(class(flist[[1]])=="filters")#reconstruct filterlist out of filterslist
				{
					#right now we consider quadgate as the only use case for filters
					gs_node_id<-NULL
#					browser()
					#clock-wise starting from top-left quadrant(assuming X followed by Y channel)											
					quadPatterns<-c(".+-.+\\+$"   #top left  -+
									,".+\\+.+\\+$" #top right ++
									,".+\\+.+-$"  #bottom right +-
									,".+-.+-$")   #bottom left	--									
					for(i in 1:length(popNames))
					{

						curAlias<-popAlias[i]
						curPop<-popNames[i]
						curPopId<-popIds[i]
						#check if popname is give as Y[*]X[*]
						YX_pattern<-paste(dims["yChannel"],".+", dims["xChannel"],".+",sep="")
						XY_pattern<-paste(dims["xChannel"],".+", dims["yChannel"],".+",sep="")
						#do the flipping if YX
						if(grepl(YX_pattern,curPop))
						{
							pos<-regexpr(dims["xChannel"],curPop)
							xterm<-substring(curPop,pos,nchar(curPop))
							yterm<-substring(curPop,1,pos-1)
							toMatch<-paste(xterm,yterm,sep="")
						}else if(grepl(XY_pattern,curPop))
							toMatch<-curPop #keep as it is if XY
						else
							stop("X,Y axis do not match between 'dims'(",dims,") and 'pop'(",curPop,")")
						
						quadInd<-which(unlist(lapply(quadPatterns,grepl,toMatch)))
						#fetch appropriate filter based on the quadrant ind
						curFlist <- lapply(flist, function(curFilters) {
									curFilters[[quadInd]]
								})
#						browser()
						#save the fcFilter if applicable
						if(extends(class(curFlist[[1]]),"fcFilter"))
						{
							fcObjList<-c(fcObjList,new("fcObject"
														,prior=prior_list
														,fcFilters=curFlist
														)
											)
						}
						curFlist <- filterList(curFlist)
						cur_gs_node_id <- add(y, curFlist, parent = parent,name=curAlias)	
						recompute(y, cur_gs_node_id)
						gs_node_id<-c(gs_node_id,cur_gs_node_id)
						
					}
					
					
					
				}else
				{
#				browser()
					#append it to fcObj if it is a fcFilter
					if(extends(class(flist[[1]]),"fcFilter"))
					{
						
						fcObjList[[1]]<-new("fcObject"
											,prior=prior_list
											,fcFilters=flist)
					}	
					
					
					gs_node_id <- add(y, flist, parent = parent,name=popAlias)
					recompute(y, gs_node_id)
					
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
			if(length(fcObjList)>0)
				names(fcObjList)<-popIds
#			gs_node_id
			list(gs_node_id=gs_node_id
				 ,fcObjList=fcObjList
		 		)
			
			
			
		})
		
setMethod("gating", signature = c("boolMethod", "GatingSet")
		, definition = function(x, y,gtPops, parent, ...) 
{
	
	args<-parameters(x)
	gm<-paste(".",names(x),sep="")
	popAlias<-unlist(lapply(gtPops,alias))
	popNames<-unlist(lapply(gtPops,names))
	gs_nodes<-getChildren(y[[1]],getNodes(y[[1]],parent))
	popIds<-unlist(lapply(gtPops,"slot","id"))

	tNodes <-args
	if (!(tNodes %in% gs_nodes)) {
		message(tNodes, " gating...")
		bf <- eval(substitute(booleanFilter(x), list(x = as.symbol(tNodes))))
		bf@filterId <- tNodes
		invisible(gs_node_id <- add(y, bf, parent = parent))
		invisible(recompute(y, gs_node_id))
		message("done.")
	}else
	{
#		browser()
		message("Skip gating!Population '",tNodes,"' already exists.")
		gs_node_id<-getChildren(y[[1]],parent)
		#select the corresponding gs node id by matching the node names
		gs_node_name<-getNodes(y[[1]])[gs_node_id]
		gs_node_id<-gs_node_id[match(tNodes,gs_node_name)]	
	}
	
	
	names(gs_node_id)<-popIds
#	gs_node_id
	list("gs_node_id"=gs_node_id)
})		

setMethod("gating", signature = c("polyFunctions", "GatingSet")
				, definition = function(x, y,gtPops, parent, ...) 
{
	
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
				bgt<-new("boolMethod",name=polyExpr,args=polyExpr)
				gating(bgt,y,parent=parent,gtPops=gtPops,...)
			})					
	
	
	
	message("done.")
		
	
	
	list()
})

## wrappers for the different gating routines
.singletGate <- function(fs, xChannel = "FSC-A", yChannel = "FSC-H", prediction_level = 0.99, ...) {
  require('flowStats')
  # Creates a list of polygon gates based on the prediction bands at the minimum and maximum
  # x_channel observation using a robust linear model trained by flowStats.
  singletGate(fs[[1]], area = xChannel, height = yChannel,
              prediction_level = prediction_level)
}

.flowClust.1d<-function(fs, xChannel = NA,yChannel,tol=1e-3,prior=NULL,filterId=""
						,usePrior="yes",split=TRUE,...)
{
	
	require("flowClust")
	require("MASS")
#	browser()
	
	fr<-fs[[1]]
#		browser()			
	if(is.na(xChannel))#1d gate			
		flowClust.1d(fr = fr
				,params = yChannel
				,tol = tol
				,filterId = filterId
				,prior = prior[[yChannel]]
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
								,prior = prior[[xChannel]]
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
#			browser()
			filter2 <- flowClust.1d(negFr, params = yChannel,
									usePrior = usePrior
									,filterId = as.character(yStain)
									,prior = prior[[yChannel]]
									, ...)
			cut.y.l <- filter2@min
			
			filter3 <- flowClust.1d(posFr, params = yChannel,
									usePrior = usePrior
									,filterId = as.character(yStain)
									,prior = prior[[yChannel]]
									, ...)
			cut.y.r <- filter3@min
		} else {
			filter2 <- flowClust.1d(fr, params = yChannel,
									usePrior = usePrior
									,filterId = as.character(yStain)
									,prior = prior[[yChannel]]
									, ...)

			cut.y.l <- cut.y.r <- filter2@min
		}
		
		###############################################################     
		#construct rectangleGates based on the cuts and popNames,clock-wise
		###############################################################
		gateList <- new("filters")
		
		chnls<-c(xChannel, yChannel)
		markers<-c(xStain,yStain)
#		browser()
		##cut.x,cut.y.l (1st,4th quad)
		postX<-posteriors(filter1)
		postY<-posteriors(filter2)
		postList<-c(postX,postY)
		
		coord <- list(c(-Inf, cut.x), c(cut.y.l, Inf))
		names(coord) <- as.character(chnls)
		thisFcFilter<-fcRectangleGate(rectangleGate(coord),postList)
		gateList[[paste(paste0(markers, c("-", "+")), collapse="")]] <-thisFcFilter 
		
		coord <- list(c(-Inf, cut.x), c(-Inf, cut.y.l))
		names(coord) <- as.character(chnls)
		thisFcFilter<-fcRectangleGate(rectangleGate(coord),postList)
		gateList[[paste(paste0(markers, c("-", "-")), collapse="")]] <- thisFcFilter
		

		##cut.x,cut.y.r(2nd,3rd quad)
		postX<-posteriors(filter1)
		if(split)
			postY<-	posteriors(filter2)
		else
			postY<-posteriors(filter3)
		postList<-c(postX,postY)
		
		coord <- list(c(cut.x, Inf), c(cut.y.r, Inf))
		names(coord) <- as.character(chnls)
		thisFcFilter<-fcRectangleGate(rectangleGate(coord),postList)
		gateList[[paste(paste0(markers, c("+", "+")), collapse="")]] <- thisFcFilter
		
		coord <- list(c(cut.x, Inf), c(-Inf, cut.y.r))
		names(coord) <- as.character(chnls)
		thisFcFilter<-fcRectangleGate(rectangleGate(coord),postList)
		gateList[[paste(paste0(markers, c("+", "-")), collapse="")]] <- thisFcFilter
		
		
		
		gateList
	}
	
}	
	


.mindensity <- function(fs, yChannel = "FSC-A", filterId = "", ...) {
  mindensity(flow_frame = fs[[1]], channel = yChannel, filter_id = filterId, ...)

}




.flowClust.2d <- function(fs, xChannel, yChannel, usePrior = "yes", prior = NULL, ...) {
  require("flowClust")
  fr <- fs[[1]]
  flowClust.2d(fr = fr, xChannel = xChannel, yChannel = yChannel,
               usePrior = usePrior, prior = prior, ...)
}

.rangeGate <- function(fs, xChannel = NA, yChannel, absolute = FALSE,
                       filterId = "", ...) {
  require('flowStats')
  fr <- fs[[1]]
  rangeGate(x = fr, stain = yChannel, inBetween = TRUE, absolute = absolute,
            filterId = filterId, ...)
}

.quantileGate <- function(fs, xChannel = NA, yChannel, probs = 0.999,
                          filterId = "", ...) {
  fr <- fs[[1]]
  quantileGate(fr = fr, probs = probs, stain = yChannel, filterId = filterId, ...)
}

.quadrantGate <- function(fs, xChannel = NA, yChannel, ...) {
  require('flowStats')
  fr<-fs[[1]]
  
  qfilter<-quadrantGate(fr, stain = c(xChannel,yChannel), absolute = FALSE, inBetween = TRUE, ...)
  
  ###############################################################     
  #construct rectangleGates based on the cuts and popNames,clock-wise
  ###############################################################
  cut.x<-qfilter@boundary[xChannel]
  cut.y<-qfilter@boundary[yChannel]
  gateList <- new("filters")
  
  chnls<-c(xChannel, yChannel)
  markers<-chnls
  
  coord <- list(c(-Inf, cut.x), c(cut.y, Inf))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("-", "+")), collapse="")]] <- rectangleGate(coord)
  
  coord <- list(c(cut.x, Inf), c(cut.y, Inf))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("+", "+")), collapse="")]] <- rectangleGate(coord)
  
  coord <- list(c(cut.x, Inf), c(-Inf, cut.y))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("+", "-")), collapse="")]] <- rectangleGate(coord)
  
  coord <- list(c(-Inf, cut.x), c(-Inf, cut.y))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("-", "-")), collapse="")]] <- rectangleGate(coord)
  
  gateList
}

