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

setMethod("gating", signature = c("gatingMethod", "GatingSet")
			, definition = function(x, gs,gtPops, parent
									,num_nodes = 1, parallel_type = c("multicore", "sock")
									,plot = FALSE, xbin = 128, ...) 
{
	
	require('parallel')
	
		
		args<-parameters(x)
		
		
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
		gs_nodes<-getNodes(gs[[1]])
		
		
#		browser()
		if (!any(grepl(popAlias, gs_nodes))) 
		{
			parent_data <- getData(gs, parent)
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
				
			}else #QuadGate
			{
				
			}
#			browser()
			#parse the string into named vector
			paired_args<-strsplit(args,split="\\,")[[1]]
			named_args<-unlist(lapply(paired_args,function(paired_arg){
													curPair<-strsplit(paired_arg,split="=")[[1]]
													curPair_value<-curPair[2]
													names(curPair_value)<-curPair[1]
													curPair_value
												}
										)
								)
			named_args<-as.list(named_args)
			#parse args for flowClust 
			#prior estimation is done separately from flowClust routine 
			#because piror_flowClust1d requires the entire parent flowSet yet flowClust only takes one flowFrame
			#####################################################################
			if(grepl("^\\.flowClust\\.[12]d$",gm))
			{
				
				if(gm==".flowClust.1d")
				{
					#get the value of neg and pos
					neg_cluster<-as.integer(named_args["neg"])			
					K<-as.integer(named_args["pos"])+neg_cluster
					
					
					# Elicitation of priors for flowClust
					prior <- list()
					if(!is.na(xChannel))
						prior$xChannel <- prior_flowClust1d(flow_set = parent_data,channel = xChannel, K = K)
					
					
					prior$yChannel <- prior_flowClust1d(flow_set = parent_data,channel = yChannel, K = K)
					
					#replace neg and pos and convert the named vector back to string
					named_args[["K"]]<-K
					named_args[["neg_cluster"]]<-neg_cluster
					named_args[["prior"]]<-prior
					named_args<-named_args[-match(c("neg","pos"),names(named_args))]	
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
			
			#we expect a filter as returned value from gm
	
			
			# Adds the list of singlet polygon gates to the workflow.
			node_id <- add(gs, flist, parent = parent,name=popAlias)
			recompute(gs, node_id)
			message("done.")
			
		}else
		{
			message("Skip gating!Population '",popAlias,"' already exists.")
#			browser()
			node_id<-getChildren(gs[[1]],parent)
		}
		
		if (plot) {
			print(plotGate(gs, node_id, xbin = xbin, pos = c(0.5, 0.8)))
		}
		node_id		
})

## wrappers for the different gating routines
.singletGate <- function(fs
		, xChannel = "FSC-A"
		,yChannel = "FSC-H"
		, prediction_level = 0.99
		,...) {
	require('flowStats')
	# Creates a list of polygon gates based on the prediction bands at the minimum and maximum
	# x_channel observation using a robust linear model trained by flowStats.
	
	
	singletGate(fs[[1]], area = xChannel
			, height = yChannel
			,prediction_level = prediction_level
	)$gate
	
}
.flowClust.1d<-function(fs
		, xChannel = NA
		,yChannel
		,tol=1e-3
		,prior=NULL
		,filterId=""
		,usePrior="yes"
		,...
		)
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
			gate_x <- flowClust.1d(fr = fr
									, params = xChannel, tol = tol
									,filterId = as.character(getChannelMarker(fr, xChannel)$desc)
									,prior = prior$xChannel
									,usePrior=usePrior
									, ...)
			gate_y <- flowClust.1d(fr = fr
									, params = yChannel
									, tol = tol
									,filterId = as.character(getChannelMarker(fr,yChannel)$desc)
									,prior = prior$yChannel
									,usePrior=usePrior
									, ...)
			
			gate_x <- ifelse(is.finite(gate_x@min), gate_x@min, gate_x@max)
			gate_y <- ifelse(is.finite(gate_y@min), gate_y@min, gate_y@max)
			
			flowClust_quadGate <- list(gate_x, gate_y)
			names(flowClust_quadGate) <- c(xChannel,yChannel)
			quadGate(filterId = filterId, flowClust_quadGate)
			
		}
}
.flowClust.2d<-function(fs
		, xChannel
		,yChannel
#		,filterId=""
		,usePrior="yes"
		,...
){
	
	require("flowClust")
	
	
	fr<-fs[[1]]
#	browser()
	flowClust.2d(fr = fr
				, xChannel = xChannel
				, yChannel = yChannel
				,usePrior=usePrior
				, ...)
	
}
.rangeGate<-function(fs
						, xChannel = NA
						,yChannel
						,absolute = FALSE
						,filterId=""
						,...
				)
{
	require('flowStats')
	rangeGate(fs[[1]], stain = yChannel, inBetween = TRUE, absolute = absolute,
			filterId = filterId, ...)	
}


.quantileGate<-function(fs
						, xChannel = NA
						,yChannel
						,probs = 0.999
						,filterId=""
						,...
						)
{
	quantileGate(fr = fs[[1]], probs = probs, stain = yChannel, filterId = filterId, ...)	
}

.quadrantGate<-function(fs
						, xChannel = NA
						,yChannel
						,...
						)
{
	require('flowStats')
	quadrantGate(fs[[1]], stain = c(xChannel,yChannel), absolute = FALSE, inBetween = TRUE, ...)			
}
